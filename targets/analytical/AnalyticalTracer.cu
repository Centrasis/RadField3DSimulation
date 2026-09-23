#include "analytical/AnalyticalTracer.hpp"
#include <cufft.h>
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RadiationSimulation::Analytical;

#define CUDA_CHECK(call) do { cudaError_t _e = (call); if (_e != cudaSuccess) \
    throw std::runtime_error(std::string("CUDA error ") + cudaGetErrorString(_e) + " at " + __FILE__ + ":" + std::to_string(__LINE__)); } while (0)
#define CUFFT_CHECK(call) do { cufftResult _e = (call); if (_e != CUFFT_SUCCESS) \
    throw std::runtime_error(std::string("cuFFT error ") + std::to_string(_e) + " at " + __FILE__ + ":" + std::to_string(__LINE__)); } while (0)

namespace {

// NIST water mass attenuation (cm^2/g) and incoherent (Compton) interaction fraction, 10-150 keV.
__constant__ float c_e_kev[10]      = {10.f, 15.f, 20.f, 30.f, 40.f, 50.f, 60.f, 80.f, 100.f, 150.f};
__constant__ float c_mu_rho[10]     = {5.329f, 1.673f, 0.8096f, 0.3756f, 0.2683f, 0.2269f, 0.2059f, 0.1837f, 0.1707f, 0.1505f};
__constant__ float c_f_compton[10]  = {0.036f, 0.10f, 0.20f, 0.43f, 0.60f, 0.70f, 0.76f, 0.83f, 0.86f, 0.90f};

__device__ float loglog(float e_kev, const float* xs, const float* ys)
{
    e_kev = fminf(fmaxf(e_kev, xs[0]), xs[9]);
    const float le = logf(e_kev);
    int i = 0;
#pragma unroll
    for (int k = 0; k < 9; ++k)
        if (le >= logf(xs[k])) i = k;
    const float t = (le - logf(xs[i])) / (logf(xs[i + 1]) - logf(xs[i]));
    return expf(logf(ys[i]) + t * (logf(ys[i + 1]) - logf(ys[i])));
}

struct DeviceParams {
    int nx, ny, nz, bins;
    float a, bin_width_ev;
    float ox, oy, oz, dx, dy, dz, e1x, e1y, e1z, e2x, e2y, e2z;
    float distance, rect_w, rect_h, cone_half_rad, focal_spot_m;
    int shape;
    float wx, wy, wz;  // world half-extent
};

// Smooth beam-edge coverage in [0,1] from a finite focal spot: the penumbra half-width at axial
// distance st is (focal_spot · st / D), floored to half a voxel so the edge is always ≥1-voxel
// smooth. `lat` is the lateral coordinate, `halfsize` the geometric edge half-extent at st.
__device__ inline float edge_coverage(float lat, float halfsize, float st, const DeviceParams& p)
{
    const float pen = fmaxf(p.focal_spot_m * st / p.distance, 0.5f * p.a);
    const float x = 0.5f + 0.5f * (halfsize - fabsf(lat)) / pen;   // signed edge distance → ramp
    const float c = fminf(fmaxf(x, 0.f), 1.f);
    return c * c * (3.f - 2.f * c);   // smoothstep
}

// Trilinear density sample at a world-space point (voxel centres sit at (i+0.5)*a - w).
__device__ float sample_density_tri(const float* __restrict__ d, const DeviceParams& p,
                                    float wx, float wy, float wz)
{
    const float gx = (wx + p.wx) / p.a - 0.5f;
    const float gy = (wy + p.wy) / p.a - 0.5f;
    const float gz = (wz + p.wz) / p.a - 0.5f;
    const int i0 = (int)floorf(gx), j0 = (int)floorf(gy), k0 = (int)floorf(gz);
    const float fx = gx - i0, fy = gy - j0, fz = gz - k0;
    float acc = 0.f;
#pragma unroll
    for (int dz = 0; dz < 2; ++dz)
    for (int dy = 0; dy < 2; ++dy)
    for (int dx = 0; dx < 2; ++dx) {
        const int i = i0 + dx, j = j0 + dy, k = k0 + dz;
        if (i < 0 || i >= p.nx || j < 0 || j >= p.ny || k < 0 || k >= p.nz) continue;
        const float w = (dx ? fx : 1.f - fx) * (dy ? fy : 1.f - fy) * (dz ? fz : 1.f - fz);
        acc += w * d[(size_t)k * p.ny * p.nx + (size_t)j * p.nx + i];
    }
    return acc;
}

// Direct beam: fixed-step radiological path from the focal spot to each voxel center, per-bin
// Beer-Lambert attenuation of the (already collimated) primary fluence.
__global__ void k_direct(DeviceParams p, const float* __restrict__ density,
                         const float* __restrict__ spectrum, float* __restrict__ direct_spec)
{
    const int v = blockIdx.x * blockDim.x + threadIdx.x;
    const int n = p.nx * p.ny * p.nz;
    if (v >= n) return;
    const int x = v % p.nx;
    const int y = (v / p.nx) % p.ny;
    const int z = v / (p.nx * p.ny);
    const float cx = (x + 0.5f) * p.a - p.wx;
    const float cy = (y + 0.5f) * p.a - p.wy;
    const float cz = (z + 0.5f) * p.a - p.wz;

    const float rx = cx - p.ox, ry = cy - p.oy, rz = cz - p.oz;
    const float t_ax = rx * p.dx + ry * p.dy + rz * p.dz;   // signed distance along the beam axis
    const float dist = sqrtf(rx * rx + ry * ry + rz * rz);

    // Collimated primary fluence with a smooth focal-spot penumbra at the beam edge (no hard step).
    float fluence = 0.f;
    if (t_ax > 0.f) {
        if (p.shape == 0) {  // rectangle: field grows linearly with axial distance
            const float lat1 = rx * p.e1x + ry * p.e1y + rz * p.e1z;
            const float lat2 = rx * p.e2x + ry * p.e2y + rz * p.e2z;
            const float sc = fmaxf(t_ax / p.distance, 1e-6f);
            const float cov = edge_coverage(lat1, 0.5f * p.rect_w * sc, t_ax, p)
                            * edge_coverage(lat2, 0.5f * p.rect_h * sc, t_ax, p);
            if (cov > 0.f) {
                const float tt = fmaxf(t_ax, p.a);
                fluence = cov * p.a * p.a * p.distance * p.distance / (p.rect_w * p.rect_h * tt * tt);
            }
        } else {  // cone — penumbra on the polar-angle edge
            const float ang = acosf(fminf(fmaxf((rx * p.dx + ry * p.dy + rz * p.dz) / fmaxf(dist, p.a), -1.f), 1.f));
            const float pen_ang = fmaxf(p.focal_spot_m / p.distance, 0.5f * p.a / fmaxf(dist, p.a));
            const float xr = 0.5f + 0.5f * (p.cone_half_rad - ang) / pen_ang;
            const float cc = fminf(fmaxf(xr, 0.f), 1.f);
            const float cov = cc * cc * (3.f - 2.f * cc);
            if (cov > 0.f) {
                const float omega = 2.f * 3.14159265358979f * (1.f - cosf(p.cone_half_rad));
                fluence = cov * p.a * p.a / (omega * fmaxf(dist * dist, p.a * p.a));
            }
        }
    }

    float* out = direct_spec + (size_t)v * p.bins;
    if (fluence <= 0.f) {
        for (int b = 0; b < p.bins; ++b) out[b] = 0.f;
        return;
    }

    // radiological path (g/cm^3 * m), fixed-step midpoint sampling with TRILINEAR density lookup
    // (smooth attenuation gradient — nearest-voxel gives staircase banding).
    const int steps = 128;
    float radio = 0.f;
    const float seg = dist / steps;
    for (int s = 0; s < steps; ++s) {
        const float f = (s + 0.5f) / steps;
        radio += sample_density_tri(density, p, p.ox + f * rx, p.oy + f * ry, p.oz + f * rz) * seg;
    }

    for (int b = 0; b < p.bins; ++b) {
        const float e_kev = (b + 0.5f) * p.bin_width_ev / 1000.f;
        const float tau = loglog(e_kev, c_e_kev, c_mu_rho) * radio * 100.f;
        out[b] = spectrum[b] * expf(-tau) * fluence;
    }
}

// Single-scatter source per voxel: local interaction rate of the direct beam that Compton-scatters,
// redistributed in energy by the precomputed Compton shift matrix M[out,in].
__global__ void k_scatter_source(DeviceParams p, const float* __restrict__ density,
                                 const float* __restrict__ direct_spec, const float* __restrict__ compton_M,
                                 float* __restrict__ q_out)
{
    const int v = blockIdx.x * blockDim.x + threadIdx.x;
    const int n = p.nx * p.ny * p.nz;
    if (v >= n) return;
    const float rho = density[v];
    const float* ds = direct_spec + (size_t)v * p.bins;
    float* q = q_out + (size_t)v * p.bins;
    // interacted-and-Compton contribution per incoming bin
    for (int in = 0; in < p.bins; ++in) {
        const float e_kev = (in + 0.5f) * p.bin_width_ev / 1000.f;
        const float mu_rho = loglog(e_kev, c_e_kev, c_mu_rho);
        const float fc = loglog(e_kev, c_e_kev, c_f_compton);
        const float interact = 1.f - expf(-mu_rho * rho * 100.f * p.a);
        q[in] = ds[in] * interact * fc;   // temporarily hold the pre-shift source
    }
    // energy redistribution: qshift[out] = sum_in M[out,in] * q[in]  (done into a register buffer)
    float shifted[256];
    for (int o = 0; o < p.bins; ++o) {
        float acc = 0.f;
        for (int in = 0; in < p.bins; ++in)
            acc += compton_M[o * p.bins + in] * q[in];
        shifted[o] = acc;
    }
    for (int o = 0; o < p.bins; ++o)
        q[o] = shifted[o];
}

// Padded, periodic 1/(4*pi*r^2) point-spread kernel for the scatter transport convolution.
__global__ void k_build_kernel(int Px, int Py, int Pz, float a, float* __restrict__ kern)
{
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    const long total = (long)Px * Py * Pz;
    if (i >= total) return;
    const int px = i % Px;
    const int py = (i / Px) % Py;
    const int pz = i / (Px * Py);
    const float ox = fminf(px, Px - px) * a;
    const float oy = fminf(py, Py - py) * a;
    const float oz = fminf(pz, Pz - pz) * a;
    const float r2 = fmaxf(ox * ox + oy * oy + oz * oz, a * a);
    kern[i] = a * a / (4.f * 3.14159265358979f * r2);
}

__global__ void k_scatter_into_pad(int nx, int ny, int nz, int Px, int Py, int bin, int bins,
                                   const float* __restrict__ q, float* __restrict__ pad)
{
    const int v = blockIdx.x * blockDim.x + threadIdx.x;
    if (v >= nx * ny * nz) return;
    const int x = v % nx, y = (v / nx) % ny, z = v / (nx * ny);
    pad[(size_t)z * Py * Px + (size_t)y * Px + x] = q[(size_t)v * bins + bin];
}

__global__ void k_pad_to_scatter(int nx, int ny, int nz, int Px, int Py, int bin, int bins,
                                 const float* __restrict__ pad, float norm, float* __restrict__ scatter_spec)
{
    const int v = blockIdx.x * blockDim.x + threadIdx.x;
    if (v >= nx * ny * nz) return;
    const int x = v % nx, y = (v / nx) % ny, z = v / (nx * ny);
    const float val = pad[(size_t)z * Py * Px + (size_t)y * Px + x] * norm;
    scatter_spec[(size_t)v * bins + bin] = fmaxf(val, 0.f);
}

__global__ void k_complex_mul(cufftComplex* __restrict__ a, const cufftComplex* __restrict__ b, long m)
{
    const long i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= m) return;
    const cufftComplex x = a[i], y = b[i];
    a[i].x = x.x * y.x - x.y * y.y;
    a[i].y = x.x * y.y + x.y * y.x;
}

// flux = sum_E spec; spectrum normalized to sum 1 per voxel (0 where flux 0); error ~ 1/sqrt(counts).
__global__ void k_finalize(int n, int bins, float particles, const float* __restrict__ spec,
                          float* __restrict__ flux, float* __restrict__ norm_spec, float* __restrict__ error)
{
    const int v = blockIdx.x * blockDim.x + threadIdx.x;
    if (v >= n) return;
    const float* s = spec + (size_t)v * bins;
    float sum = 0.f;
    for (int b = 0; b < bins; ++b) sum += s[b];
    flux[v] = sum;
    float* o = norm_spec + (size_t)v * bins;
    if (sum > 0.f)
        for (int b = 0; b < bins; ++b) o[b] = s[b] / sum;
    else
        for (int b = 0; b < bins; ++b) o[b] = 0.f;
    error[v] = fminf(1.f, rsqrtf(fmaxf(sum * particles, 1.f)));
}

std::vector<float> build_compton_matrix(int bins, float bin_width_ev)
{
    // M[out,in]: angle-averaged Compton energy redistribution (cos(theta) uniform in [-1,1]).
    const int n_angles = 256;
    std::vector<float> m((size_t)bins * bins, 0.f);
    const float width_kev = bin_width_ev / 1000.f;
    for (int j = 0; j < bins; ++j) {
        const float e = (j + 0.5f) * width_kev;
        for (int s = 0; s < n_angles; ++s) {
            const float cos_t = -1.f + 2.f * (s + 0.5f) / n_angles;
            const float e_out = e / (1.f + (e / 511.f) * (1.f - cos_t));
            int o = (int)std::floor(e_out / width_kev - 0.5f);
            if (o < 0) o = 0;
            if (o >= bins) o = bins - 1;
            m[(size_t)o * bins + j] += 1.f / n_angles;
        }
    }
    return m;
}

}  // namespace

namespace RadiationSimulation {
namespace Analytical {

void run_analytical(const AnalyticalParams& params,
                    const std::vector<float>& density,
                    const std::vector<float>& spectrum,
                    const AnalyticalOutput& out)
{
    const int nx = params.nx, ny = params.ny, nz = params.nz, bins = params.bins;
    const int n = nx * ny * nz;
    if (bins > 256) throw std::runtime_error("analytical tracer supports at most 256 spectrum bins");

    DeviceParams dp{};
    dp.nx = nx; dp.ny = ny; dp.nz = nz; dp.bins = bins; dp.a = params.voxel_m; dp.bin_width_ev = params.bin_width_ev;
    dp.ox = params.origin[0]; dp.oy = params.origin[1]; dp.oz = params.origin[2];
    dp.dx = params.direction[0]; dp.dy = params.direction[1]; dp.dz = params.direction[2];
    dp.e1x = params.e1[0]; dp.e1y = params.e1[1]; dp.e1z = params.e1[2];
    dp.e2x = params.e2[0]; dp.e2y = params.e2[1]; dp.e2z = params.e2[2];
    dp.distance = params.distance; dp.rect_w = params.rect_w; dp.rect_h = params.rect_h;
    dp.cone_half_rad = params.cone_half_rad; dp.shape = params.shape; dp.focal_spot_m = params.focal_spot_m;
    dp.wx = 0.5f * nx * params.voxel_m; dp.wy = 0.5f * ny * params.voxel_m; dp.wz = 0.5f * nz * params.voxel_m;

    float *d_density, *d_spectrum, *d_direct_spec, *d_q, *d_scatter_spec, *d_compton;
    CUDA_CHECK(cudaMalloc(&d_density, (size_t)n * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_spectrum, (size_t)bins * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_direct_spec, (size_t)n * bins * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_q, (size_t)n * bins * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_scatter_spec, (size_t)n * bins * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_compton, (size_t)bins * bins * sizeof(float)));
    CUDA_CHECK(cudaMemcpy(d_density, density.data(), (size_t)n * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_spectrum, spectrum.data(), (size_t)bins * sizeof(float), cudaMemcpyHostToDevice));
    const std::vector<float> compton = build_compton_matrix(bins, params.bin_width_ev);
    CUDA_CHECK(cudaMemcpy(d_compton, compton.data(), (size_t)bins * bins * sizeof(float), cudaMemcpyHostToDevice));

    const int threads = 256;
    const int blocks = (n + threads - 1) / threads;
    k_direct<<<blocks, threads>>>(dp, d_density, d_spectrum, d_direct_spec);
    CUDA_CHECK(cudaGetLastError());
    k_scatter_source<<<blocks, threads>>>(dp, d_density, d_direct_spec, d_compton, d_q);
    CUDA_CHECK(cudaGetLastError());

    // --- scatter transport: batched-per-bin FFT convolution with the 1/(4 pi r^2) kernel ---
    const int Px = 2 * nx, Py = 2 * ny, Pz = 2 * nz;
    const long pad_real = (long)Px * Py * Pz;
    const long pad_cplx = (long)Pz * Py * (Px / 2 + 1);
    float* d_pad; cufftComplex *d_specF, *d_kernF;
    CUDA_CHECK(cudaMalloc(&d_pad, pad_real * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_specF, pad_cplx * sizeof(cufftComplex)));
    CUDA_CHECK(cudaMalloc(&d_kernF, pad_cplx * sizeof(cufftComplex)));
    cufftHandle plan_r2c, plan_c2r;
    CUFFT_CHECK(cufftPlan3d(&plan_r2c, Pz, Py, Px, CUFFT_R2C));  // x is the fastest (last) dim
    CUFFT_CHECK(cufftPlan3d(&plan_c2r, Pz, Py, Px, CUFFT_C2R));

    k_build_kernel<<<(pad_real + threads - 1) / threads, threads>>>(Px, Py, Pz, params.voxel_m, d_pad);
    CUDA_CHECK(cudaGetLastError());
    CUFFT_CHECK(cufftExecR2C(plan_r2c, d_pad, d_kernF));

    const float fft_norm = 1.f / (float)pad_real;
    for (int b = 0; b < bins; ++b) {
        CUDA_CHECK(cudaMemset(d_pad, 0, pad_real * sizeof(float)));
        k_scatter_into_pad<<<blocks, threads>>>(nx, ny, nz, Px, Py, b, bins, d_q, d_pad);
        CUFFT_CHECK(cufftExecR2C(plan_r2c, d_pad, d_specF));
        k_complex_mul<<<(pad_cplx + threads - 1) / threads, threads>>>(d_specF, d_kernF, pad_cplx);
        CUFFT_CHECK(cufftExecC2R(plan_c2r, d_specF, d_pad));
        k_pad_to_scatter<<<blocks, threads>>>(nx, ny, nz, Px, Py, b, bins, d_pad, fft_norm, d_scatter_spec);
    }
    CUDA_CHECK(cudaGetLastError());

    // --- finalize both channels ---
    float *d_dflux, *d_dnorm, *d_derr, *d_sflux, *d_snorm, *d_serr;
    CUDA_CHECK(cudaMalloc(&d_dflux, n * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_dnorm, (size_t)n * bins * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_derr, n * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_sflux, n * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_snorm, (size_t)n * bins * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_serr, n * sizeof(float)));
    k_finalize<<<blocks, threads>>>(n, bins, params.particles, d_direct_spec, d_dflux, d_dnorm, d_derr);
    k_finalize<<<blocks, threads>>>(n, bins, params.particles, d_scatter_spec, d_sflux, d_snorm, d_serr);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // Copy straight into the field-owned destination buffers (x-fastest flat order matches VoxelLayer).
    CUDA_CHECK(cudaMemcpy(out.direct.flux, d_dflux, n * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(out.scatter.flux, d_sflux, n * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(out.direct.error, d_derr, n * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(out.scatter.error, d_serr, n * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(out.direct.spectrum, d_dnorm, (size_t)n * bins * sizeof(float), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(out.scatter.spectrum, d_snorm, (size_t)n * bins * sizeof(float), cudaMemcpyDeviceToHost));

    cufftDestroy(plan_r2c); cufftDestroy(plan_c2r);
    for (void* p : {(void*)d_density, (void*)d_spectrum, (void*)d_direct_spec, (void*)d_q, (void*)d_scatter_spec,
                    (void*)d_compton, (void*)d_pad, (void*)d_specF, (void*)d_kernF, (void*)d_dflux, (void*)d_dnorm,
                    (void*)d_derr, (void*)d_sflux, (void*)d_snorm, (void*)d_serr})
        cudaFree(p);
}

}  // namespace Analytical
}  // namespace RadiationSimulation
