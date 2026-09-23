#pragma once
#include <cstdint>
#include <vector>

namespace RadiationSimulation {
    namespace Analytical {

        /**
         * @brief Beam / grid configuration for one analytical field, in world (field-centered) metres.
         *
         * The grid is x-fastest flat-indexed (idx = z*ny*nx + y*nx + x), matching RadFiled3D's
         * VoxelGrid::get_voxel_idx, so host buffers map directly onto the stored field.
         */
        struct AnalyticalParams {
            int nx = 0, ny = 0, nz = 0;   ///< voxel counts per axis
            float voxel_m = 0.f;          ///< isotropic voxel size (m)
            int bins = 0;                 ///< field-spectrum histogram bins
            float bin_width_ev = 0.f;     ///< field-spectrum bin width (eV)
            float origin[3] = {0, 0, 0};  ///< tube focal spot, world-centered (m)
            float direction[3] = {0, 0, -1}; ///< beam axis (unit, points toward the field center)
            float e1[3] = {1, 0, 0};      ///< in-plane collimation axis (beam frame x)
            float e2[3] = {0, 1, 0};      ///< in-plane collimation axis (beam frame y)
            float distance = 1.f;         ///< source distance (m); rectangle field size is defined here
            int shape = 0;                ///< 0 = rectangle, 1 = cone
            float rect_w = 0.f, rect_h = 0.f; ///< rectangle field size at the isocenter plane (m)
            float cone_half_rad = 0.f;    ///< cone half opening angle (rad)
            float focal_spot_m = 0.f;     ///< effective focal-spot size (m) → beam-edge penumbra blur
            float particles = 1.f;        ///< primaries (only scales the mimicked statistical error)
        };

        /**
         * @brief Destination buffers for one channel — pointers INTO the memory a CartesianRadiationField
         * already allocated for its layers (VoxelLayer data buffers are contiguous in x-fastest flat
         * order, matching the tracer). The GPU results are copied straight into them; no intermediate
         * host allocation. Spectra are voxel-major, bin-minor, normalized per voxel to sum 1.
         */
        struct AnalyticalChannelBuffers {
            float* flux = nullptr;      ///< n = nx*ny*nz
            float* error = nullptr;     ///< n
            float* spectrum = nullptr;  ///< n * bins
        };

        struct AnalyticalOutput {
            AnalyticalChannelBuffers scatter;
            AnalyticalChannelBuffers direct;
        };

        /**
         * @brief Run the analytical (non-Monte-Carlo) direct-beam + single-scatter estimate on the GPU,
         * writing the result directly into the caller's (field-owned) channel buffers.
         *
         * @param params     Beam/grid configuration.
         * @param density    World density grid in g/cm^3, length nx*ny*nz, x-fastest flat order.
         * @param spectrum   Tube spectrum rebinned to `bins`, normalized to sum 1.
         * @param out        Destination flux/error/spectrum buffers for the scatter and direct channels.
         */
        void run_analytical(const AnalyticalParams& params,
                            const std::vector<float>& density,
                            const std::vector<float>& spectrum,
                            const AnalyticalOutput& out);
    }
}
