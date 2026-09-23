#include "analytical/AnalyticalVoxelizer.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <queue>
#include <glm/gtc/quaternion.hpp>

using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Analytical;


float AnalyticalVoxelizer::material_density(const std::string& name)
{
    // (prefix, g/cm^3); first matching prefix wins. Mirrors the analytical Python reference.
    static const std::array<std::pair<const char*, float>, 6> table = {{
        {"G4_TISSUE_SOFT", 1.00f},
        {"G4_LUNG",        0.26f},
        {"G4_BONE",        1.85f},
        {"G4_A-150_TISSUE",1.13f},
        {"G4_WATER",       1.00f},
        {"G4_AIR",         0.0012f},
    }};
    for (const auto& [prefix, rho] : table)
        if (name.rfind(prefix, 0) == 0)
            return rho;
    return 1.0f;
}

namespace {

    struct Solid {
        std::vector<glm::vec3> tris;  // 3 vertices per triangle, world space
        float density;
    };

    // Depth-first: parents before children, so a child (e.g. lung inside patient) overrides the
    // parent's density when composited in order.
    void collect_solids(const std::shared_ptr<Mesh>& mesh, const glm::vec3& parent_pos, std::vector<Solid>& out)
    {
        const glm::vec3 world_pos = parent_pos + mesh->getPosition();  // additive through parents (G4 nesting)
        const glm::mat3 rot = glm::mat3_cast(mesh->getRotation());     // about the mesh origin
        const glm::vec3 scale = mesh->getScale();

        Solid solid;
        solid.density = AnalyticalVoxelizer::material_density(mesh->getMaterialName());
        auto to_world = [&](const glm::vec3& v) { return world_pos + rot * (v * scale); };

        const std::vector<glm::vec3>& verts = mesh->getVertices();
        for (const Face* f : mesh->getFaces()) {
            if (f->getType() == FaceType::Tri) {
                const glm::uvec3 i = static_cast<const TriFace*>(f)->getIndices();
                solid.tris.push_back(to_world(verts[i.x]));
                solid.tris.push_back(to_world(verts[i.y]));
                solid.tris.push_back(to_world(verts[i.z]));
            } else {
                const glm::uvec4 i = static_cast<const QuadFace*>(f)->getIndices();
                const glm::vec3 v0 = to_world(verts[i.x]), v1 = to_world(verts[i.y]);
                const glm::vec3 v2 = to_world(verts[i.z]), v3 = to_world(verts[i.w]);
                solid.tris.insert(solid.tris.end(), {v0, v1, v2, v0, v2, v3});
            }
        }
        out.push_back(std::move(solid));
        for (const auto& child : mesh->getChildren())
            collect_solids(child, world_pos, out);
    }

}  // namespace

std::vector<float> AnalyticalVoxelizer::voxelize_density(
    const std::vector<std::shared_ptr<Mesh>>& root_meshes,
    const glm::ivec3& counts, float voxel_m)
{
    // Surface-voxelize + outside flood fill, per solid. Robust to NON-WATERTIGHT meshes
    // (the RAF phantom is open triangle soup — an even-odd parity fill would leave interior
    // holes), and to solids clipped by the world box (each solid gets its own padded local
    // grid, so the flood fill cannot leak into the interior through the world boundary).
    const int nx = counts.x, ny = counts.y, nz = counts.z;
    std::vector<float> density(static_cast<size_t>(nx) * ny * nz, 0.f);
    const glm::vec3 world_lo = -0.5f * glm::vec3(nx, ny, nz) * voxel_m;

    std::vector<Solid> solids;
    for (const auto& mesh : root_meshes)
        collect_solids(mesh, glm::vec3(0.f), solids);

    for (const Solid& solid : solids) {
        if (solid.tris.empty())
            continue;
        glm::vec3 mn(1e30f), mx(-1e30f);
        for (const glm::vec3& v : solid.tris) {
            mn = glm::min(mn, v);
            mx = glm::max(mx, v);
        }
        // local grid, aligned to the world grid, padded by 1 voxel of guaranteed-outside border
        const glm::ivec3 lo(
            (int)std::floor((mn.x - world_lo.x) / voxel_m) - 1,
            (int)std::floor((mn.y - world_lo.y) / voxel_m) - 1,
            (int)std::floor((mn.z - world_lo.z) / voxel_m) - 1);
        const glm::ivec3 hi(
            (int)std::floor((mx.x - world_lo.x) / voxel_m) + 1,
            (int)std::floor((mx.y - world_lo.y) / voxel_m) + 1,
            (int)std::floor((mx.z - world_lo.z) / voxel_m) + 1);
        const int lx = hi.x - lo.x + 1, ly = hi.y - lo.y + 1, lz = hi.z - lo.z + 1;
        if ((long long)lx * ly * lz > (long long)512 * 512 * 512)
            continue;  // degenerate transform guard
        auto lidx = [&](int x, int y, int z) { return ((size_t)z * ly + y) * lx + x; };

        // 0 = unknown, 1 = surface, 2 = outside
        std::vector<uint8_t> state((size_t)lx * ly * lz, 0);

        // mark surface voxels: sample each triangle at sub-voxel spacing
        for (size_t t = 0; t + 2 < solid.tris.size(); t += 3) {
            const glm::vec3 &A = solid.tris[t], &B = solid.tris[t + 1], &C = solid.tris[t + 2];
            const float max_edge = std::max({glm::length(B - A), glm::length(C - A), glm::length(C - B)});
            const int n = std::max(1, (int)std::ceil(max_edge / (0.33f * voxel_m)));
            for (int i = 0; i <= n; ++i) {
                for (int j = 0; j <= n - i; ++j) {
                    const glm::vec3 pnt = A + (B - A) * (float(i) / n) + (C - A) * (float(j) / n);
                    const int x = (int)std::floor((pnt.x - world_lo.x) / voxel_m) - lo.x;
                    const int y = (int)std::floor((pnt.y - world_lo.y) / voxel_m) - lo.y;
                    const int z = (int)std::floor((pnt.z - world_lo.z) / voxel_m) - lo.z;
                    if (x >= 0 && x < lx && y >= 0 && y < ly && z >= 0 && z < lz)
                        state[lidx(x, y, z)] = 1;
                }
            }
        }

        // flood fill "outside" from the padded border across non-surface voxels
        std::queue<glm::ivec3> queue;
        auto push_outside = [&](int x, int y, int z) {
            const size_t id = lidx(x, y, z);
            if (state[id] == 0) {
                state[id] = 2;
                queue.push({x, y, z});
            }
        };
        for (int y = 0; y < ly; ++y)
            for (int x = 0; x < lx; ++x) { push_outside(x, y, 0); push_outside(x, y, lz - 1); }
        for (int z = 0; z < lz; ++z)
            for (int x = 0; x < lx; ++x) { push_outside(x, 0, z); push_outside(x, ly - 1, z); }
        for (int z = 0; z < lz; ++z)
            for (int y = 0; y < ly; ++y) { push_outside(0, y, z); push_outside(lx - 1, y, z); }
        while (!queue.empty()) {
            const glm::ivec3 c = queue.front(); queue.pop();
            if (c.x > 0) push_outside(c.x - 1, c.y, c.z);
            if (c.x < lx - 1) push_outside(c.x + 1, c.y, c.z);
            if (c.y > 0) push_outside(c.x, c.y - 1, c.z);
            if (c.y < ly - 1) push_outside(c.x, c.y + 1, c.z);
            if (c.z > 0) push_outside(c.x, c.y, c.z - 1);
            if (c.z < lz - 1) push_outside(c.x, c.y, c.z + 1);
        }

        // solid = surface + enclosed interior; write the part inside the world grid
        for (int z = 0; z < lz; ++z) {
            const int wz = z + lo.z;
            if (wz < 0 || wz >= nz) continue;
            for (int y = 0; y < ly; ++y) {
                const int wy = y + lo.y;
                if (wy < 0 || wy >= ny) continue;
                for (int x = 0; x < lx; ++x) {
                    const int wx = x + lo.x;
                    if (wx < 0 || wx >= nx) continue;
                    if (state[lidx(x, y, z)] != 2)
                        density[((size_t)wz * ny + wy) * nx + wx] = solid.density;
                }
            }
        }
    }

    return density;
}
