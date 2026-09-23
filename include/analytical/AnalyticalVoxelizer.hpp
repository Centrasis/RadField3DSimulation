#pragma once
#include <memory>
#include <string>
#include <vector>
#include <glm/vec3.hpp>
#include "Geometry.hpp"

namespace RadiationSimulation {
    namespace Analytical {

        /**
         * @brief Rasterizes loaded meshes into a material-density grid (g/cm^3) for the analytical
         * attenuation model.
         *
         * Placement matches the Geant4 path (Geant4::Mesh / Geant4::SceneConstructor): scale is baked into the
         * vertices, rotation is applied about the mesh origin, translation composes ADDITIVELY through
         * parents (a child follows its parent's world position). Interior filling uses a per-column
         * z-scanline parity test (even-odd rule) — O(triangles), not O(voxels * triangles).
         *
         * The grid is x-fastest flat-indexed (idx = z*ny*nx + y*nx + x); the world origin is the grid
         * center. Density is looked up from the mesh MaterialName via a small tissue table (soft tissue
         * default, lung, bone, water, air); overlapping meshes: the last one placed wins per voxel.
         */
        class AnalyticalVoxelizer {
        public:
            /**
             * @param root_meshes Meshes as returned by GeometryLoader::Load (with desc transforms).
             * @param counts      Voxel counts (nx, ny, nz).
             * @param voxel_m     Isotropic voxel size in metres.
             * @return Density grid (g/cm^3), length nx*ny*nz, x-fastest flat order.
             */
            static std::vector<float> voxelize_density(
                const std::vector<std::shared_ptr<Geometry::Mesh>>& root_meshes,
                const glm::ivec3& counts,
                float voxel_m);

            /** @brief Density (g/cm^3) for a Geant4 material name; 1.0 if unknown. */
            static float material_density(const std::string& material_name);
        };
    }
}
