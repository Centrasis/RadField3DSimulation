#pragma once
#include <cstddef>
#include <cstdint>
#include <glm/vec3.hpp>
#include <glm/mat4x4.hpp>
#include "Geometry.hpp"


namespace RadiationSimulation::Geometry::Voxelization {
	/** Axis-aligned grid of cubic voxels in world coordinates (metres).
	* Voxels are flat-indexed x-fastest (idx = (z * ny + y) * nx + x), the same order as RadFiled3D's VoxelGrid.
	*/
	struct VoxelGrid {
		glm::uvec3 voxel_counts;
		double voxel_size;
		/// World position of the lower corner of voxel (0, 0, 0).
		glm::dvec3 min_corner;

		inline size_t voxel_count() const {
			return static_cast<size_t>(this->voxel_counts.x) * this->voxel_counts.y * this->voxel_counts.z;
		}
	};

	/** Read-only view on the closed surface of a solid, made of triangles and quads, in world coordinates (metres).
	* Implementations are queried concurrently from several threads.
	*/
	class ISurface {
	public:
		virtual ~ISurface() = default;
		virtual size_t polygon_count() const = 0;
		/// 3 for a triangle, 4 for a quad.
		virtual size_t vertex_count(size_t polygon) const = 0;
		virtual glm::dvec3 vertex(size_t polygon, size_t vertex) const = 0;
	};

	/** Surface of a Geometry::Mesh, whose vertices are mapped to world coordinates (metres) by `local_to_world`. */
	class TransformedMeshSurface : public ISurface {
	protected:
		const Mesh& mesh;
		const glm::dmat4 local_to_world;
	public:
		TransformedMeshSurface(const Mesh& mesh, const glm::dmat4& local_to_world);
		virtual size_t polygon_count() const override;
		virtual size_t vertex_count(size_t polygon) const override;
		virtual glm::dvec3 vertex(size_t polygon, size_t vertex) const override;
	};

	/** Marks every voxel of `grid` that overlaps the solid enclosed by `surface` with 255 in `occupancy`
	* (grid.voxel_count() entries). All other entries stay untouched, so the solids of one kind can be accumulated
	* into the same mask.
	*
	* A voxel overlaps the solid when a part of it with non-zero volume lies inside the solid; voxels the surface only
	* touches are not marked. This is decided on an OpenVDB signed distance field of the solid: voxels whose centre lies
	* inside or closer than half a voxel to the surface overlap, voxels further away than half a voxel diagonal do not,
	* and the voxels in between are tested exactly against the surface polygon closest to their centre. Solids thinner
	* than a voxel are therefore still marked. The surface must be closed up to gaps smaller than a voxel; the inside of
	* a surface with larger holes is not defined. Solids reaching beyond the grid are clipped to it.
	*
	* Memory is bounded by the narrow band of the one solid being processed and freed before returning.
	*
	* @param surface The closed surface of the solid in world coordinates.
	* @param grid The voxel grid to mark.
	* @param occupancy Mask of grid.voxel_count() bytes in the grid's flat voxel order.
	* @param max_threads Maximum number of worker threads. -1 uses all available cores.
	* @return The number of voxels that overlap the solid.
	*/
	size_t mark_overlapping_voxels(const ISurface& surface, const VoxelGrid& grid, uint8_t* occupancy, int max_threads = -1);
}
