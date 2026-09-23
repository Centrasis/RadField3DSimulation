#include "Voxelization.hpp"
#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <tbb/global_control.h>
#include <glm/glm.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>


using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Geometry::Voxelization;


TransformedMeshSurface::TransformedMeshSurface(const Mesh& mesh, const glm::dmat4& local_to_world)
	: mesh(mesh),
	  local_to_world(local_to_world)
{
}

size_t TransformedMeshSurface::polygon_count() const
{
	return this->mesh.faceCount();
}

size_t TransformedMeshSurface::vertex_count(size_t polygon) const
{
	return this->mesh.getFaces()[polygon]->getType() == FaceType::Quad ? 4 : 3;
}

glm::dvec3 TransformedMeshSurface::vertex(size_t polygon, size_t vertex) const
{
	const Face* face = this->mesh.getFaces()[polygon];
	const size_t index = (face->getType() == FaceType::Quad)
		? static_cast<const QuadFace*>(face)->getIndices()[static_cast<glm::length_t>(vertex)]
		: static_cast<const TriFace*>(face)->getIndices()[static_cast<glm::length_t>(vertex)];
	return glm::dvec3(this->local_to_world * glm::dvec4(glm::dvec3(this->mesh.getVertices()[index]), 1.0));
}

namespace {
	// Band widths (in voxels) of the signed distance field. The exterior band must exceed half a voxel diagonal
	// (~0.87 voxels) so every voxel the surface can reach carries a distance and a closest polygon.
	constexpr float EXTERIOR_BAND_VOXELS = 2.f;
	constexpr float INTERIOR_BAND_VOXELS = 2.f;
	// Relative tolerance (in voxels) around the certain-overlap distance, absorbing the float precision of the
	// distance field. Voxels within it are decided by the exact polygon test.
	constexpr double DISTANCE_TOLERANCE_VOXELS = 1e-3;
	// The exact test shrinks the voxel by this fraction of its size, so polygons only touching a voxel face do not
	// count as overlap.
	constexpr double TOUCH_TOLERANCE_VOXELS = 1e-6;

	// OpenVDB MeshDataAdapter over an ISurface, reporting vertices in the index space of the distance field.
	class IndexSpaceSurface {
		const ISurface& surface;
		const openvdb::math::Transform& transform;
	public:
		IndexSpaceSurface(const ISurface& surface, const openvdb::math::Transform& transform)
			: surface(surface), transform(transform) {}

		size_t polygonCount() const { return this->surface.polygon_count(); }
		size_t pointCount() const { return this->surface.polygon_count() * 4; }
		size_t vertexCount(size_t n) const { return this->surface.vertex_count(n); }
		void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const {
			const glm::dvec3 world = this->surface.vertex(n, v);
			pos = this->transform.worldToIndex(openvdb::Vec3d(world.x, world.y, world.z));
		}
	};

	// Separating axis test of a triangle against an axis-aligned box centred at the origin with the given half
	// extents (Akenine-Moeller). Touching counts as overlap; callers shrink the box to exclude it.
	bool triangle_overlaps_box(const std::array<glm::dvec3, 3>& triangle, const glm::dvec3& half_extents)
	{
		const glm::dvec3& v0 = triangle[0];
		const glm::dvec3& v1 = triangle[1];
		const glm::dvec3& v2 = triangle[2];
		const std::array<glm::dvec3, 3> edges = { v1 - v0, v2 - v1, v0 - v2 };

		for (const glm::dvec3& edge : edges) {
			for (int axis = 0; axis < 3; axis++) {
				glm::dvec3 unit(0.0);
				unit[axis] = 1.0;
				const glm::dvec3 separating_axis = glm::cross(unit, edge);
				const double p0 = glm::dot(v0, separating_axis);
				const double p1 = glm::dot(v1, separating_axis);
				const double p2 = glm::dot(v2, separating_axis);
				const double radius = half_extents.x * std::abs(separating_axis.x) + half_extents.y * std::abs(separating_axis.y) + half_extents.z * std::abs(separating_axis.z);
				if (std::min({ p0, p1, p2 }) > radius || std::max({ p0, p1, p2 }) < -radius)
					return false;
			}
		}

		for (int axis = 0; axis < 3; axis++) {
			if (std::min({ v0[axis], v1[axis], v2[axis] }) > half_extents[axis] || std::max({ v0[axis], v1[axis], v2[axis] }) < -half_extents[axis])
				return false;
		}

		const glm::dvec3 normal = glm::cross(edges[0], edges[1]);
		const double radius = half_extents.x * std::abs(normal.x) + half_extents.y * std::abs(normal.y) + half_extents.z * std::abs(normal.z);
		return std::abs(glm::dot(normal, v0)) <= radius;
	}

	bool polygon_overlaps_voxel(const ISurface& surface, size_t polygon, const glm::dvec3& voxel_center, const glm::dvec3& half_extents)
	{
		const glm::dvec3 a = surface.vertex(polygon, 0) - voxel_center;
		const glm::dvec3 b = surface.vertex(polygon, 1) - voxel_center;
		const glm::dvec3 c = surface.vertex(polygon, 2) - voxel_center;
		if (triangle_overlaps_box({ a, b, c }, half_extents))
			return true;
		if (surface.vertex_count(polygon) == 4) {
			const glm::dvec3 d = surface.vertex(polygon, 3) - voxel_center;
			return triangle_overlaps_box({ a, c, d }, half_extents);
		}
		return false;
	}
}

size_t RadiationSimulation::Geometry::Voxelization::mark_overlapping_voxels(const ISurface& surface, const VoxelGrid& grid, uint8_t* occupancy, int max_threads)
{
	if (grid.voxel_size <= 0.0)
		throw std::invalid_argument("The voxel size must be positive.");
	if (surface.polygon_count() == 0 || grid.voxel_count() == 0)
		return 0;

	// Skip solids that cannot reach the grid before paying for a distance field.
	glm::dvec3 surface_min(std::numeric_limits<double>::max());
	glm::dvec3 surface_max(std::numeric_limits<double>::lowest());
	for (size_t polygon = 0; polygon < surface.polygon_count(); polygon++) {
		for (size_t v = 0; v < surface.vertex_count(polygon); v++) {
			const glm::dvec3 p = surface.vertex(polygon, v);
			surface_min = glm::min(surface_min, p);
			surface_max = glm::max(surface_max, p);
		}
	}
	const glm::dvec3 grid_max = grid.min_corner + glm::dvec3(grid.voxel_counts) * grid.voxel_size;
	if (glm::any(glm::greaterThan(surface_min, grid_max)) || glm::any(glm::lessThan(surface_max, grid.min_corner)))
		return 0;

	openvdb::initialize();
	std::unique_ptr<tbb::global_control> thread_limit;
	if (max_threads > 0)
		thread_limit = std::make_unique<tbb::global_control>(tbb::global_control::max_allowed_parallelism, static_cast<size_t>(max_threads));

	// Index-space coordinate (i, j, k) is the centre of grid voxel (i, j, k).
	openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(grid.voxel_size);
	const glm::dvec3 first_center = grid.min_corner + 0.5 * grid.voxel_size;
	transform->postTranslate(openvdb::Vec3d(first_center.x, first_center.y, first_center.z));

	openvdb::Int32Grid closest_polygon;
	const IndexSpaceSurface mesh_adapter(surface, *transform);
	openvdb::FloatGrid::Ptr distance = openvdb::tools::meshToVolume<openvdb::FloatGrid>(
		mesh_adapter, *transform, EXTERIOR_BAND_VOXELS, INTERIOR_BAND_VOXELS, 0, &closest_polygon
	);

	// The narrow band encloses the whole solid: its interior beyond the band is stored as inside tiles.
	const openvdb::CoordBBox band = distance->evalActiveVoxelBoundingBox();
	if (band.empty())
		return 0;
	const openvdb::Coord first(std::max(band.min().x(), 0), std::max(band.min().y(), 0), std::max(band.min().z(), 0));
	const openvdb::Coord last(
		std::min(band.max().x(), static_cast<int>(grid.voxel_counts.x) - 1),
		std::min(band.max().y(), static_cast<int>(grid.voxel_counts.y) - 1),
		std::min(band.max().z(), static_cast<int>(grid.voxel_counts.z) - 1)
	);

	const double certain_overlap_distance = (0.5 - DISTANCE_TOLERANCE_VOXELS) * grid.voxel_size;
	const double possible_overlap_distance = (0.5 * std::sqrt(3.0) + DISTANCE_TOLERANCE_VOXELS) * grid.voxel_size;
	const glm::dvec3 half_extents(0.5 * grid.voxel_size - TOUCH_TOLERANCE_VOXELS * grid.voxel_size);

	openvdb::FloatGrid::ConstAccessor distance_accessor = distance->getConstAccessor();
	openvdb::Int32Grid::ConstAccessor polygon_accessor = closest_polygon.getConstAccessor();
	size_t overlapping = 0;
	for (int z = first.z(); z <= last.z(); z++) {
		for (int y = first.y(); y <= last.y(); y++) {
			for (int x = first.x(); x <= last.x(); x++) {
				const openvdb::Coord ijk(x, y, z);
				const double d = distance_accessor.getValue(ijk);
				bool overlaps = d < certain_overlap_distance;
				if (!overlaps && d <= possible_overlap_distance) {
					int32_t polygon = 0;
					if (polygon_accessor.probeValue(ijk, polygon) && polygon >= 0) {
						const glm::dvec3 center = grid.min_corner + (glm::dvec3(x, y, z) + 0.5) * grid.voxel_size;
						overlaps = polygon_overlaps_voxel(surface, static_cast<size_t>(polygon), center, half_extents);
					}
				}
				if (overlaps) {
					occupancy[(static_cast<size_t>(z) * grid.voxel_counts.y + y) * grid.voxel_counts.x + x] = 255;
					overlapping++;
				}
			}
		}
	}
	return overlapping;
}
