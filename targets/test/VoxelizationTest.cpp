#include "Voxelization.hpp"
#include "Geometry.hpp"
#include <gtest/gtest.h>
#include <glm/gtc/matrix_transform.hpp>
#include <array>
#include <cmath>
#include <memory>
#include <vector>


using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Geometry::Voxelization;

namespace {
	// Closed axis-aligned box [min, max] made of six outward facing quads.
	std::shared_ptr<Mesh> make_box(const glm::vec3& min, const glm::vec3& max)
	{
		const std::vector<glm::vec3> vertices = {
			{ min.x, min.y, min.z }, { max.x, min.y, min.z }, { max.x, max.y, min.z }, { min.x, max.y, min.z },
			{ min.x, min.y, max.z }, { max.x, min.y, max.z }, { max.x, max.y, max.z }, { min.x, max.y, max.z }
		};
		const std::vector<Face*> faces = {
			new QuadFace(glm::uvec4(0, 3, 2, 1)), new QuadFace(glm::uvec4(4, 5, 6, 7)),
			new QuadFace(glm::uvec4(0, 1, 5, 4)), new QuadFace(glm::uvec4(2, 3, 7, 6)),
			new QuadFace(glm::uvec4(1, 2, 6, 5)), new QuadFace(glm::uvec4(0, 4, 7, 3))
		};
		return std::make_shared<Mesh>(vertices, faces, "box");
	}

	// 10^3 voxels of 0.1 m starting at the origin: voxel i spans [0.1 i, 0.1 (i + 1)] on every axis.
	VoxelGrid make_grid()
	{
		return VoxelGrid{ glm::uvec3(10), 0.1, glm::dvec3(0.0) };
	}

	size_t voxel_index(const VoxelGrid& grid, int x, int y, int z)
	{
		return (static_cast<size_t>(z) * grid.voxel_counts.y + y) * grid.voxel_counts.x + x;
	}

	std::vector<uint8_t> voxelize(const Mesh& mesh, const glm::dmat4& transform, const VoxelGrid& grid)
	{
		std::vector<uint8_t> occupancy(grid.voxel_count(), 0);
		mark_overlapping_voxels(TransformedMeshSurface(mesh, transform), grid, occupancy.data());
		return occupancy;
	}

	// Exact reference: does the box [-half, half] placed by `rotation` and `center` overlap the voxel with non-zero volume?
	bool oriented_box_overlaps_voxel(const glm::dmat3& rotation, const glm::dvec3& center, const glm::dvec3& half, const glm::dvec3& voxel_center, double voxel_half)
	{
		const std::array<glm::dvec3, 3> a = { glm::dvec3(1, 0, 0), glm::dvec3(0, 1, 0), glm::dvec3(0, 0, 1) };
		const std::array<glm::dvec3, 3> b = { rotation[0], rotation[1], rotation[2] };
		const glm::dvec3 offset = center - voxel_center;
		auto separated = [&](const glm::dvec3& axis) {
			if (glm::dot(axis, axis) < 1e-12)
				return false;
			double r_voxel = 0.0, r_box = 0.0;
			for (int i = 0; i < 3; i++) {
				r_voxel += voxel_half * std::abs(glm::dot(a[i], axis));
				r_box += half[i] * std::abs(glm::dot(b[i], axis));
			}
			// touching (equality) is not an overlap with non-zero volume
			return std::abs(glm::dot(offset, axis)) >= r_voxel + r_box - 1e-12;
		};
		for (int i = 0; i < 3; i++) {
			if (separated(a[i]) || separated(b[i]))
				return false;
			for (int j = 0; j < 3; j++)
				if (separated(glm::cross(a[i], b[j])))
					return false;
		}
		return true;
	}
}

TEST(Voxelization, VoxelAlignedBoxMarksExactlyItsVoxels) {
	const VoxelGrid grid = make_grid();
	const auto occupancy = voxelize(*make_box(glm::vec3(0.2f), glm::vec3(0.6f)), glm::dmat4(1.0), grid);
	for (int z = 0; z < 10; z++)
		for (int y = 0; y < 10; y++)
			for (int x = 0; x < 10; x++) {
				const bool inside = x >= 2 && x < 6 && y >= 2 && y < 6 && z >= 2 && z < 6;
				EXPECT_EQ(occupancy[voxel_index(grid, x, y, z)], inside ? 255 : 0) << x << " " << y << " " << z;
			}
}

TEST(Voxelization, ShiftedBoxMarksEveryTouchedVoxelVolume) {
	const VoxelGrid grid = make_grid();
	// spans voxels 2..5 partially on every axis
	const auto occupancy = voxelize(*make_box(glm::vec3(0.225f), glm::vec3(0.525f)), glm::dmat4(1.0), grid);
	size_t marked = 0;
	for (uint8_t v : occupancy)
		marked += v == 255;
	EXPECT_EQ(marked, 4u * 4u * 4u);
	EXPECT_EQ(occupancy[voxel_index(grid, 2, 2, 2)], 255);
	EXPECT_EQ(occupancy[voxel_index(grid, 5, 5, 5)], 255);
	EXPECT_EQ(occupancy[voxel_index(grid, 6, 5, 5)], 0);
}

TEST(Voxelization, PlateThinnerThanAVoxelIsMarked) {
	const VoxelGrid grid = make_grid();
	// 1 mm plate inside voxel layer z = 4, covering the whole grid in x and y
	const auto occupancy = voxelize(*make_box(glm::vec3(-0.5f, -0.5f, 0.445f), glm::vec3(1.5f, 1.5f, 0.446f)), glm::dmat4(1.0), grid);
	for (int z = 0; z < 10; z++)
		for (int y = 0; y < 10; y++)
			for (int x = 0; x < 10; x++)
				EXPECT_EQ(occupancy[voxel_index(grid, x, y, z)], z == 4 ? 255 : 0) << x << " " << y << " " << z;
}

TEST(Voxelization, PlateOnAVoxelBoundaryMarksBothLayers) {
	const VoxelGrid grid = make_grid();
	// 2 cm plate centred on the boundary between voxel layers z = 4 and z = 5
	const auto occupancy = voxelize(*make_box(glm::vec3(-0.5f, -0.5f, 0.49f), glm::vec3(1.5f, 1.5f, 0.51f)), glm::dmat4(1.0), grid);
	size_t missing = 0;
	for (int z = 0; z < 10; z++)
		for (int y = 0; y < 10; y++)
			for (int x = 0; x < 10; x++) {
				const bool expected = z == 4 || z == 5;
				EXPECT_EQ(occupancy[voxel_index(grid, x, y, z)], expected ? 255 : 0) << x << " " << y << " " << z;
				missing += expected && occupancy[voxel_index(grid, x, y, z)] != 255;
			}
	std::cout << "boundary plate: " << missing << " of 200 voxels missing" << std::endl;
}

TEST(Voxelization, TiltedThinPlateLeavesNoColumnEmpty) {
	const VoxelGrid grid = make_grid();
	const glm::dmat4 tilt = glm::translate(glm::dmat4(1.0), glm::dvec3(0.5)) * glm::rotate(glm::dmat4(1.0), glm::radians(30.0), glm::dvec3(1.0, 0.0, 0.0));
	const auto occupancy = voxelize(*make_box(glm::vec3(-2.f, -0.3f, -0.0005f), glm::vec3(2.f, 0.3f, 0.0005f)), tilt, grid);
	for (int x = 0; x < 10; x++) {
		for (int y = 2; y < 8; y++) {
			bool column_marked = false;
			for (int z = 0; z < 10; z++)
				column_marked |= occupancy[voxel_index(grid, x, y, z)] == 255;
			EXPECT_TRUE(column_marked) << x << " " << y;
		}
	}
}

TEST(Voxelization, SolidBeyondTheGridIsClipped) {
	const VoxelGrid grid = make_grid();
	const auto occupancy = voxelize(*make_box(glm::vec3(-0.3f), glm::vec3(0.3f)), glm::dmat4(1.0), grid);
	size_t marked = 0;
	for (uint8_t v : occupancy)
		marked += v == 255;
	EXPECT_EQ(marked, 3u * 3u * 3u);

	const auto outside = voxelize(*make_box(glm::vec3(2.f), glm::vec3(3.f)), glm::dmat4(1.0), grid);
	for (uint8_t v : outside)
		EXPECT_EQ(v, 0);
}

TEST(Voxelization, SolidsAccumulateIntoOneMask) {
	const VoxelGrid grid = make_grid();
	std::vector<uint8_t> occupancy(grid.voxel_count(), 0);
	const auto outer = make_box(glm::vec3(0.1f), glm::vec3(0.4f));
	const auto disjoint = make_box(glm::vec3(0.7f), glm::vec3(0.9f));
	const auto nested = make_box(glm::vec3(0.2f), glm::vec3(0.3f));
	EXPECT_EQ(mark_overlapping_voxels(TransformedMeshSurface(*outer, glm::dmat4(1.0)), grid, occupancy.data()), 27u);
	EXPECT_EQ(mark_overlapping_voxels(TransformedMeshSurface(*disjoint, glm::dmat4(1.0)), grid, occupancy.data()), 8u);
	EXPECT_EQ(mark_overlapping_voxels(TransformedMeshSurface(*nested, glm::dmat4(1.0)), grid, occupancy.data()), 1u);
	size_t marked = 0;
	for (uint8_t v : occupancy)
		marked += v == 255;
	EXPECT_EQ(marked, 27u + 8u);
}

TEST(Voxelization, RotatedBoxMatchesExactOverlap) {
	const VoxelGrid grid = make_grid();
	const glm::dvec3 center(0.47, 0.52, 0.49);
	const glm::dvec3 half(0.27, 0.17, 0.11);
	const glm::dmat4 rotation = glm::rotate(glm::rotate(glm::dmat4(1.0), glm::radians(35.0), glm::dvec3(0, 0, 1)), glm::radians(20.0), glm::dvec3(1, 0, 0));
	const glm::dmat4 placement = glm::translate(glm::dmat4(1.0), center) * rotation;
	const auto occupancy = voxelize(*make_box(-glm::vec3(half), glm::vec3(half)), placement, grid);

	size_t expected = 0, false_positive = 0, false_negative = 0;
	for (int z = 0; z < 10; z++)
		for (int y = 0; y < 10; y++)
			for (int x = 0; x < 10; x++) {
				const glm::dvec3 voxel_center = grid.min_corner + (glm::dvec3(x, y, z) + 0.5) * grid.voxel_size;
				const bool overlaps = oriented_box_overlaps_voxel(glm::dmat3(rotation), center, half, voxel_center, 0.5 * grid.voxel_size);
				const bool marked = occupancy[voxel_index(grid, x, y, z)] == 255;
				expected += overlaps;
				false_positive += marked && !overlaps;
				false_negative += !marked && overlaps;
			}
	EXPECT_GT(expected, 0u);
	EXPECT_EQ(false_positive, 0u);
	// voxels only cut by a polygon that is not the closest one to their centre may be missed
	EXPECT_LE(false_negative, expected / 50);
	std::cout << "rotated box: " << expected << " overlapping voxels, " << false_negative << " missed" << std::endl;
}
