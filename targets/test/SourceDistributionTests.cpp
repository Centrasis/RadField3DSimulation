#include "RadiationSimulation.hpp"
#include "RadiationSource.hpp"
#include <gtest/gtest.h>
#include <Randomize.hh>
#if defined _WIN32 || defined _WIN64
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#include <fstream>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/constants.hpp>

namespace {
	const RadiationSimulation::UniformRandom uniform = [] { return G4UniformRand(); };

	TEST(RectangleShape, Sampling)
	{
		RadiationSimulation::RectangleSourceShape shape(glm::vec2(1.0f, 0.5f), 1.0f);

		for (size_t i = 0; i < 100; i++) {
			glm::vec3 direction = shape.drawRayDirection(uniform);
			EXPECT_GE(direction.x, -0.5f);
			EXPECT_LE(direction.x, 0.5f);
			EXPECT_GE(direction.y, -0.25f);
			EXPECT_LE(direction.y, 0.25f);
		}
	}

	// Slopes (x / -z, y / -z) of a direction: its point on the plane at unit distance.
	glm::dvec2 slopes_of(const glm::vec3& direction)
	{
		return glm::dvec2(direction.x / -direction.z, direction.y / -direction.z);
	}

	TEST(RectangleShape, IsAPointSourceWithCutOutEdges)
	{
		// a wide field (corner at 25°), where the per-area fall-off of a point source is large (cos³ = 0.74 at the corner)
		const glm::vec2 size(0.8f, 0.5f);
		const double half_x = 0.4, half_y = 0.25;
		constexpr size_t bins = 8;
		const size_t n = 1000000;
		auto bin_of = [&](const glm::dvec2& s) {
			const size_t i = std::min(bins - 1, static_cast<size_t>((s.x + half_x) / (2.0 * half_x) * bins));
			const size_t j = std::min(bins - 1, static_cast<size_t>((s.y + half_y) / (2.0 * half_y) * bins));
			return i * bins + j;
		};

		G4Random::setTheSeed(4711);
		RadiationSimulation::RectangleSourceShape rectangle(size, 1.0f);
		std::vector<double> collimated(bins * bins, 0.0);
		for (size_t i = 0; i < n; i++) {
			const glm::dvec2 s = slopes_of(rectangle.drawRayDirection(uniform));
			ASSERT_LE(std::abs(s.x), half_x * (1.0 + 1e-5));
			ASSERT_LE(std::abs(s.y), half_y * (1.0 + 1e-5));
			collimated[bin_of(s)] += 1.0 / n;
		}

		// reference: an isotropic point source (the cone around the corner), with everything outside the rectangle cut away
		RadiationSimulation::ConeSourceShape cone(static_cast<float>(glm::degrees(std::atan(std::hypot(half_x, half_y)))) + 0.01f);
		std::vector<double> point_source(bins * bins, 0.0);
		for (size_t kept = 0; kept < n;) {
			const glm::dvec2 s = slopes_of(cone.drawRayDirection(uniform));
			if (std::abs(s.x) > half_x || std::abs(s.y) > half_y)
				continue;
			point_source[bin_of(s)] += 1.0 / n;
			kept++;
		}

		for (size_t b = 0; b < bins * bins; b++) {
			// five standard deviations of the difference of two binomial shares
			const double tolerance = 5.0 * std::sqrt(2.0 * point_source[b] / n);
			EXPECT_NEAR(collimated[b], point_source[b], tolerance) << "bin " << b / bins << ", " << b % bins;
		}

		// the edges are really darker per area: corner against centre bins, cos³ of the corner bin centre ≈ 0.78
		const double centre = 0.25 * (collimated[(bins / 2 - 1) * bins + bins / 2 - 1] + collimated[(bins / 2 - 1) * bins + bins / 2]
			+ collimated[(bins / 2) * bins + bins / 2 - 1] + collimated[(bins / 2) * bins + bins / 2]);
		const double corner_tan = std::hypot(half_x * (1.0 - 1.0 / bins), half_y * (1.0 - 1.0 / bins));
		const double corner_cos3 = std::pow(1.0 / std::sqrt(1.0 + corner_tan * corner_tan), 3.0);
		const double centre_cos3 = std::pow(1.0 / std::sqrt(1.0 + std::pow(std::hypot(half_x, half_y) / bins, 2.0)), 3.0);
		EXPECT_NEAR(collimated[0] / centre, corner_cos3 / centre_cos3, 0.02);
	}

	TEST(RectangleShape, FillsTheRectangleUpToItsEdges)
	{
		G4Random::setTheSeed(4712);
		RadiationSimulation::RectangleSourceShape shape(glm::vec2(0.3f, 0.02f), 2.0f);
		glm::dvec2 max_slope(0.0);
		for (size_t i = 0; i < 100000; i++) {
			const glm::dvec2 s = slopes_of(shape.drawRayDirection(uniform));
			max_slope = glm::max(max_slope, glm::abs(s));
		}
		EXPECT_GT(max_slope.x, 0.99 * 0.075);
		EXPECT_LE(max_slope.x, 0.075 * (1.0 + 1e-5));
		EXPECT_GT(max_slope.y, 0.99 * 0.005);
		EXPECT_LE(max_slope.y, 0.005 * (1.0 + 1e-5));
		EXPECT_THROW(RadiationSimulation::RectangleSourceShape(glm::vec2(0.3f), 0.f), std::invalid_argument);
	}

	TEST(RectangleShape, OutputSampling)
	{
		glm::vec2 rect_size(0.3f, 0.2f);
		RadiationSimulation::RectangleSourceShape shape(rect_size, 1.0f);
		std::array<std::array<float, 1000>, 500>* image = new std::array<std::array<float, 1000>, 500>();
		glm::uvec2 img_dim(
			image->size(),
			(*image)[0].size()
		);
		glm::uvec2 half_img_dim = glm::uvec2(img_dim.x / 2, img_dim.y / 2);

		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				(*image)[i][j] = 0.0f;
			}
		}

		size_t samples = 1000000;
		for (size_t i = 0; i < samples; i++) {
			glm::vec3 direction = shape.drawRayDirection(uniform);

			glm::uvec2 idx(
				direction.y * static_cast<float>(img_dim.y),
				direction.x * static_cast<float>(img_dim.y)
			);
			idx += half_img_dim;

			if (idx.x >= 0 && idx.x < half_img_dim.x * 2 && idx.y >= 0 && idx.y < half_img_dim.y * 2)
				(*image)[idx.x][idx.y] += 1.0f;
		}

		float max_count = 0.0f;
		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				if ((*image)[i][j] > max_count) {
					max_count = (*image)[i][j];
				}
			}
		}

		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				(*image)[i][j] /= max_count;
			}
		}

		auto file_path = fs::absolute(fs::path("rectangle_sampling.bmp"));
		if (fs::exists(file_path)) {
			fs::remove(file_path);
		}

		// write image array to bitmap file and create if it does not exist
		std::ofstream file(file_path, std::ios::binary);

		// write header
		file << "BM";
		uint32_t fileSize = 54 + 4 * img_dim.x * img_dim.y;
		file.write(reinterpret_cast<char*>(&fileSize), sizeof(uint32_t));
		uint32_t reserved = 0;
		file.write(reinterpret_cast<char*>(&reserved), sizeof(uint32_t));
		uint32_t offset = 54;
		file.write(reinterpret_cast<char*>(&offset), sizeof(uint32_t));
		uint32_t headerSize = 40;
		file.write(reinterpret_cast<char*>(&headerSize), sizeof(uint32_t));
		uint32_t width = img_dim.y;
		file.write(reinterpret_cast<char*>(&width), sizeof(uint32_t));
		uint32_t height = img_dim.x;
		file.write(reinterpret_cast<char*>(&height), sizeof(uint32_t));
		uint16_t planes = 1;
		file.write(reinterpret_cast<char*>(&planes), sizeof(uint16_t));
		uint16_t bitsPerPixel = 32;
		file.write(reinterpret_cast<char*>(&bitsPerPixel), sizeof(uint16_t));
		uint32_t compression = 0;
		file.write(reinterpret_cast<char*>(&compression), sizeof(uint32_t));
		uint32_t imageSize = 4 * img_dim.x * img_dim.y;
		file.write(reinterpret_cast<char*>(&imageSize), sizeof(uint32_t));
		uint32_t xPixelsPerMeter = 0;
		file.write(reinterpret_cast<char*>(&xPixelsPerMeter), sizeof(uint32_t));
		uint32_t yPixelsPerMeter = 0;
		file.write(reinterpret_cast<char*>(&yPixelsPerMeter), sizeof(uint32_t));
		uint32_t colorsUsed = 0;
		file.write(reinterpret_cast<char*>(&colorsUsed), sizeof(uint32_t));
		uint32_t importantColors = 0;
		file.write(reinterpret_cast<char*>(&importantColors), sizeof(uint32_t));

		// write image data
		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				uint32_t color = static_cast<uint32_t>((*image)[i][j] * 255.0f);
				file.write(reinterpret_cast<char*>(&color), sizeof(uint32_t));
			}
		}

		delete image;
	}

	TEST(RectangleShape, OutputSourceSampling)
	{
		glm::vec2 rect_size(0.3f, 0.2f);
		RadiationSimulation::XRaySource source(10.f, std::make_unique<RadiationSimulation::RectangleSourceShape>(rect_size, 2.0f));
		glm::quat rotation = glm::angleAxis(glm::radians(0.f), glm::vec3(0.f, 1.f, 0.f)) * glm::angleAxis(glm::radians(0.f), glm::vec3(1.f, 0.f, 0.f));
		const glm::vec3 source_dir = rotation * glm::vec3(0.f, 0.f, -1.f);
		source.setTransform(glm::vec3(0, 0, -2.f), source_dir);
		std::array<std::array<float, 1000>, 500>* image = new std::array<std::array<float, 1000>, 500>();
		glm::uvec2 img_dim(
			image->size(),
			(*image)[0].size()
		);
		glm::uvec2 half_img_dim = glm::uvec2(img_dim.x / 2, img_dim.y / 2);

		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				(*image)[i][j] = 0.0f;
			}
		}

		size_t samples = 1000000;
		for (size_t i = 0; i < samples; i++) {
			glm::vec3 direction = source.drawRayDirection(uniform);
			glm::vec3 position = direction * 2.f + source.getLocation();

			glm::uvec2 idx(
				position.y * static_cast<float>(img_dim.y),
				position.x * static_cast<float>(img_dim.y)
			);
			idx += half_img_dim;

			if (idx.x >= 0 && idx.x < half_img_dim.x * 2 && idx.y >= 0 && idx.y < half_img_dim.y * 2)
				(*image)[idx.x][idx.y] += 1.0f;
		}

		float max_count = 0.0f;
		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				if ((*image)[i][j] > max_count) {
					max_count = (*image)[i][j];
				}
			}
		}

		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				(*image)[i][j] /= max_count;
			}
		}

		auto file_path = fs::absolute(fs::path("rectangle_source_sampling.bmp"));
		if (fs::exists(file_path)) {
			fs::remove(file_path);
		}

		// write image array to bitmap file and create if it does not exist
		std::ofstream file(file_path, std::ios::binary);

		// write header
		file << "BM";
		uint32_t fileSize = 54 + 4 * img_dim.x * img_dim.y;
		file.write(reinterpret_cast<char*>(&fileSize), sizeof(uint32_t));
		uint32_t reserved = 0;
		file.write(reinterpret_cast<char*>(&reserved), sizeof(uint32_t));
		uint32_t offset = 54;
		file.write(reinterpret_cast<char*>(&offset), sizeof(uint32_t));
		uint32_t headerSize = 40;
		file.write(reinterpret_cast<char*>(&headerSize), sizeof(uint32_t));
		uint32_t width = img_dim.y;
		file.write(reinterpret_cast<char*>(&width), sizeof(uint32_t));
		uint32_t height = img_dim.x;
		file.write(reinterpret_cast<char*>(&height), sizeof(uint32_t));
		uint16_t planes = 1;
		file.write(reinterpret_cast<char*>(&planes), sizeof(uint16_t));
		uint16_t bitsPerPixel = 32;
		file.write(reinterpret_cast<char*>(&bitsPerPixel), sizeof(uint16_t));
		uint32_t compression = 0;
		file.write(reinterpret_cast<char*>(&compression), sizeof(uint32_t));
		uint32_t imageSize = 4 * img_dim.x * img_dim.y;
		file.write(reinterpret_cast<char*>(&imageSize), sizeof(uint32_t));
		uint32_t xPixelsPerMeter = 0;
		file.write(reinterpret_cast<char*>(&xPixelsPerMeter), sizeof(uint32_t));
		uint32_t yPixelsPerMeter = 0;
		file.write(reinterpret_cast<char*>(&yPixelsPerMeter), sizeof(uint32_t));
		uint32_t colorsUsed = 0;
		file.write(reinterpret_cast<char*>(&colorsUsed), sizeof(uint32_t));
		uint32_t importantColors = 0;
		file.write(reinterpret_cast<char*>(&importantColors), sizeof(uint32_t));

		// write image data
		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				uint32_t color = static_cast<uint32_t>((*image)[i][j] * 255.0f);
				file.write(reinterpret_cast<char*>(&color), sizeof(uint32_t));
			}
		}

		delete image;
	}

	TEST(ConeShape, OutputSampling)
	{
		glm::vec2 rect_size(0.3f, 0.2f);
		RadiationSimulation::ConeSourceShape shape(5.0f);
		std::array<std::array<float, 1000>, 500>* image = new std::array<std::array<float, 1000>, 500>();
		glm::uvec2 img_dim(
			image->size(),
			(*image)[0].size()
		);
		glm::uvec2 half_img_dim = glm::uvec2(img_dim.x / 2, img_dim.y / 2);

		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				(*image)[i][j] = 0.0f;
			}
		}

		size_t samples = 1000000;
		for (size_t i = 0; i < samples; i++) {
			glm::vec3 direction = shape.drawRayDirection(uniform);

			glm::uvec2 idx(
				direction.y * static_cast<float>(img_dim.y),
				direction.x * static_cast<float>(img_dim.y)
			);
			idx += half_img_dim;

			if (idx.x >= 0 && idx.x < half_img_dim.x * 2 && idx.y >= 0 && idx.y < half_img_dim.y * 2)
				(*image)[idx.x][idx.y] += 1.0f;
		}

		float max_count = 0.0f;
		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				if ((*image)[i][j] > max_count) {
					max_count = (*image)[i][j];
				}
			}
		}

		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				(*image)[i][j] /= max_count;
			}
		}

		auto file_path = fs::absolute(fs::path("cone_sampling.bmp"));
		if (fs::exists(file_path)) {
			fs::remove(file_path);
		}

		// write image array to bitmap file and create if it does not exist
		std::ofstream file(file_path, std::ios::binary);

		// write header
		file << "BM";
		uint32_t fileSize = 54 + 4 * img_dim.x * img_dim.y;
		file.write(reinterpret_cast<char*>(&fileSize), sizeof(uint32_t));
		uint32_t reserved = 0;
		file.write(reinterpret_cast<char*>(&reserved), sizeof(uint32_t));
		uint32_t offset = 54;
		file.write(reinterpret_cast<char*>(&offset), sizeof(uint32_t));
		uint32_t headerSize = 40;
		file.write(reinterpret_cast<char*>(&headerSize), sizeof(uint32_t));
		uint32_t width = img_dim.y;
		file.write(reinterpret_cast<char*>(&width), sizeof(uint32_t));
		uint32_t height = img_dim.x;
		file.write(reinterpret_cast<char*>(&height), sizeof(uint32_t));
		uint16_t planes = 1;
		file.write(reinterpret_cast<char*>(&planes), sizeof(uint16_t));
		uint16_t bitsPerPixel = 32;
		file.write(reinterpret_cast<char*>(&bitsPerPixel), sizeof(uint16_t));
		uint32_t compression = 0;
		file.write(reinterpret_cast<char*>(&compression), sizeof(uint32_t));
		uint32_t imageSize = 4 * img_dim.x * img_dim.y;
		file.write(reinterpret_cast<char*>(&imageSize), sizeof(uint32_t));
		uint32_t xPixelsPerMeter = 0;
		file.write(reinterpret_cast<char*>(&xPixelsPerMeter), sizeof(uint32_t));
		uint32_t yPixelsPerMeter = 0;
		file.write(reinterpret_cast<char*>(&yPixelsPerMeter), sizeof(uint32_t));
		uint32_t colorsUsed = 0;
		file.write(reinterpret_cast<char*>(&colorsUsed), sizeof(uint32_t));
		uint32_t importantColors = 0;
		file.write(reinterpret_cast<char*>(&importantColors), sizeof(uint32_t));

		// write image data
		for (size_t i = 0; i < img_dim.x; i++) {
			for (size_t j = 0; j < img_dim.y; j++) {
				uint32_t color = static_cast<uint32_t>((*image)[i][j] * 255.0f);
				file.write(reinterpret_cast<char*>(&color), sizeof(uint32_t));
			}
		}

		delete image;
	}

	TEST(EllipsoidShape, StaysInsideAndFillsTheEllipse)
	{
		G4Random::setTheSeed(2024);
		const double tan_x = std::tan(glm::radians(10.0)), tan_y = std::tan(glm::radians(25.0));
		RadiationSimulation::EllipsoidSourceShape shape(glm::vec2(10.f, 25.f));
		double max_x = 0.0, max_y = 0.0;
		for (size_t i = 0; i < 200000; i++) {
			const glm::vec3 direction = shape.drawRayDirection(uniform);
			ASSERT_LT(direction.z, 0.f);
			const double slope_x = direction.x / -direction.z, slope_y = direction.y / -direction.z;
			EXPECT_LE((slope_x / tan_x) * (slope_x / tan_x) + (slope_y / tan_y) * (slope_y / tan_y), 1.0 + 1e-5);
			max_x = std::max(max_x, std::abs(slope_x));
			max_y = std::max(max_y, std::abs(slope_y));
		}
		// the half opening angles are reached, not halved
		EXPECT_GT(max_x, 0.99 * tan_x);
		EXPECT_GT(max_y, 0.99 * tan_y);
	}

	TEST(EllipsoidShape, CircularCaseIsUniformPerSolidAngle)
	{
		G4Random::setTheSeed(2025);
		const double alpha = glm::radians(20.0);
		RadiationSimulation::EllipsoidSourceShape shape(glm::vec2(20.f, 20.f));
		// the cone of polar angles with cos above this value holds half of the solid angle
		const double half_solid_angle_cos = 0.5 * (1.0 + std::cos(alpha));
		const size_t n = 200000;
		size_t inner = 0;
		for (size_t i = 0; i < n; i++)
			inner += -shape.drawRayDirection(uniform).z > half_solid_angle_cos;
		EXPECT_NEAR(static_cast<double>(inner) / n, 0.5, 0.005);
	}

	TEST(EllipsoidShape, ReportsDegreesAndRejectsInvalidAngles)
	{
		RadiationSimulation::EllipsoidSourceShape shape(glm::vec2(10.f, 25.f));
		EXPECT_EQ(shape.getOpeningAnglesDegrees(), glm::vec2(10.f, 25.f));
		EXPECT_THROW(RadiationSimulation::EllipsoidSourceShape(glm::vec2(0.f, 10.f)), std::invalid_argument);
		EXPECT_THROW(RadiationSimulation::EllipsoidSourceShape(glm::vec2(10.f, 90.f)), std::invalid_argument);
	}
}

TEST(CArmRotation, TurnsTheBasePoseOntoTheBeam) {
	RadiationSimulation::XRaySource source(60e3f, std::make_unique<RadiationSimulation::RectangleSourceShape>(glm::vec2(0.2f), 1.f));
	// base pose: tube below the isocentre, beam along +Y — no turn
	source.setTransform(glm::vec3(0.f, -0.6f, 0.f), glm::vec3(0.f, 1.f, 0.f));
	const glm::quat base = source.getCArmRotation();
	EXPECT_NEAR(glm::angle(glm::normalize(base)), 0.f, 1e-5f);
	// any pose (the angles RadField3D samples): +Y goes onto the beam, and the field's own axes turn the same way
	for (float phi : { 0.f, 30.f, 135.f, 270.f }) {
		for (float theta : { 0.f, 45.f, 90.f, 150.f, 180.f }) {
			const glm::quat angles = glm::angleAxis(glm::radians(theta), glm::vec3(1.f, 0.f, 0.f)) * glm::angleAxis(glm::radians(phi), glm::vec3(0.f, 1.f, 0.f));
			const glm::vec3 dir = angles * glm::vec3(0.f, 0.f, -1.f);
			source.setTransform(-dir * 0.6f, dir);
			const glm::quat c = source.getCArmRotation();
			ASSERT_TRUE(std::isfinite(c.w) && std::isfinite(c.x) && std::isfinite(c.y) && std::isfinite(c.z)) << phi << " " << theta;
			EXPECT_NEAR(glm::dot(c * glm::vec3(0.f, 1.f, 0.f), dir), 1.f, 1e-5f) << phi << " " << theta;
		}
	}
}

TEST(SourceShapes, HalfTangentsGiveTheBeamHalfSizeAtUnitDistance) {
	const glm::vec2 rect = RadiationSimulation::RectangleSourceShape(glm::vec2(0.2f, 0.1f), 0.785f).getHalfTangents();
	EXPECT_NEAR(rect.x, 0.1f / 0.785f, 1e-6f);
	EXPECT_NEAR(rect.y, 0.05f / 0.785f, 1e-6f);
	const glm::vec2 cone = RadiationSimulation::ConeSourceShape(10.f).getHalfTangents();
	EXPECT_NEAR(cone.x, std::tan(glm::radians(10.f)), 1e-6f);
	EXPECT_EQ(cone.x, cone.y);
	EXPECT_TRUE(std::isinf(RadiationSimulation::ConeSourceShape(90.f).getHalfTangents().x));
	const glm::vec2 ellipse = RadiationSimulation::EllipsoidSourceShape(glm::vec2(10.f, 20.f)).getHalfTangents();
	EXPECT_NEAR(ellipse.x, std::tan(glm::radians(10.f)), 1e-6f);
	EXPECT_NEAR(ellipse.y, std::tan(glm::radians(20.f)), 1e-6f);
}
