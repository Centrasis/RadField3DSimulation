#include "RadiationSimulation.hpp"
#include "RadiationSource.hpp"
#include <gtest/gtest.h>
#if defined _WIN32 || defined _WIN64
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#include <fstream>
#include <memory>
#include <glm/gtc/quaternion.hpp>

namespace {
	TEST(RectangleShape, Sampling)
	{
		RadiationSimulation::RectangleSourceShape shape(glm::vec2(1.0f, 0.5f), 1.0f);

		for (size_t i = 0; i < 100; i++) {
			glm::vec3 direction = shape.drawRayDirection();
			EXPECT_GE(direction.x, -0.5f);
			EXPECT_LE(direction.x, 0.5f);
			EXPECT_GE(direction.y, -0.25f);
			EXPECT_LE(direction.y, 0.25f);
		}
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
			glm::vec3 direction = shape.drawRayDirection();

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
			glm::vec3 direction = source.drawRayDirection();
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
			glm::vec3 direction = shape.drawRayDirection();

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
}