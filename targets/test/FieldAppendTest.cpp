#include "FieldAppend.hpp"
#include <gtest/gtest.h>
#include <radfiled3d/storage/radiation_field_store.hpp>
#include <radfiled3d/voxel_grid.hpp>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <string>


using namespace RadiationSimulation;

namespace {
	constexpr size_t BINS = 4;

	struct RunValues {
		float flux;
		float spectrum[BINS];
		float error;
		float angular;
	};

	// A field laid out like the ones RadField3D writes; voxel 0 gets `hit`, all other voxels `rest`.
	std::shared_ptr<radfiled3d::CartesianRadiationField> make_field(const RunValues& hit, const RunValues& rest, bool with_geometry = false)
	{
		auto field = std::make_shared<radfiled3d::CartesianRadiationField>(glm::vec3(0.2f), glm::vec3(0.1f));
		for (const char* name : { "direct_beam", "scatter_field" }) {
			auto channel = std::static_pointer_cast<radfiled3d::VoxelGridBuffer>(field->add_channel(name));
			channel->add_layer<float>("error", 1.f, "Variance");
			channel->add_layer<float>("flux", 0.f, "counts / primary_particles");
			channel->add_custom_layer<radfiled3d::HistogramVoxel<float>>("spectrum", radfiled3d::HistogramVoxel<float>(BINS, 1000.f, nullptr), 0.f, "eV");
			channel->add_custom_layer<radfiled3d::AngularResolvedVoxel<float>>("angular_flux", radfiled3d::AngularResolvedVoxel<float>(glm::uvec2(2, 1), nullptr), 0.f, "counts / primary_particles");
			for (size_t i = 0; i < channel->get_voxel_count(); i++) {
				const RunValues& v = i == 0 ? hit : rest;
				channel->get_layer<float>("flux")[i] = v.flux;
				channel->get_layer<float>("error")[i] = v.error;
				for (size_t b = 0; b < BINS; b++)
					channel->get_layer<float>("spectrum")[i * BINS + b] = v.spectrum[b];
				channel->get_layer<float>("angular_flux")[i * 2] = v.angular;
				channel->get_layer<float>("angular_flux")[i * 2 + 1] = v.angular;
			}
			channel->set_statistical_error("spectrum", 0.2f);
		}
		if (with_geometry) {
			auto geometry = std::static_pointer_cast<radfiled3d::VoxelGridBuffer>(field->add_channel("geometry"));
			geometry->add_layer<uint8_t>("patient", 0, "occupancy");
			geometry->get_layer<uint8_t>("patient")[3] = 255;
		}
		return field;
	}

	std::shared_ptr<radfiled3d::storage::v1::RadiationFieldMetadata> make_metadata(uint64_t primaries, uint64_t duration, const std::string& geometry = "phantom.obj")
	{
		auto metadata = std::make_shared<radfiled3d::storage::v1::RadiationFieldMetadata>(
			radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation(
				primaries, geometry, "physics",
				radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation::XRayTube(glm::vec3(0.f, 0.f, -1.f), glm::vec3(0.f, 0.f, 1.f), 60000.f, "spectrum.csv")
			),
			radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Software("RadField3D", "1.0", "repo", "commit")
		);
		metadata->add_dynamic_metadata<uint64_t>("simulation_duration_s", duration);
		return metadata;
	}

	std::string read_all(const std::string& file)
	{
		std::ifstream stream(file, std::ios::binary);
		return std::string((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
	}

	std::string test_file(const std::string& name)
	{
		const std::string path = (std::filesystem::temp_directory_path() / name).string();
		std::remove(path.c_str());
		return path;
	}
}

TEST(FieldAppend, CombinesRunsWeightedByPrimariesAndFlux) {
	const std::string file = test_file("rf3_append_combine.rf3");
	// run A (300 primaries) hit voxel 0; run B (100 primaries) hit voxel 0 and, unlike A, voxel 1..7
	const RunValues a_hit{ 2.f, { 1.f, 0.f, 0.f, 0.f }, 0.1f, 4.f };
	const RunValues a_rest{ 0.f, { 0.f, 0.f, 0.f, 0.f }, 1.f, 0.f };
	const RunValues b_hit{ 6.f, { 0.f, 0.5f, 0.5f, 0.f }, 0.2f, 8.f };
	const RunValues b_rest{ 1.f, { 0.f, 0.f, 0.f, 1.f }, 0.3f, 1.f };
	append_radiation_field(make_field(a_hit, a_rest, true), make_metadata(300, 10), file);
	append_radiation_field(make_field(b_hit, b_rest), make_metadata(100, 5), file);

	auto field = std::static_pointer_cast<radfiled3d::CartesianRadiationField>(radfiled3d::storage::FieldStore::load(file));
	auto metadata = std::dynamic_pointer_cast<radfiled3d::storage::v1::RadiationFieldMetadata>(radfiled3d::storage::FieldStore::load_metadata(file));
	EXPECT_EQ(static_cast<uint64_t>(metadata->get_header().simulation.primary_particle_count), 400u);
	EXPECT_EQ(metadata->get_dynamic_metadata<radfiled3d::ScalarVoxel<uint64_t>>("simulation_duration_s").get_data(), 15u);

	const double r = 0.25;
	for (const char* name : { "direct_beam", "scatter_field" }) {
		auto channel = field->get_channel(name);
		const float* flux = channel->get_layer<float>("flux");
		const float* spectrum = channel->get_layer<float>("spectrum");
		const float* error = channel->get_layer<float>("error");
		const float* angular = channel->get_layer<float>("angular_flux");

		// voxel 0: flux-weighted spectrum mix, weights 0.75 * 2 = 1.5 and 0.25 * 6 = 1.5
		EXPECT_FLOAT_EQ(flux[0], 0.75f * 2.f + 0.25f * 6.f);
		EXPECT_FLOAT_EQ(spectrum[0], 0.5f);
		EXPECT_FLOAT_EQ(spectrum[1], 0.25f);
		EXPECT_FLOAT_EQ(spectrum[2], 0.25f);
		EXPECT_FLOAT_EQ(spectrum[3], 0.f);
		EXPECT_FLOAT_EQ(error[0], static_cast<float>(std::sqrt(std::pow(1.5 * 0.1, 2) + std::pow(1.5 * 0.2, 2)) / 3.0));
		EXPECT_FLOAT_EQ(angular[0], 0.75f * 4.f + 0.25f * 8.f);

		// voxel 1: only run B hit it, its spectrum stays normalized
		EXPECT_FLOAT_EQ(flux[1], static_cast<float>(r * 1.0));
		EXPECT_FLOAT_EQ(spectrum[BINS + 3], 1.f);
		EXPECT_FLOAT_EQ(spectrum[BINS + 0], 0.f);
		EXPECT_FLOAT_EQ(error[1], 0.3f);

		EXPECT_FLOAT_EQ(channel->get_statistical_error("spectrum"), static_cast<float>(std::sqrt(std::pow(0.75 * 0.2, 2) + std::pow(0.25 * 0.2, 2))));
	}
	// the geometry of the file is kept although the appended field has none
	ASSERT_TRUE(field->has_channel("geometry"));
	EXPECT_EQ(field->get_channel("geometry")->get_layer<uint8_t>("patient")[3], 255);
	std::remove(file.c_str());
}

TEST(FieldAppend, HugeValuesAreCombinedInDoublePrecision) {
	const std::string file = test_file("rf3_append_huge.rf3");
	const float huge = std::numeric_limits<float>::max() * 0.9f;
	const RunValues v{ huge, { 1.f, 0.f, 0.f, 0.f }, 0.1f, huge };
	append_radiation_field(make_field(v, v), make_metadata(100, 1), file);
	append_radiation_field(make_field(v, v), make_metadata(100, 1), file);

	auto field = std::static_pointer_cast<radfiled3d::CartesianRadiationField>(radfiled3d::storage::FieldStore::load(file));
	EXPECT_FLOAT_EQ(field->get_channel("scatter_field")->get_layer<float>("flux")[0], huge);
	EXPECT_FLOAT_EQ(field->get_channel("scatter_field")->get_layer<float>("angular_flux")[0], huge);
	EXPECT_FLOAT_EQ(field->get_channel("scatter_field")->get_layer<float>("spectrum")[0], 1.f);
	std::remove(file.c_str());
}

TEST(FieldAppend, RefusesWhatItCannotCombineAndLeavesTheFileUnchanged) {
	const std::string file = test_file("rf3_append_refuse.rf3");
	const RunValues v{ 1.f, { 1.f, 0.f, 0.f, 0.f }, 0.1f, 1.f };
	append_radiation_field(make_field(v, v, true), make_metadata(100, 1), file);
	const std::string before = read_all(file);

	auto unknown_layer = make_field(v, v);
	std::static_pointer_cast<radfiled3d::VoxelGridBuffer>(unknown_layer->get_channel("scatter_field"))->add_layer<float>("kerma", 0.f, "Gy");
	EXPECT_THROW(append_radiation_field(unknown_layer, make_metadata(100, 1), file), std::runtime_error);
	EXPECT_EQ(read_all(file), before);

	EXPECT_THROW(append_radiation_field(make_field(v, v), make_metadata(100, 1, "other.obj"), file), std::runtime_error);
	EXPECT_EQ(read_all(file), before);

	auto joined = std::make_shared<radfiled3d::CartesianRadiationField>(glm::vec3(0.2f), glm::vec3(0.1f));
	std::static_pointer_cast<radfiled3d::VoxelGridBuffer>(joined->add_channel("radiation"))->add_layer<float>("flux", 1.f, "counts / primary_particles");
	EXPECT_THROW(append_radiation_field(joined, make_metadata(100, 1), file), std::runtime_error);
	EXPECT_EQ(read_all(file), before);
	std::remove(file.c_str());
}

TEST(FieldAppend, RefusesTheFilesSeedAndRecordsZeroForCombinedRuns) {
	const std::string file = test_file("rf3_append_seed.rf3");
	const RunValues v{ 1.f, { 1.f, 0.f, 0.f, 0.f }, 0.1f, 1.f };
	auto with_seed = [](uint64_t seed) {
		auto metadata = make_metadata(100, 1);
		metadata->add_dynamic_metadata<uint64_t>("random_seed", seed);
		return metadata;
	};
	append_radiation_field(make_field(v, v), with_seed(42), file);
	const std::string before = read_all(file);

	EXPECT_THROW(append_radiation_field(make_field(v, v), with_seed(42), file), std::runtime_error);
	EXPECT_EQ(read_all(file), before);

	append_radiation_field(make_field(v, v), with_seed(43), file);
	auto metadata = std::dynamic_pointer_cast<radfiled3d::storage::v1::RadiationFieldMetadata>(radfiled3d::storage::FieldStore::load_metadata(file));
	EXPECT_EQ(metadata->get_dynamic_metadata<radfiled3d::ScalarVoxel<uint64_t>>("random_seed").get_data(), 0u);
	EXPECT_EQ(static_cast<uint64_t>(metadata->get_header().simulation.primary_particle_count), 200u);
	std::remove(file.c_str());
}

TEST(FieldAppend, AddsTheTubeSpectraOfBothRuns) {
	const std::string file = test_file("rf3_append_tube.rf3");
	const RunValues v{ 1.f, { 1.f, 0.f, 0.f, 0.f }, 0.1f, 1.f };
	auto with_tube = [](std::vector<float> counts, float bin_width) {
		auto metadata = make_metadata(100, 1);
		metadata->set_dynamic_custom_metadata<radfiled3d::HistogramVoxel<float>>("tube_spectrum", radfiled3d::HistogramVoxel<float>(counts.size(), bin_width, nullptr));
		auto histogram = metadata->get_dynamic_metadata<radfiled3d::HistogramVoxel<float>>("tube_spectrum").get_histogram();
		std::copy(counts.begin(), counts.end(), histogram.begin());
		return metadata;
	};
	append_radiation_field(make_field(v, v), with_tube({ 16777216.f, 3.f, 0.f }, 1000.f), file);
	append_radiation_field(make_field(v, v), with_tube({ 16777216.f, 2.f, 7.f }, 1000.f), file);

	auto metadata = std::dynamic_pointer_cast<radfiled3d::storage::v1::RadiationFieldMetadata>(radfiled3d::storage::FieldStore::load_metadata(file));
	auto tube = metadata->get_dynamic_metadata<radfiled3d::HistogramVoxel<float>>("tube_spectrum").get_histogram();
	EXPECT_FLOAT_EQ(tube[0], 33554432.f);
	EXPECT_FLOAT_EQ(tube[1], 5.f);
	EXPECT_FLOAT_EQ(tube[2], 7.f);

	const std::string before = read_all(file);
	EXPECT_THROW(append_radiation_field(make_field(v, v), with_tube({ 1.f, 1.f, 1.f }, 500.f), file), std::runtime_error);
	EXPECT_EQ(read_all(file), before);
	std::remove(file.c_str());
}
