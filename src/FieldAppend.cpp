#include "FieldAppend.hpp"
#include <radfiled3d/storage/radiation_field_store.hpp>
#include <radfiled3d/helpers/file_lock.hpp>
#include <radfiled3d/voxel.hpp>
#include <radfiled3d/voxel_grid.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <limits>
#include <atomic>
#include <chrono>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#if defined _WIN32 || defined _WIN64
#include <process.h>
#else
#include <unistd.h>
#endif


namespace fs = std::filesystem;

namespace {
	uint64_t current_process_id()
	{
#if defined _WIN32 || defined _WIN64
		return static_cast<uint64_t>(_getpid());
#else
		return static_cast<uint64_t>(getpid());
#endif
	}

	// The layers RadField3D writes; any other layer has no known combination rule.
	void require_known_layer(const std::string& channel, const std::string& layer)
	{
		if (layer != "flux" && layer != "error" && layer != "spectrum" && layer != "angular_flux")
			throw std::runtime_error("Cannot append layer '" + layer + "' of channel '" + channel + "': its meaning is unknown, so it cannot be combined without risking wrong data.");
	}

	float to_stored(double value, const std::string& what)
	{
		if (!std::isfinite(value) || std::abs(value) > static_cast<double>(std::numeric_limits<float>::max()))
			throw std::runtime_error("Appending yields " + std::to_string(value) + " for " + what + ", which is not representable as float32. The file was left unchanged.");
		return static_cast<float>(value);
	}

	void require(bool condition, const std::string& message)
	{
		if (!condition)
			throw std::runtime_error("Cannot append to the field file: " + message);
	}

	void check_same_simulation(const radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation& existing, const radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation& added)
	{
		require(std::string(existing.geometry) == std::string(added.geometry), "the geometry differs ('" + std::string(existing.geometry) + "' vs. '" + std::string(added.geometry) + "').");
		require(std::string(existing.physics_list) == std::string(added.physics_list), "the physics list differs.");
		require(existing.tube.max_energy_eV == added.tube.max_energy_eV, "the tube's maximum energy differs.");
		require(existing.tube.radiation_direction == added.tube.radiation_direction, "the radiation direction differs.");
		require(existing.tube.radiation_origin == added.tube.radiation_origin, "the radiation origin differs.");
		require(std::string(existing.tube.tube_id) == std::string(added.tube.tube_id), "the tube spectrum differs.");
	}

	// Combines one radiation channel of the new run (`added`, share `r` of all primaries) into `existing` in place.
	void combine_channel(const std::string& name, radfiled3d::VoxelBuffer& existing, radfiled3d::VoxelBuffer& added, double r)
	{
		const size_t voxels = existing.get_voxel_count();
		require(voxels == added.get_voxel_count(), "channel '" + name + "' has a different voxel count.");
		const std::vector<std::string> layers = existing.get_layers();
		require(layers == added.get_layers(), "channel '" + name + "' has different layers.");
		require(existing.has_layer("flux"), "channel '" + name + "' has no flux layer to weight its spectra with.");
		for (const std::string& layer : layers) {
			const radfiled3d::IVoxel& existing_voxel = existing.get_voxel_flat(layer, 0);
			const radfiled3d::IVoxel& added_voxel = added.get_voxel_flat(layer, 0);
			require(existing_voxel.get_type() == added_voxel.get_type() && existing_voxel.get_bytes() == added_voxel.get_bytes(), "layer '" + layer + "' of channel '" + name + "' has a different type or size.");
			require(existing.get_layer_unit(layer) == added.get_layer_unit(layer), "layer '" + layer + "' of channel '" + name + "' has a different unit.");
			require_known_layer(name, layer);
		}
		if (existing.has_layer("spectrum")) {
			const auto& existing_hist = existing.get_voxel_flat<radfiled3d::HistogramVoxel<float>>("spectrum", 0);
			const auto& added_hist = added.get_voxel_flat<radfiled3d::HistogramVoxel<float>>("spectrum", 0);
			require(existing_hist.get_bins() == added_hist.get_bins() && existing_hist.get_histogram_bin_width() == added_hist.get_histogram_bin_width(), "the spectra of channel '" + name + "' have different bins.");
		}

		const double q = 1.0 - r;
		float* flux_a = existing.get_layer<float>("flux");
		const float* flux_b = added.get_layer<float>("flux");

		// Per-voxel weights of both runs' contributions, taken before the flux is overwritten.
		std::vector<double> weight_a(voxels), weight_b(voxels);
		for (size_t i = 0; i < voxels; i++) {
			weight_a[i] = q * static_cast<double>(flux_a[i]);
			weight_b[i] = r * static_cast<double>(flux_b[i]);
		}

		if (existing.has_layer("spectrum")) {
			const size_t bins = existing.get_voxel_flat<radfiled3d::HistogramVoxel<float>>("spectrum", 0).get_bins();
			float* spec_a = existing.get_layer<float>("spectrum");
			const float* spec_b = added.get_layer<float>("spectrum");
			for (size_t i = 0; i < voxels; i++) {
				const double w = weight_a[i] + weight_b[i];
				double sum = 0.0;
				std::vector<double> mixed(bins, 0.0);
				if (w > 0.0) {
					for (size_t b = 0; b < bins; b++) {
						mixed[b] = (weight_a[i] * spec_a[i * bins + b] + weight_b[i] * spec_b[i * bins + b]) / w;
						sum += mixed[b];
					}
				}
				for (size_t b = 0; b < bins; b++)
					spec_a[i * bins + b] = sum > 0.0 ? to_stored(mixed[b] / sum, "the spectrum of channel '" + name + "'") : 0.f;
			}
		}

		if (existing.has_layer("error")) {
			float* error_a = existing.get_layer<float>("error");
			const float* error_b = added.get_layer<float>("error");
			for (size_t i = 0; i < voxels; i++) {
				const double w = weight_a[i] + weight_b[i];
				if (w > 0.0) {
					const double ea = weight_a[i] * error_a[i];
					const double eb = weight_b[i] * error_b[i];
					error_a[i] = to_stored(std::sqrt(ea * ea + eb * eb) / w, "the error of channel '" + name + "'");
				}
				else {
					error_a[i] = std::min(error_a[i], error_b[i]);
				}
			}
		}

		for (size_t i = 0; i < voxels; i++)
			flux_a[i] = to_stored(q * flux_a[i] + r * static_cast<double>(flux_b[i]), "the flux of channel '" + name + "'");

		if (existing.has_layer("angular_flux")) {
			const size_t values = existing.get_voxel_flat("angular_flux", 0).get_bytes() / sizeof(float) * voxels;
			float* angular_a = existing.get_layer<float>("angular_flux");
			const float* angular_b = added.get_layer<float>("angular_flux");
			for (size_t i = 0; i < values; i++)
				angular_a[i] = to_stored(q * angular_a[i] + r * static_cast<double>(angular_b[i]), "the angular flux of channel '" + name + "'");
		}

		for (const std::string& layer : layers) {
			const double ea = q * existing.get_statistical_error(layer);
			const double eb = r * added.get_statistical_error(layer);
			existing.set_statistical_error(layer, to_stored(std::sqrt(ea * ea + eb * eb), "the statistical error of layer '" + layer + "'"));
		}
	}
}

void RadiationSimulation::append_radiation_field(std::shared_ptr<radfiled3d::CartesianRadiationField> field, std::shared_ptr<radfiled3d::storage::v1::RadiationFieldMetadata> metadata, const std::string& file)
{
	radfiled3d::FileLock lock(file);

	if (!fs::exists(file)) {
		store_radiation_field_atomically(field, metadata, file);
		return;
	}

	auto existing_metadata = std::dynamic_pointer_cast<radfiled3d::storage::v1::RadiationFieldMetadata>(radfiled3d::storage::FieldStore::load_metadata(file));
	require(existing_metadata != nullptr, "'" + file + "' has no V1 metadata.");
	auto existing = std::dynamic_pointer_cast<radfiled3d::CartesianRadiationField>(radfiled3d::storage::FieldStore::load(file));
	require(existing != nullptr, "'" + file + "' does not hold a cartesian field.");
	require(existing->get_voxel_counts() == field->get_voxel_counts() && existing->get_voxel_dimensions() == field->get_voxel_dimensions(), "the voxel grid differs.");

	auto header = metadata->get_header();
	const auto& existing_header = existing_metadata->get_header();
	check_same_simulation(existing_header.simulation, header.simulation);

	const double existing_primaries = static_cast<double>(existing_header.simulation.primary_particle_count);
	const double added_primaries = static_cast<double>(header.simulation.primary_particle_count);
	require(existing_primaries + added_primaries > 0.0, "neither field holds any primary particles.");
	const double r = added_primaries / (existing_primaries + added_primaries);

	for (const auto& [name, added_channel] : field->get_channels()) {
		require(existing->has_channel(name), "the file has no channel '" + name + "' (was it post-processed, e.g. its channels joined?).");
		combine_channel(name, *existing->get_channel(name), *added_channel, r);
	}

	header.simulation.primary_particle_count = existing_header.simulation.primary_particle_count + header.simulation.primary_particle_count;
	metadata->set_header(header);
	const std::vector<std::string> existing_keys = existing_metadata->get_dynamic_metadata_keys();
	const std::vector<std::string> added_keys = metadata->get_dynamic_metadata_keys();
	const std::string duration_key = "simulation_duration_s";
	if (std::find(existing_keys.begin(), existing_keys.end(), duration_key) != existing_keys.end() && std::find(added_keys.begin(), added_keys.end(), duration_key) != added_keys.end())
		metadata->get_dynamic_metadata<radfiled3d::ScalarVoxel<uint64_t>>(duration_key) += existing_metadata->get_dynamic_metadata<radfiled3d::ScalarVoxel<uint64_t>>(duration_key).get_data();
	// the tube spectrum counts the generated primaries per energy: the combined run holds the counts of both
	const std::string tube_key = "tube_spectrum";
	if (std::find(existing_keys.begin(), existing_keys.end(), tube_key) != existing_keys.end() && std::find(added_keys.begin(), added_keys.end(), tube_key) != added_keys.end()) {
		auto& existing_tube = existing_metadata->get_dynamic_metadata<radfiled3d::HistogramVoxel<float>>(tube_key);
		auto& added_tube = metadata->get_dynamic_metadata<radfiled3d::HistogramVoxel<float>>(tube_key);
		require(existing_tube.get_bins() == added_tube.get_bins() && existing_tube.get_histogram_bin_width() == added_tube.get_histogram_bin_width(), "the tube spectra have different bins.");
		auto existing_counts = existing_tube.get_histogram();
		auto added_counts = added_tube.get_histogram();
		for (size_t i = 0; i < added_counts.size(); i++)
			added_counts[i] = to_stored(static_cast<double>(existing_counts[i]) + static_cast<double>(added_counts[i]), "tube spectrum");
	}
	const std::string seed_key = "random_seed";
	if (std::find(added_keys.begin(), added_keys.end(), seed_key) != added_keys.end()) {
		auto& added_seed = metadata->get_dynamic_metadata<radfiled3d::ScalarVoxel<uint64_t>>(seed_key);
		if (std::find(existing_keys.begin(), existing_keys.end(), seed_key) != existing_keys.end()) {
			const uint64_t existing_seed = existing_metadata->get_dynamic_metadata<radfiled3d::ScalarVoxel<uint64_t>>(seed_key).get_data();
			require(existing_seed == 0 || existing_seed != added_seed.get_data(), "the run uses the random seed " + std::to_string(existing_seed) + " of the file, so it would repeat its random numbers.");
		}
		added_seed = static_cast<uint64_t>(0);
	}

	store_radiation_field_atomically(existing, metadata, file);
}

void RadiationSimulation::store_radiation_field_atomically(std::shared_ptr<radfiled3d::IRadiationField> field, std::shared_ptr<radfiled3d::storage::RadiationFieldMetadata> metadata, const std::string& file)
{
	const std::string temporary = file + ".store-" + unique_file_token();
	try {
		radfiled3d::storage::FieldStore::store(field, metadata, temporary);
		fs::rename(temporary, file);
	}
	catch (...) {
		std::error_code ignored;
		fs::remove(temporary, ignored);
		throw;
	}
}

std::string RadiationSimulation::unique_file_token()
{
	static std::atomic<uint64_t> counter{ 0 };
	std::ostringstream token;
	token << std::hex << current_process_id() << '-' << std::chrono::system_clock::now().time_since_epoch().count() << '-' << counter++;
	return token.str();
}
