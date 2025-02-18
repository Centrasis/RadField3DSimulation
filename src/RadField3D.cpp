#include "RadiationSimulation.hpp"
#include "GeometryLoader.hpp"
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <chrono>
#include <G4ios.hh>
#include <RadFiled3D/storage/RadiationFieldStore.hpp>
#include <RadFiled3D/GridTracer.hpp>
#include <Geant4/G4World.hpp>
#include <G4SystemOfUnits.hh>
#include <RadiationFieldDetector.hpp>
#if defined _WIN32 || defined _WIN64
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#include <glm/gtc/quaternion.hpp>
#include <stdexcept>
#include <fstream>


using namespace RadiationSimulation;


void store_radiation_field(std::shared_ptr<RadFiled3D::IRadiationField> field, fs::path out_path, size_t n_particles, std::shared_ptr<XRaySource> source, const std::string& geometry_file, const std::string& spectrum_file, const glm::vec3& source_dir, float source_distance, float xray_energy, bool should_append_to_file, long long start_time) {
	long long end_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	auto metadata = std::make_shared<RadFiled3D::Storage::V1::RadiationFieldMetadata>(
		RadFiled3D::Storage::FiledTypes::V1::RadiationFieldMetadataHeader::Simulation(
			n_particles,
			geometry_file,
			"QGSP_BIC_HP+StdPhysics_Option4",
			RadFiled3D::Storage::FiledTypes::V1::RadiationFieldMetadataHeader::Simulation::XRayTube(
				source_dir,
				-source_dir * source_distance,
				xray_energy,
				spectrum_file
			)
		),
		RadFiled3D::Storage::FiledTypes::V1::RadiationFieldMetadataHeader::Software(
			"RadField3D",
			"1.0.0",
			"https://github.com/Centrasis/RadField3DSimulation",
			"HEAD"
		)
	);

	//fs::path new_path = out_path;
	//new_path.replace_extension("." + std::to_string(n_particles) + ".rf");

	RadFiled3D::HistogramVoxel spectrum = static_cast<XRaySpectrumSource*>(source.get())->getGeneratedSpectrum();
	metadata->set_dynamic_metadata<RadFiled3D::HistogramVoxel>("tube_spectrum", spectrum);
	uint64_t duration = end_time - start_time;
	metadata->set_dynamic_metadata<RadFiled3D::ScalarVoxel<uint64_t>>("simulation_duration_s", RadFiled3D::ScalarVoxel<uint64_t>(&duration));

	RadFiled3D::CartesianRadiationField* cfield = (RadFiled3D::CartesianRadiationField*)field.get();
	if (should_append_to_file) {
		RadFiled3D::Storage::FieldStore::join(field, metadata, out_path.string(), RadFiled3D::Storage::FieldJoinMode::Add, RadFiled3D::Storage::FieldJoinCheckMode::MetadataSimulationSimilar);
	}
	else {
		RadFiled3D::Storage::FieldStore::store(field, metadata, out_path.string());
	}
}


int main(int argc, char* argv[]) {
	fs::path geometry_file = "";
	fs::path spectrum_file;
	double xray_energy = 0.0;
	double max_energy = 0.0;
	float source_angle_alpha = 0.f;
	float source_angle_beta = 0.f;
	float source_distance = 1.f;
	glm::vec3 source_opening_angle = glm::vec3(20.f);
	glm::vec3 world_dim = glm::vec3(1.f);
#ifdef WITH_GEANT4_UIVIS
	bool show_gui = false;
#endif
	bool should_append_to_file = false;
	fs::path out_path;
	size_t particle_count = 1e+6;
	float voxel_dim = 0.1f;
	float energy_resolution = 1e+3;
	std::string source_shape = "cone";
	RadFiled3D::GridTracerAlgorithm tracing_algorithm = RadFiled3D::GridTracerAlgorithm::SAMPLING;


	if (argc <= 1 || argv[1] == "-h" || argv[1] == "--help") {
		G4cout << "RadiationField Calculator Help:\nThe following parameters should be passed to this program\n" << G4endl;
		G4cout << "  --geom: Path to a geometry file" << G4endl;
		G4cout << "  --out: Path where the radiation field should be stored" << G4endl;
		G4cout << "  --max-energy: Maximum Energy of the X-Raytube in eV to expect" << G4endl;
		G4cout << "  --source-alpha: Y-Rotation of the X-Raytube in deg" << G4endl;
		G4cout << "  --source-beta:  Z-Rotation of the X-Raytube in deg" << G4endl;
		G4cout << "  --particles: Maximum number of particles to process" << G4endl;
		G4cout << "  --voxel-dim: Voxel dimensions in m" << G4endl;
		G4cout << "  --world-dim: World dimensions in m 'x y z'" << G4endl;
		G4cout << "  --source-distance: Distance of the source in m" << G4endl;
		G4cout << "  --source-shape: Type of the radiation field shape. Must be one of ['cone', 'rectangle', 'ellipsoid']" << G4endl;
		G4cout << "  --spectrum: Path to an spectrum file when energy is not explicitly set" << G4endl;
		G4cout << "  --source-opening-angle: Opening angle of the source in deg" << G4endl;
		G4cout << "  --energy-resolution: Resolution of the energy scroring. Effectively equals the bin width of the spectra histograms in eV. Default: 1 keV" << G4endl;
		G4cout << "  --tracing-algorithm: Algorithm to use for the grid tracing. Must be one of ['sampling', 'bresenham', 'linetracing']" << G4endl;
#ifdef WITH_GEANT4_UIVIS
		G4cout << "  --gui: Flag if the Geant4 GUI should be shown" << G4endl;
#else
		G4cout << "  --gui: Flag if the Geant4 GUI should be shown. Not available in this build and will be ignored." << G4endl;
#endif
		G4cout << "  --append: Flag if this simulation data should be appended to an potentially existing file using the SimulationSimilar policy" << G4endl;
		return 0;
	}

	size_t i = 0;
	size_t increment = 1;
	while (i + increment < argc) {
		i += increment;
		increment = 2;
		std::string arg = argv[i];
		std::string value = "";
		if (i < argc - 1)
			value = argv[i + 1];

		if (arg == "--geom") {
			geometry_file = value;
			if (!geometry_file.is_absolute())
				geometry_file = fs::absolute(geometry_file);
			continue;
		}
		if (arg == "--source-shape") {
			source_shape = value;
			continue;
		}
		if (arg == "--spectrum") {
			spectrum_file = value;
			if (!spectrum_file.is_absolute())
				spectrum_file = fs::absolute(spectrum_file);
			continue;
		}
		if (arg == "--out") {
			out_path = value;
			if (!out_path.is_absolute())
				out_path = fs::absolute(out_path);
			continue;
		}
		if (arg == "--tracing-algorithm") {
			if (value == "sampling")
				tracing_algorithm = RadFiled3D::GridTracerAlgorithm::SAMPLING;
			else if (value == "bresenham")
				tracing_algorithm = RadFiled3D::GridTracerAlgorithm::BRESENHAM;
			else if (value == "linetracing")
				tracing_algorithm = RadFiled3D::GridTracerAlgorithm::LINETRACING;
			else {
				G4cerr << "Unknown tracing algorithm: " << value << ". Aborting..." << G4endl;
				return -1;
			}
			continue;
		}
		if (arg == "--max-energy") {
			xray_energy = std::stod(value);
			max_energy = xray_energy;
			continue;
		}
		if (arg == "--source-opening-angle") {
			if (source_shape == "cone") {
				try {
					source_opening_angle.x = std::stod(value);
				}
				catch (std::exception& e) {
					G4cerr << "Invalid argument for source opening angle: " << value << G4endl;
					throw e;
				}
				continue;
			}
			if (source_shape == "rectangle" || source_shape == "ellipsoid") {
				try {
					std::stringstream ssin(value);
					float x, y;
					ssin >> x;
					ssin >> y;
					source_opening_angle = glm::vec3(x, y, 0.f);
				}
				catch (std::exception& e) {
					G4cerr << "Invalid argument for source opening angle: " << value << G4endl;
					throw e;
				}
				continue;
			}
			continue;
		}
		if (arg == "--particles") {
			try {
				particle_count = static_cast<size_t>(std::stof(value));
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for particle count: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--source-alpha") {
			try {
				source_angle_alpha = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for source alpha: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--source-beta") {
			try {
				source_angle_beta = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for source beta: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--voxel-dim") {
			try {
				voxel_dim = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for voxel dimension: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--world-dim") {
			try {
				std::stringstream ssin(value);
				float x, y, z;
				ssin >> x;
				ssin >> y;
				ssin >> z;
				world_dim = glm::vec3(x, y, z);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for world dimension: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--source-distance") {
			try {
				source_distance = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for source distance: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--energy-resolution") {
			try {
				energy_resolution = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for energy resolution: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--gui") {
#ifdef WITH_GEANT4_UIVIS
			show_gui = true;
#endif
			increment = 1;
			continue;
		}
		if (arg == "--append") {
			should_append_to_file = true;
			increment = 1;
			continue;
		}
	}

	// if out_path does not end with .rf3 append it
	if (out_path.extension() != ".rf3") {
		out_path.replace_extension(".rf3");
	}

	G4RadiationSimulationHandler* simulation_handler = static_cast<G4RadiationSimulationHandler*>(RadiationSimulator::initialize(RadiationHandlerType::Geant4MedicalXRay).get());
	RadiationSimulator::set_world_info(std::make_unique<WorldInfo>("Air", world_dim));

	if (!geometry_file.empty()) {
		G4cout << "Attempt to load geometry from: " << geometry_file.string() << G4endl;
		auto meshes = GeometryLoader::Load(geometry_file.string());
		G4cout << "Meshes loaded: " << meshes.size() << G4endl;
		RadiationSimulator::add_geometry(meshes);
	}
	else {
		RadiationSimulator::add_geometry(std::vector<std::shared_ptr<Mesh>>());
		G4cout << "No geometry file specified! Perfoming empty simulation..." << G4endl;
	}

	std::shared_ptr<XRaySource> source = std::shared_ptr<XRaySource>(NULL);

	std::unique_ptr<ISourceShape> shape = NULL;

	if (source_shape == "cone") {
		shape = std::make_unique<ConeSourceShape>(source_opening_angle.x);
		G4cout << "Using source opening angle: " << source_opening_angle.x << "°" << G4endl;
	} else if (source_shape == "rectangle") {
		shape = std::make_unique<RectangleSourceShape>(glm::vec2(source_opening_angle.x, source_opening_angle.y), source_distance);
		G4cout << "Using source rectangular dimensions: " << source_opening_angle.x << "m x " << source_opening_angle.y << "m" << G4endl;
	} else if (source_shape == "ellipsoid") {
		shape = std::make_unique<EllipsoidSourceShape>(glm::vec2(source_opening_angle.x, source_opening_angle.y));
		G4cout << "Using source ellipsoid shape with angles: " << source_opening_angle.x << "° x " << source_opening_angle.y << "°" << G4endl;
	} else {
		G4cout << "Unknown source shape: " << source_shape << ". Aborting..." << G4endl;
		return -1;
	}

	G4cout << "Using source shape: " << source_shape << G4endl;

	if (xray_energy != 0.0 && spectrum_file.empty()) {
		source = std::make_shared<XRaySource>(xray_energy, std::move(shape));
	} else {
		if (spectrum_file.empty()) {
			G4cout << "No spectrum file specified and no energy! Aborting..." << G4endl;
			return -1;
		} else {
			G4cout << "Attempt to load spectrum from: " << spectrum_file.string() << G4endl;
			auto spectrum_probabilities = SpectrumLoader::LoadSpectrum(spectrum_file.string());
			source = std::make_shared<XRaySpectrumSource>(
				spectrum_probabilities,
				std::move(shape),
				500,
				max_energy
			);
			xray_energy = spectrum_probabilities->max();
		}
	}

	G4cout << "Set X-Ray source energy to: " << xray_energy / 1e+3 << "keV" << G4endl;
	G4cout << "Set X-Ray source rotation (" << source_angle_alpha << "°, " << source_angle_beta << "°)" << G4endl;
	G4cout << "Set X-Ray source distance to " << source_distance << " m" << G4endl;
	G4cout << "Set tracking energy maximum to " << max_energy / 1e+3 << "keV" << G4endl;
	
	glm::quat rotation = glm::angleAxis(glm::radians(source_angle_alpha), glm::vec3(0.f, 1.f, 0.f)) * glm::angleAxis(glm::radians(source_angle_beta), glm::vec3(1.f, 0.f, 0.f));
	const glm::vec3 source_dir = rotation * glm::vec3(0.f, 0.f, -1.f);
	source->setTransform(-source_dir * source_distance, source_dir);
	RadiationSimulator::add_radiation_source(source);

	RadiationSimulator::set_radiation_field_resolution(world_dim, glm::vec3(voxel_dim), max_energy * eV, energy_resolution * eV);

	G4cout << "Using world dimensions: " << world_dim.x << "m x " << world_dim.y << "m x " << world_dim.z << "m" << G4endl;
	G4cout << "Using voxel grid tracer: " << (tracing_algorithm == RadFiled3D::GridTracerAlgorithm::SAMPLING ? "SAMPLING" : (tracing_algorithm == RadFiled3D::GridTracerAlgorithm::BRESENHAM ? "BRESENHAM" : "LINETRACING")) << G4endl;
	G4cout << "Start simulating with voxel dimension: " << voxel_dim << "m and an particle count of " << particle_count << G4endl << G4endl;

#ifdef WITH_GEANT4_UIVIS
	if (show_gui) {
		RadiationSimulator::display_gui();
		return 0;
	}
#endif

	if (!should_append_to_file && fs::exists(out_path)) {
		fs::remove(out_path);
	}

	if (should_append_to_file)
		RadFiled3D::Storage::FieldStore::enable_file_lock_syncronization(true);

	long long start_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	G4cout << "Simulating...";
	G4cout << G4endl << "Writing field to: " << out_path.string() << G4endl;
	size_t last_particle_count = 0;
	RadiationSimulator::add_callback_every_n_particles([&](std::shared_ptr<RadFiled3D::IRadiationField> field, size_t n_particles) {
		G4cout << "Simulation auto save after " << n_particles << " particles" << G4endl;
		last_particle_count = n_particles;
		store_radiation_field(field, out_path, n_particles, source, geometry_file.string(), spectrum_file.string(), source_dir, source_distance, xray_energy, should_append_to_file, start_time);
	}, 1e+6);

	auto field_promise = RadiationSimulator::simulate_radiation_field(particle_count, tracing_algorithm);
	auto field = field_promise.get_future().get();
	last_particle_count = G4World::Get()->get_radiation_field_detector()->get_number_of_tracked_particles();
	store_radiation_field(field, out_path, last_particle_count, source, geometry_file.string(), spectrum_file.string(), source_dir, source_distance, xray_energy, should_append_to_file, start_time);

	G4cout << G4endl << "Wrote field to: " << out_path.string() << G4endl;

	std::quick_exit(EXIT_SUCCESS);

	// RadiationSimulator::deinitialize(); double free error.

	return 0;
}