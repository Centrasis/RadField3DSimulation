#include "RadiationSimulation.hpp"
#include "GeometryLoader.hpp"
#include "FieldAppend.hpp"
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <chrono>
#include <G4ios.hh>
#include <radfiled3d/storage/radiation_field_store.hpp>
#include <radfiled3d/grid_tracer.hpp>
#include <radfiled3d/helpers/typing.hpp>
#include <Geant4/G4World.hpp>
#include <G4SystemOfUnits.hh>
#include <Geant4/G4RadiationFieldDetector.hpp>
#include <Geant4/G4PhysicsList.hpp>
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
namespace G4 = RadiationSimulation::Geant4;


// Seed of the run from the system time; the process id is mixed in (splitmix64 finalizer) so that runs started in the
// same clock tick, e.g. on several cluster nodes, still get unrelated seeds.
uint64_t time_based_seed() {
	uint64_t z = static_cast<uint64_t>(std::hash<std::string>{}(unique_file_token())) + 0x9e3779b97f4a7c15ULL;
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
	return z ^ (z >> 31);
}

// Stores `field` into `out_path`. With `join_into_file`, the field is appended to an existing file (see
// append_radiation_field); otherwise the file is replaced.
void store_radiation_field(std::shared_ptr<radfiled3d::IRadiationField> field, fs::path out_path, size_t n_particles, std::shared_ptr<XRaySource> source, const std::string& geometry_file, const std::string& spectrum_file, const glm::vec3& source_dir, float source_distance, float xray_energy, bool join_into_file, long long start_time, radfiled3d::typing::FieldShape field_shape, uint64_t random_seed) {
	long long end_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

	auto metadata = std::make_shared<radfiled3d::storage::v1::RadiationFieldMetadata>(
		radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation(
			n_particles,
			geometry_file,
			G4::MedicalPhysicsList::getName(),
			radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation::XRayTube(
				source_dir,
				-source_dir * source_distance,
				xray_energy,
				spectrum_file
			)
		),
		radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Software(
			"RadField3D",
			"1.2.0",
			"https://github.com/Centrasis/RadField3DSimulation",
			"HEAD"
		)
	);

	const RadiationSimulation::ISourceShape* actual_field_shape = source->getShape();
	XRaySpectrumSource* spectrum_source = dynamic_cast<XRaySpectrumSource*>(source.get());
	// counted in 64 bit; float only rounds large counts here, it does not clip them
	const std::vector<uint64_t> generated_counts = spectrum_source->getGeneratedCounts();
	metadata->set_dynamic_custom_metadata<radfiled3d::HistogramVoxel<float>>("tube_spectrum", radfiled3d::HistogramVoxel<float>(generated_counts.size(), spectrum_source->getGeneratedSpectrumBinWidth_eV(), nullptr));
	auto tube_spectrum = metadata->get_dynamic_metadata<radfiled3d::HistogramVoxel<float>>("tube_spectrum").get_histogram();
	for (size_t i = 0; i < generated_counts.size(); i++)
		tube_spectrum[i] = static_cast<float>(generated_counts[i]);
	uint64_t duration = end_time - start_time;
	metadata->add_dynamic_metadata<uint64_t>("simulation_duration_s", duration);
	metadata->add_dynamic_metadata<uint8_t>("xray_field_shape", static_cast<uint8_t>(field_shape));
	metadata->add_dynamic_metadata<uint64_t>("random_seed", random_seed);

	float angle_deg = 0.f;
	glm::vec2 field_dims_or_angles;
	if (field_shape == radfiled3d::typing::FieldShape::Cone) {
		angle_deg = dynamic_cast<const RadiationSimulation::ConeSourceShape*>(actual_field_shape)->getOpeningAngleDegrees();
		metadata->add_dynamic_metadata<float>("xray_tube_opening_angle_deg", angle_deg);
	}
	if (field_shape == radfiled3d::typing::FieldShape::Rectangle) {
		field_dims_or_angles = dynamic_cast<const RadiationSimulation::RectangleSourceShape*>(actual_field_shape)->getFieldSizeMeters();
		metadata->add_dynamic_metadata<glm::vec2>("xray_tube_field_rect_dimensions_m", field_dims_or_angles);
	}
	if (field_shape == radfiled3d::typing::FieldShape::Ellipsis) {
		field_dims_or_angles = dynamic_cast<const RadiationSimulation::EllipsoidSourceShape*>(actual_field_shape)->getOpeningAnglesDegrees();
		metadata->add_dynamic_metadata<glm::vec2>("xray_tube_field_ellipsis_opening_angles_deg", field_dims_or_angles);
	}

	// The voxelized geometry is written only when the file is created; later stores update the radiation channels and
	// metadata and keep the file's geometry channel.
	const bool file_exists = fs::exists(out_path);
	if (!file_exists)
		World::Get()->get_radiation_field_detector()->add_geometry_channel_to(*std::static_pointer_cast<radfiled3d::CartesianRadiationField>(field));

	if (join_into_file)
		append_radiation_field(std::static_pointer_cast<radfiled3d::CartesianRadiationField>(field), metadata, out_path.string());
	else if (file_exists)
		radfiled3d::storage::FieldStore::replace(field, metadata, out_path.string());
	else
		store_radiation_field_atomically(field, metadata, out_path.string());
}


int main(int argc, char* argv[]) {
try {
	fs::path geometry_file = "";
	fs::path geometry_desc_file = "";
	fs::path spectrum_file = "";
	double xray_energy = 0.0;
	double max_energy = 0.0;
	float source_angle_phi = 0.f;
	float source_angle_theta = 0.f;
	float source_distance = 1.f;
	radfiled3d::typing::FieldShape field_shape;
	glm::vec3 source_opening_angle = glm::vec3(20.f);
	glm::vec3 world_dim = glm::vec3(1.f);
#ifdef WITH_GEANT4_UIVIS
	bool show_gui = false;
#endif
	bool should_append_to_file = false;
	fs::path out_path;
	size_t particle_count = 1e+6;
	size_t autosave_interval = 1e+7;
	float voxel_dim = 0.1f;
	float energy_resolution = 1e+3;
	float statistical_error_threshold = 0.1f;
	float statistical_error_enforcement_ratio = 0.95f;
	std::string source_shape = "cone";
	std::string world_material = "Air";
	int cpu_count = -1;

	size_t angular_phi_segments = 0;
	uint32_t directional_lobes = 0;
	size_t angular_theta_segments = 0;
	// every voxel a step touches counts; SAMPLING misses voxels an oblique step only clips (about 7 % fewer hits at 10°)
	radfiled3d::GridTracerAlgorithm tracing_algorithm = radfiled3d::GridTracerAlgorithm::LINETRACING;
	bool path_length_weighting = true;


	if (argc <= 1 || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help") {
		G4cout << "RadiationField Calculator Help:\nThe following parameters should be passed to this program\n" << G4endl;
		G4cout << "  --geom: Path to a geometry file" << G4endl;
		G4cout << "  --geom-desc: Path to a geometry description file. Default: <geom[noExt]>.desc" << G4endl;
		G4cout << "  --out: Path where the radiation field should be stored" << G4endl;
		G4cout << "  --max-energy: Maximum Energy of the X-Raytube in eV to expect" << G4endl;
		G4cout << "  --source-phi: Y-Rotation of the X-Raytube in deg" << G4endl;
		G4cout << "  --source-theta:  Z-Rotation of the X-Raytube in deg" << G4endl;
		G4cout << "  --particles: Maximum number of particles to process" << G4endl;
		G4cout << "  --voxel-dim: Voxel dimensions in m" << G4endl;
		G4cout << "  --world-dim: World dimensions in m 'x y z'" << G4endl;
		G4cout << "  --angular-resolution: Optionally enables an extra layer that captures the angular distribution of the flux in each voxel." << G4endl;
		G4cout << "  --path-length-weighting: 'on' (default) or 'off'. With the line tracer, flux, spectrum, angular bins and vMF lobes weight every voxel a photon enters by its path length inside it (track-length estimate, flux unit: path length in voxel edges per primary); 'off' counts every touched voxel once (as before 1.2.0)" << G4endl;
		G4cout << "  --directional-lobes: Maximum number of von Mises-Fisher lobes per voxel (1-8) learned during the run and stored as layer 'vmf_lobes' of scatter_field. Default: 0 (off)" << G4endl;
		G4cout << "  --source-distance: Distance of the source in m" << G4endl;
		G4cout << "  --source-shape: Type of the radiation field shape. Must be one of ['cone', 'rectangle', 'ellipsoid']" << G4endl;
		G4cout << "  --spectrum: Path to an spectrum file when energy is not explicitly set" << G4endl;
		G4cout << "  --source-opening-angle: Opening angle of the source in deg" << G4endl;
		G4cout << "  --energy-resolution: Resolution of the energy scroring. Effectively equals the bin width of the spectra histograms in eV. Default: 1 keV" << G4endl;
		G4cout << "  --tracing-algorithm: Algorithm to use for the grid tracing. Must be one of ['sampling', 'bresenham', 'linetracing']. Default: linetracing" << G4endl;
		G4cout << "  --world-material: Material of the world. Default: Air" << G4endl;
		G4cout << "  --statistical-error-threshold: Threshold of the statistical error to stop the simulation early. Default: 0.1 (10%)" << G4endl;
		G4cout << "  --statistical-error-enforcement-ratio: Ratio of voxels that must be below the statistical error threshold to stop the simulation early. Default: 0.95 (95%)" << G4endl;
#ifdef WITH_GEANT4_UIVIS
		G4cout << "  --gui: Flag if the Geant4 GUI should be shown" << G4endl;
#else
		G4cout << "  --gui: Flag if the Geant4 GUI should be shown. Not available in this build and will be ignored." << G4endl;
#endif
		G4cout << "  --append: Flag if this simulation data should be appended to an potentially existing file using the SimulationSimilar policy" << G4endl;
		G4cout << "  --cpu-count: Number of CPU cores to use. Default: -1 (all available cores)" << G4endl;
		G4cout << "  --autosave-interval: Store the field every N particles. Default: 1e7, 0 disables auto-saves" << G4endl;
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
		if (arg == "--geom-desc") {
			geometry_desc_file = value;
			if (geometry_desc_file.empty()) {
				geometry_desc_file = geometry_file;
				geometry_desc_file.replace_extension(".desc");
			}
			if (!geometry_desc_file.is_absolute())
				geometry_desc_file = fs::absolute(geometry_desc_file);
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
		if (arg == "--world-material") {
			world_material = value;
			continue;
		}
		if (arg == "--out") {
			out_path = value;
			if (!out_path.is_absolute())
				out_path = fs::absolute(out_path);
			continue;
		}
		if (arg == "--path-length-weighting") {
			if (value != "on" && value != "off") {
				G4cerr << "--path-length-weighting must be 'on' or 'off', got: " << value << ". Aborting..." << G4endl;
				return -1;
			}
			path_length_weighting = value == "on";
			continue;
		}
		if (arg == "--tracing-algorithm") {
			if (value == "sampling")
				tracing_algorithm = radfiled3d::GridTracerAlgorithm::SAMPLING;
			else if (value == "bresenham")
				tracing_algorithm = radfiled3d::GridTracerAlgorithm::BRESENHAM;
			else if (value == "linetracing")
				tracing_algorithm = radfiled3d::GridTracerAlgorithm::LINETRACING;
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
		if (arg == "--cpu-count") {
			cpu_count = std::stoi(value);
			continue;
		}
		if (arg == "--autosave-interval") {
			const double interval = std::stod(value);
			if (!(interval >= 0.0))
				throw std::invalid_argument("--autosave-interval must be >= 0, got " + value);
			autosave_interval = static_cast<size_t>(interval);
			continue;
		}

		if (arg == "--directional-lobes") {
			const int lobes = std::stoi(value);
			if (lobes < 0 || lobes > 8)
				throw std::invalid_argument("--directional-lobes must lie in [0, 8], got " + value);
			directional_lobes = static_cast<uint32_t>(lobes);
			continue;
		}
		if (arg == "--angular-resolution") {
			try {
				std::stringstream ssin(value);
				size_t phi, theta;
				ssin >> phi;
				ssin >> theta;
				angular_phi_segments = phi;
				angular_theta_segments = theta;
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for angular resolution: " << value << G4endl;
				throw e;
			}
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
		if (arg == "--source-phi") {
			try {
				source_angle_phi = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for source alpha: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--source-theta") {
			try {
				source_angle_theta = std::stof(value);
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
		if (arg == "--statistical-error-threshold") {
			try {
				statistical_error_threshold = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for statistical error threshold: " << value << G4endl;
				throw e;
			}
			continue;
		}
		if (arg == "--statistical-error-enforcement-ratio") {
			try {
				statistical_error_enforcement_ratio = std::stof(value);
			}
			catch (std::exception& e) {
				G4cerr << "Invalid argument for statistical error enforcement ratio: " << value << G4endl;
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

	// Verify the Geant4 datasets this simulation needs are reachable before initialising Geant4; a missing
	// dataset otherwise aborts deep inside Geant4 with a cryptic G4Exception. Report clearly and exit.
	{
		const char* required_datasets[] = { "G4LEDATA", "G4LEVELGAMMADATA", "G4ENSDFSTATEDATA", "G4PARTICLEXSDATA" };
		std::vector<std::string> missing;
		for (const char* var : required_datasets) {
			const char* path = std::getenv(var);
			if (path == nullptr || !fs::exists(path))
				missing.push_back(var);
		}
		if (!missing.empty()) {
			std::cerr << "ERROR: required Geant4 data is not available. The following dataset environment variables are"
			             " unset or point to a missing directory:" << std::endl;
			for (const std::string& var : missing)
				std::cerr << "  - " << var << std::endl;
			std::cerr << "Install the Geant4 datasets and source the Geant4 environment (e.g. 'source"
			             " <geant4-install>/bin/geant4.sh', or export the G4*DATA variables; run 'geant4-config"
			             " --datasets' to list them), then retry." << std::endl;
			return 1;
		}
	}

	if (cpu_count > 0) {
		G4cout << "Initialize radiation simulation handler to use " << cpu_count << " threads." << G4endl;
	}
	else {
		G4cout << "Initialize radiation simulation handler to use all available threads." << G4endl;
	}
	G4::RadiationSimulationHandler* simulation_handler = RadiationSimulator::initialize(cpu_count).get();
	const uint64_t random_seed = time_based_seed();
	RadiationSimulator::set_random_seed(random_seed);
	G4cout << "Random seed: " << random_seed << G4endl;
	RadiationSimulator::set_world_info(std::make_unique<WorldInfo>(world_material, world_dim));

	std::vector<std::shared_ptr<Geometry::Mesh>> meshes;
	if (!geometry_file.empty()) {
		G4cout << "Attempt to load geometry from: " << geometry_file.string() << G4endl;
		if (!geometry_desc_file.empty())
			G4cout << "Attempt to load geometry description from: " << geometry_desc_file.string() << G4endl;
		meshes = GeometryLoader::Load(geometry_file.string(), geometry_desc_file.string());
		G4cout << "Meshes loaded: " << meshes.size() << G4endl;
	}
	else {
		G4cout << "No geometry file specified! Perfoming empty simulation..." << G4endl;
	}

	std::shared_ptr<XRaySource> source = std::shared_ptr<XRaySource>(NULL);

	std::unique_ptr<ISourceShape> shape = NULL;

	if (source_shape == "cone") {
		shape = std::make_unique<ConeSourceShape>(source_opening_angle.x);
		field_shape = radfiled3d::typing::FieldShape::Cone;
		G4cout << "Using source opening angle: " << source_opening_angle.x << "°" << G4endl;
	} else if (source_shape == "rectangle") {
		shape = std::make_unique<RectangleSourceShape>(glm::vec2(source_opening_angle.x, source_opening_angle.y), source_distance);
		field_shape = radfiled3d::typing::FieldShape::Rectangle;
		G4cout << "Using source rectangular dimensions: " << source_opening_angle.x << "m x " << source_opening_angle.y << "m" << G4endl;
	} else if (source_shape == "ellipsoid") {
		shape = std::make_unique<EllipsoidSourceShape>(glm::vec2(source_opening_angle.x, source_opening_angle.y));
		field_shape = radfiled3d::typing::FieldShape::Ellipsis;
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

	if (source_angle_theta > 180.f || source_angle_theta < 0.f) {
		G4cerr << "Source theta angle out of bounds. Only [0°, 180°] allowed. Aborting..." << G4endl;
		return -1;
	}
	if (source_angle_phi < 0.f || source_angle_phi > 360.f) {
		G4cerr << "Source phi angle out of bounds. Must be in [0°, 360°]. Aborting..." << G4endl;
		return -1;
	}

	G4cout << "Set X-Ray source energy to: " << xray_energy / 1e+3 << "keV" << G4endl;
	G4cout << "Set X-Ray source rotation (" << source_angle_phi << "°, " << source_angle_theta << "°)" << G4endl;
	G4cout << "Set X-Ray source distance to " << source_distance << " m" << G4endl;
	G4cout << "Set tracking energy maximum to " << max_energy / 1e+3 << "keV" << G4endl;

	// Create rotation quaternions:
	// - phi: rotation in XZ plane (around Y axis)
	// - theta: rotation in YZ plane (around X axis)
	// Apply alpha first (XZ plane), then beta (YZ plane)
	glm::quat alpha_rotation = glm::angleAxis(glm::radians(source_angle_phi), glm::vec3(0.f, 1.f, 0.f));
	glm::quat beta_rotation = glm::angleAxis(glm::radians(source_angle_theta), glm::vec3(1.f, 0.f, 0.f));
	glm::quat combined_rotation = beta_rotation * alpha_rotation;
	
	// Start with direction pointing toward center (negative Z direction)
	// and apply the combined rotation
	const glm::vec3 source_dir = combined_rotation * glm::vec3(0.f, 0.f, -1.f);
	source->setTransform(-source_dir * source_distance, source_dir);
	RadiationSimulator::add_radiation_source(source);
	RadiationSimulator::add_geometry(meshes);

	RadiationSimulator::set_radiation_field_resolution(
		world_dim,
		glm::vec3(voxel_dim),
		max_energy * eV,
		energy_resolution * eV,
		statistical_error_threshold,
		statistical_error_enforcement_ratio,
		glm::uvec2(angular_phi_segments, angular_theta_segments),
		directional_lobes
	);

	G4cout << "Using world dimensions: " << world_dim.x << "m x " << world_dim.y << "m x " << world_dim.z << "m" << G4endl;
	G4cout << "Using world material: " << world_material << G4endl;
	G4cout << "Using voxel grid tracer: " << (tracing_algorithm == radfiled3d::GridTracerAlgorithm::SAMPLING ? "SAMPLING" : (tracing_algorithm == radfiled3d::GridTracerAlgorithm::BRESENHAM ? "BRESENHAM" : "LINETRACING"))
		<< ((tracing_algorithm == radfiled3d::GridTracerAlgorithm::LINETRACING && path_length_weighting) ? " (path-length weighted)" : "") << G4endl;
	G4cout << "Start simulating with voxel dimension: " << voxel_dim << "m and an particle count of " << particle_count << G4endl << G4endl;

	if (!should_append_to_file && fs::exists(out_path)) {
		fs::remove(out_path);
	}

	// Auto-saves hold this run's cumulative result. When appending, they go to a checkpoint only this run writes, so the
	// output file receives this run exactly once (joined at the end) and a crashed run leaves its progress there.
	const fs::path checkpoint_path = should_append_to_file ? fs::path(out_path.string() + ".partial-" + unique_file_token()) : out_path;

	G4cout << "Simulating...";
	G4cout << G4endl << "Writing field to: " << out_path.string() << G4endl;
	if (should_append_to_file)
		G4cout << "Auto-saves of this run go to: " << checkpoint_path.string() << " (joined into the field when the run finishes)" << G4endl;
	size_t last_particle_count = 0;
	long long start_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	if (autosave_interval > 0) {
		RadiationSimulator::add_callback_every_n_particles([&](std::shared_ptr<radfiled3d::IRadiationField> field, size_t n_particles) {
			G4cout << "Simulation auto save after " << n_particles << " particles" << G4endl;
			last_particle_count = n_particles;
			store_radiation_field(field, checkpoint_path, n_particles, source, geometry_file.string(), spectrum_file.string(), source_dir, source_distance, xray_energy, false, start_time, field_shape, random_seed);
		}, autosave_interval);
	}

#ifdef WITH_GEANT4_UIVIS
	if (show_gui) {
		RadiationSimulator::display_gui();
		std::quick_exit(EXIT_SUCCESS);
	}
#endif

	auto field = RadiationSimulator::simulate_radiation_field(particle_count, tracing_algorithm, path_length_weighting);
	last_particle_count = G4::World::Get()->get_radiation_field_detector()->get_number_of_tracked_particles();
	try {
		store_radiation_field(field, out_path, last_particle_count, source, geometry_file.string(), spectrum_file.string(), source_dir, source_distance, xray_energy, should_append_to_file, start_time, field_shape, random_seed);
	}
	catch (const std::exception&) {
		if (checkpoint_path == out_path)
			throw;
		// the join failed and left the output file unchanged: keep this run's result instead of losing it
		store_radiation_field(field, checkpoint_path, last_particle_count, source, geometry_file.string(), spectrum_file.string(), source_dir, source_distance, xray_energy, false, start_time, field_shape, random_seed);
		std::cerr << "Could not append to " << out_path.string() << "; this run's field was kept in " << checkpoint_path.string() << std::endl;
		throw;
	}
	if (checkpoint_path != out_path && fs::exists(checkpoint_path))
		fs::remove(checkpoint_path);

	G4cout << G4endl << "Wrote field to: " << out_path.string() << G4endl;

	// Tear down Geant4 while all its globals are still alive; destroying the run manager after main() returns
	// (via the static handler's destructor) would touch already-freed Geant4 globals.
	RadiationSimulator::deinitialize();

	return 0;
}
catch (const std::exception& e) {
	// Turn any configuration/parse/simulation error into a clear message + non-zero exit instead of an
	// uncaught-exception abort (core dump). Tear Geant4 down here too (no-op if it was never initialised) so
	// the run manager is not destroyed at static-exit time, when Geant4's globals are already gone.
	std::cerr << "ERROR: " << e.what() << std::endl;
	RadiationSimulator::deinitialize();
	return 1;
}
}