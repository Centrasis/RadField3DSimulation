#include "RadiationSimulationHandler.hpp"
#include <G4RunManagerFactory.hh>
#include <thread>
#include <cmath>
#include "Geant4/G4SceneConstructor.hpp"
#include <QGSP_BIC_HP.hh>
#include <G4UIExecutive.hh>
#include <G4UIsession.hh>
#include "Geant4/G4RadiationSource.hpp"
#include "World.hpp"
#include "Geant4/G4World.hpp"
#include "Geant4/G4RadiationFieldDetector.hpp"
#include "radfiled3d/storage/radiation_field_store.hpp"
#include <G4SteppingVerbose.hh>
#include "G4StepLimiterPhysics.hh"
#include "Randomize.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Gamma.hh"
#include "Geant4/G4PhysicsList.hpp"
#include <glm/gtc/quaternion.hpp>


using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Geant4;

RadiationSimulation::Geant4::RadiationSimulationHandler::RadiationSimulationHandler(const int cpu_count)
	: cpu_count(cpu_count),
	  physics(nullptr),
	  run_mgr_initialized(false),
	  has_ui(false),
	  field_detector(nullptr)
{
}

bool Geant4::RadiationSimulationHandler::initialize()
{
	G4SteppingVerbose::UseBestUnit(4);
#ifdef WITH_GEANT4_UIVIS
	if (cpu_count > 1)
		throw std::runtime_error("Cannot use multiple threads with GUI enabled!");
	// Switch to Serial Mode to be able to use the GUI to draw tracks
	this->G4mgr = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial));
#else
	this->G4mgr = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT));
	G4int nThreads = std::max((this->cpu_count > 0) ? std::min(G4Threading::G4GetNumberOfCores(), this->cpu_count) : G4Threading::G4GetNumberOfCores(), 1);
	this->G4mgr->SetNumberOfThreads(nThreads);
#endif

	this->physics = new MedicalPhysicsList();
	this->G4mgr->SetUserInitialization(this->physics);

	return true;
}

void Geant4::RadiationSimulationHandler::set_random_seed(uint64_t seed)
{
	// Geant4's default engine (MixMaxRng) takes the seed as two 32-bit words (long is 32 bit on Windows).
	const long seeds[2] = { static_cast<long>(seed & 0xffffffffu), static_cast<long>(seed >> 32) };
	G4Random::setTheSeeds(seeds);
}

void RadiationSimulation::Geant4::RadiationSimulationHandler::finalize()
{
	if (this->meshes.size() > 0) {
		if (this->run_mgr_initialized) {
			//this->G4mgr->ReinitializeGeometry();
		}
		else {
			this->G4mgr->SetUserInitialization(new Geant4::SceneConstructor(this->meshes));
		}
	}

	if (RadiationSimulation::World::Get()->get_radiation_field_detector().get() == NULL) {
		// The app is the sole owner of the detector, which is shared across MT workers as one field. Geant4
		// owns only the per-worker Geant4::RadiationFieldSteppingAction forwarders that route steps into it (see
		// Geant4::RadiationFieldAction::Build).
		auto rad_det = std::make_shared<Geant4::RadiationFieldDetector>(
			this->radiation_field_resolution.radiation_field_dimensions,
			this->radiation_field_resolution.radiation_field_voxel_dimensions,
			// round, don't truncate: 0.12f/0.001f = 119.999992 in float would yield 119 bins for 120 keV / 1 keV
			static_cast<size_t>(std::round(this->radiation_field_resolution.radiation_field_max_energy / this->radiation_field_resolution.energy_resolution)),
			static_cast<double>(this->radiation_field_resolution.energy_resolution),
			this->radiation_field_resolution.statistical_error.threshold,
			this->radiation_field_resolution.statistical_error.enforcement_ratio,
			this->radiation_field_resolution.statistical_error.enforcement_resolution,
			this->radiation_field_resolution.angular_resolution,
			this->radiation_field_resolution.directional_lobes
		);
		rad_det->register_on_new_particle([=, this](size_t evt_count, const G4Step* step) {
			for (auto& cb : this->callbacks) {
				if (evt_count > 0 && evt_count % cb.first == 0) {
					std::shared_ptr<radfiled3d::IRadiationField> field = rad_det->get_normalized_field_copy();
					cb.second(field, evt_count);
				}
			}
		});
		RadiationSimulation::World::Get()->set_radiation_field_detector(
			rad_det
		);

		this->G4mgr->SetUserInitialization(
			new Geant4::RadiationFieldAction(
				rad_det,
				RadiationSimulation::World::Get()->get_radiation_source()   // physics source; each worker builds its OWN gun in Build()
			)
		);
	}

	this->G4mgr->SetVerboseLevel(2);
	this->physics->SetDefaultCutValue(0.1 * mm);
	this->physics->SetVerboseLevel(1);

	if (!this->run_mgr_initialized) {
		this->G4mgr->Initialize();
		this->run_mgr_initialized = true;
		auto processes = G4Gamma::Definition()->GetProcessManager()->GetProcessList();
		std::cout << "Processes involved: " << std::endl;
		for(size_t i = 0; i < processes->size(); i++)
			std::cout << (*processes)[i]->GetProcessName() << std::endl;
	}

	this->field_detector = RadiationSimulation::World::Get()->get_radiation_field_detector();
}

void RadiationSimulation::Geant4::RadiationSimulationHandler::display_gui()
{
#ifdef WITH_GEANT4_UIVIS
	G4cout << "Displaying GUI..." << G4endl;
	int dummy_argc = 1;
	char* dummy_argv[] = { (char*)"RadField3D" };
	auto ui = std::unique_ptr<G4UIExecutive>(new G4UIExecutive(dummy_argc, dummy_argv));
	this->has_ui = true;
	this->G4VisManager = std::make_unique<G4VisExecutive>();
	this->G4VisManager->Initialize();
	// NON-OWNING: G4UImanager is a Geant4-owned singleton — an owning shared_ptr would delete it (double-free).
	this->G4UIManager = std::shared_ptr<G4UImanager>(G4UImanager::GetUIpointer(), [](G4UImanager*) {});
	
	this->G4UIManager->ApplyCommand("/vis/scene/add/trajectories 0");
	this->G4UIManager->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");
	this->G4UIManager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set all false");
#else
	G4cout << "Can't display GUI as no OpenGL was linked! Skipping!" << G4endl;
#endif
	this->finalize();
#ifdef WITH_GEANT4_UIVIS
	// Show the OGL Context and define how to display the particles
	this->G4UIManager->ApplyCommand("/vis/open OGL 600x600-0+0");
	this->G4UIManager->ApplyCommand("/vis/viewer/set/autoRefresh false");
	this->G4UIManager->ApplyCommand("/vis/verbose errors");
	this->G4UIManager->ApplyCommand("/vis/drawVolume");
	this->G4UIManager->ApplyCommand("/vis/viewer/set/viewpointThetaPhi -90. 0.");
	this->G4UIManager->ApplyCommand("/vis/viewer/zoom 1.4");
	this->G4UIManager->ApplyCommand("/vis/scene/add/hits");
	this->G4UIManager->ApplyCommand("/vis/scene/endOfEventAction accumulate");
	this->G4UIManager->ApplyCommand("/vis/viewer/set/autoRefresh true");
	this->G4UIManager->ApplyCommand("/vis/verbose warnings");

	this->G4UIManager->ApplyCommand("/vis/enable");
	this->update_gui();
	ui->SessionStart();
	// No need to delete ui, unique_ptr will handle it
#endif
}

void RadiationSimulation::Geant4::RadiationSimulationHandler::update_gui()
{
#ifdef WITH_GEANT4_UIVIS
	if (has_ui) {
		this->G4UIManager->ApplyCommand("/vis/viewer/flush");
		this->G4UIManager->ApplyCommand("/vis/viewer/rebuild");
	}
#endif
}

void RadiationSimulation::Geant4::RadiationSimulationHandler::add_geometry(const std::vector<std::shared_ptr<Geometry::Mesh>>& meshes)
{
	for (auto& m : meshes)
		this->meshes.push_back(m);
}

std::shared_ptr<radfiled3d::IRadiationField> RadiationSimulation::Geant4::RadiationSimulationHandler::simulate_radiation_field(size_t n_particles, radfiled3d::GridTracerAlgorithm tracing_algorithm, bool path_length_weighting)
{
	G4cout << "Particles to calculate: " << n_particles << G4endl;
	if (this->field_detector) {
		switch (tracing_algorithm) {
		case radfiled3d::GridTracerAlgorithm::SAMPLING:
			this->field_detector->define_grid_tracer<radfiled3d::SamplingGridTracer>(path_length_weighting);
			break;
		case radfiled3d::GridTracerAlgorithm::BRESENHAM:
			this->field_detector->define_grid_tracer<radfiled3d::BresenhamGridTracer>(path_length_weighting);
			break;
		case radfiled3d::GridTracerAlgorithm::LINETRACING:
			this->field_detector->define_grid_tracer<radfiled3d::LinetracingGridTracer>(path_length_weighting);
			break;
		}
		this->field_detector->finalize(n_particles);

		// Voxelized before the run, so that every stored field (auto-saves included) carries the geometry and a
		// failing voxelization cannot discard simulated particles.
		const std::shared_ptr<Geant4::World> world = Geant4::World::Get();
		if (!this->geometry_voxelized && world && world->get_volume()) {
			this->field_detector->voxelize_geometry(*world->get_volume(), this->cpu_count);
			this->geometry_voxelized = true;
		}
	}
	this->G4mgr->BeamOn(n_particles);
	this->update_gui();

	return (RadiationSimulation::World::Get()->get_radiation_field_detector()) ? RadiationSimulation::World::Get()->get_radiation_field_detector()->evaluate() : std::shared_ptr<radfiled3d::IRadiationField>(NULL);
}

void Geant4::RadiationSimulationHandler::deinitialize()
{
#ifdef WITH_GEANT4_UIVIS
	this->G4VisManager.reset();
	this->G4UIManager.reset();
#endif
	this->G4mgr.reset();
}

void RadiationSimulation::Geant4::RadiationSimulationHandler::add_callback_every_n_particles(std::function<void(std::shared_ptr<radfiled3d::IRadiationField>, size_t)> callback, size_t n_particles)
{
	callbacks.push_back({ n_particles, callback });
}

void RadiationSimulation::Geant4::RadiationSimulationHandler::set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution, float statistical_error_threshold, float statistical_error_enforcement_ratio, glm::uvec2 angular_resolution, uint32_t directional_lobes)
{
	this->radiation_field_resolution.radiation_field_dimensions = radiation_field_dimensions;
	this->radiation_field_resolution.radiation_field_voxel_dimensions = radiation_field_voxel_dimensions;
	this->radiation_field_resolution.radiation_field_max_energy = radiation_field_max_energy;
	this->radiation_field_resolution.energy_resolution = energy_resolution;
	this->radiation_field_resolution.statistical_error.threshold = statistical_error_threshold;
	this->radiation_field_resolution.statistical_error.enforcement_ratio = statistical_error_enforcement_ratio;
	this->radiation_field_resolution.angular_resolution = angular_resolution;
	this->radiation_field_resolution.directional_lobes = directional_lobes;
}
