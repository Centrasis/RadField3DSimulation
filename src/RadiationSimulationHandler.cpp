#include "RadiationSimulationHandler.hpp"
#include <G4RunManagerFactory.hh>
#ifdef WITH_GEANT4_UIVIS
#include <G4UImanager.hh>
#endif
#include <thread>
#include "Geant4/G4SceneConstructor.hpp"
#include <QGSP_BIC_HP.hh>
#include <G4UIExecutive.hh>
#include <G4UIsession.hh>
#include "Geant4/G4RadiationSource.hpp"
#include "World.hpp"
#include "Geant4/G4World.hpp"
#include "Geant4/G4RadiationFieldDetector.hpp"
#include "RadFiled3D/storage/RadiationFieldStore.hpp"
#include <G4SteppingVerbose.hh>
#include "G4StepLimiterPhysics.hh"
#include "Randomize.hh"
#include <random>
#include "G4EmStandardPhysics_option4.hh"
#include "Geant4/G4PhysicsList.hpp"


using namespace RadiationSimulation;

bool G4RadiationSimulationHandler::initialize()
{
	CLHEP::HepRandom::setTheEngine(&this->random_generator);
	G4SteppingVerbose::UseBestUnit(4);
#ifdef WITH_GEANT4_UIVIS
	// Switch to Serial Mode to be able to use the GUI to draw tracks
	this->G4mgr = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial));
#else
	this->G4mgr = std::unique_ptr<G4RunManager>(G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT));
	G4int nThreads = std::max(G4Threading::G4GetNumberOfCores(), 1);
	this->G4mgr->SetNumberOfThreads(nThreads);
#endif

	this->physics = new MedicalPhysicsList();
	this->G4mgr->SetUserInitialization(this->physics);

	return true;
}

void RadiationSimulation::G4RadiationSimulationHandler::finalize()
{
	if (World::Get()->get_radiation_field_detector().get() == NULL) {
		auto rad_det = std::make_shared<G4RadiationFieldDetector>(
			this->radiation_field_resolution.radiation_field_dimensions,
			this->radiation_field_resolution.radiation_field_voxel_dimensions,
			static_cast<size_t>(this->radiation_field_resolution.radiation_field_max_energy / this->radiation_field_resolution.energy_resolution),
			static_cast<double>(this->radiation_field_resolution.energy_resolution),
			0.1f
		);
		rad_det->register_on_new_particle([=](size_t evt_count, const G4Step* step) {
			for (auto& cb : this->callbacks) {
				if (evt_count > 0 && evt_count % cb.first == 0) {
					std::shared_ptr<RadFiled3D::IRadiationField> field = rad_det->get_normalized_field_copy();
					cb.second(field, evt_count);
				}
			}
		});
		World::Get()->set_radiation_field_detector(
			rad_det
		);

		this->G4mgr->SetUserInitialization(
			new G4RadiationFieldAction(
				rad_det,
				std::make_shared<G4RadiationSource>(
					World::Get()->get_radiation_source()
				)
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

	this->field_detector = std::dynamic_pointer_cast<G4RadiationFieldDetector>(World::Get()->get_radiation_field_detector());
}

void RadiationSimulation::G4RadiationSimulationHandler::display_gui()
{
#ifdef WITH_GEANT4_UIVIS
	auto ui = std::make_unique<G4UIExecutive>(0, NULL); // Use unique_ptr
	this->has_ui = true;
	this->G4VisManager = std::make_unique<G4VisExecutive>();
	this->G4VisManager->Initialize();
	this->G4UIManager = std::shared_ptr<G4UImanager>(G4UImanager::GetUIpointer());
	
	this->G4UIManager->ApplyCommand("/vis/scene/add/trajectories 0");
	this->G4UIManager->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID");
	this->G4UIManager->ApplyCommand("/vis/modeling/trajectories/drawByParticleID-0/set all false");
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

void RadiationSimulation::G4RadiationSimulationHandler::update_gui()
{
#ifdef WITH_GEANT4_UIVIS
	if (has_ui) {
		this->G4UIManager->ApplyCommand("/vis/viewer/flush");
		this->G4UIManager->ApplyCommand("/vis/viewer/rebuild");
	}
#endif
}

void RadiationSimulation::G4RadiationSimulationHandler::add_geometry(const std::vector<std::shared_ptr<Mesh>>& meshes)
{
	this->G4mgr->SetUserInitialization(new G4SceneConstructor(meshes));
	if (this->run_mgr_initialized) {
		this->G4mgr->ReinitializeGeometry();
	}
}

std::shared_ptr<RadFiled3D::IRadiationField> RadiationSimulation::G4RadiationSimulationHandler::simulate_radiation_field(size_t n_particles, RadFiled3D::GridTracerAlgorithm tracing_algorithm)
{
	G4cout << "Particles to calculate: " << n_particles << G4endl;
	if (this->field_detector) {
		switch (tracing_algorithm) {
		case RadFiled3D::GridTracerAlgorithm::SAMPLING:
			this->field_detector->define_grid_tracer<RadFiled3D::SamplingGridTracer>();
			break;
		case RadFiled3D::GridTracerAlgorithm::BRESENHAM:
			this->field_detector->define_grid_tracer<RadFiled3D::BresenhamGridTracer>();
			break;
		case RadFiled3D::GridTracerAlgorithm::LINETRACING:
			this->field_detector->define_grid_tracer<RadFiled3D::LinetracingGridTracer>();
			break;
		}
		this->field_detector->finalize(n_particles);
	}
	this->G4mgr->BeamOn(n_particles);
	this->update_gui();

	return (World::Get()->get_radiation_field_detector()) ? World::Get()->get_radiation_field_detector()->evaluate() : std::shared_ptr<RadFiled3D::IRadiationField>(NULL);
}

void G4RadiationSimulationHandler::deinitialize()
{
#ifdef WITH_GEANT4_UIVIS
	this->G4VisManager.reset();
	this->G4UIManager.reset();
#endif
	this->G4mgr.reset();
}

void RadiationSimulation::G4RadiationSimulationHandler::add_callback_every_n_particles(std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)> callback, size_t n_particles)
{
	callbacks.push_back({ n_particles, callback });
}

void RadiationSimulation::RadiationSimulationHandler::set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution)
{
	this->radiation_field_resolution.radiation_field_dimensions = radiation_field_dimensions;
	this->radiation_field_resolution.radiation_field_voxel_dimensions = radiation_field_voxel_dimensions;
	this->radiation_field_resolution.radiation_field_max_energy = radiation_field_max_energy;
	this->radiation_field_resolution.energy_resolution = energy_resolution;
}
