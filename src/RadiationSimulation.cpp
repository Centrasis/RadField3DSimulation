#include "RadiationSimulation.hpp"
#include "World.hpp"
#include "RadFiled3D/storage/RadiationFieldStore.hpp"
#include <stdexcept>

using namespace RadiationSimulation;

std::shared_ptr<RadiationSimulationHandler> RadiationSimulator::handler;
bool RadiationSimulator::bIsBusy = false;


std::shared_ptr<RadiationSimulationHandler> RadiationSimulator::initialize(RadiationHandlerType handler_type)
{
	World::instance = std::make_shared<World>();
	World::world_info = std::make_unique<WorldInfo>("Air", glm::vec3(1.f));
	switch (handler_type) {
		case RadiationHandlerType::Geant4MedicalXRay:
			RadiationSimulator::handler = std::make_shared<G4RadiationSimulationHandler>();
			break;
		default:
			throw std::runtime_error("Unknown RadiationHandlerType!");
	}
	if (!RadiationSimulator::handler->initialize())
		throw std::runtime_error("Handler initialization failed!");
	return RadiationSimulator::handler;
}

void RadiationSimulator::deinitialize()
{
	RadiationSimulator::handler->deinitialize();
}

std::promise<std::shared_ptr<RadFiled3D::IRadiationField>> RadiationSimulator::simulate_radiation_field(size_t n_particles, RadFiled3D::GridTracerAlgorithm tracing_algorithm)
{
	auto promise = std::promise<std::shared_ptr<RadFiled3D::IRadiationField>>();
	RadiationSimulator::handler->finalize();
	RadiationSimulator::bIsBusy = true;
	auto field = RadiationSimulator::handler->simulate_radiation_field(n_particles, tracing_algorithm);
	promise.set_value(field);
	RadiationSimulator::bIsBusy = false;
	return promise;
}

void RadiationSimulation::RadiationSimulator::add_callback_every_n_particles(std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)> callback, size_t n_particles)
{
	RadiationSimulator::handler->add_callback_every_n_particles(callback, n_particles);
}

void RadiationSimulation::RadiationSimulator::display_gui()
{
	RadiationSimulator::handler->finalize();
	RadiationSimulator::handler->display_gui();
}

void RadiationSimulation::RadiationSimulator::add_geometry(const std::vector<std::shared_ptr<Mesh>>& meshes)
{
	RadiationSimulator::handler->add_geometry(meshes);
	World::Get()->set_geometries(meshes);
}

void RadiationSimulation::RadiationSimulator::add_geometry(std::shared_ptr<Mesh> mesh)
{
	RadiationSimulator::add_geometry(std::vector<std::shared_ptr<Mesh>>({ mesh }));
}

void RadiationSimulation::RadiationSimulator::add_radiation_source(std::shared_ptr<RadiationSource> source)
{
	World::Get()->radiation_source = source;
}

void RadiationSimulation::RadiationSimulator::set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution)
{
	RadiationSimulator::handler->set_radiation_field_resolution(radiation_field_dimensions, radiation_field_voxel_dimensions, radiation_field_max_energy, energy_resolution);
}

void RadiationSimulation::RadiationSimulator::set_world_info(std::unique_ptr<RadiationSimulation::WorldInfo> info)
{
	World::world_info = std::move(info);
}
