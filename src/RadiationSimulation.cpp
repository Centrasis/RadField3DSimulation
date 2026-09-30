#include "RadiationSimulation.hpp"
#include "World.hpp"
#include "radfiled3d/storage/radiation_field_store.hpp"
#include <stdexcept>

using namespace RadiationSimulation;

std::shared_ptr<Geant4::RadiationSimulationHandler> RadiationSimulator::handler;
bool RadiationSimulator::bIsBusy = false;


std::shared_ptr<Geant4::RadiationSimulationHandler> RadiationSimulator::initialize(const int cpu_count)
{
	World::instance = std::make_shared<World>();
	World::world_info = std::make_unique<WorldInfo>("Air", glm::vec3(1.f));
	RadiationSimulator::handler = std::make_shared<Geant4::RadiationSimulationHandler>(cpu_count);
	if (!RadiationSimulator::handler->initialize())
		throw std::runtime_error("Handler initialization failed!");
	return RadiationSimulator::handler;
}

void RadiationSimulator::deinitialize()
{
	if (RadiationSimulator::handler)
		RadiationSimulator::handler->deinitialize();
}

std::shared_ptr<radfiled3d::IRadiationField> RadiationSimulator::simulate_radiation_field(size_t n_particles, radfiled3d::GridTracerAlgorithm tracing_algorithm, bool path_length_weighting)
{
	RadiationSimulator::handler->finalize();
	RadiationSimulator::bIsBusy = true;
	auto field = RadiationSimulator::handler->simulate_radiation_field(n_particles, tracing_algorithm, path_length_weighting);
	RadiationSimulator::bIsBusy = false;
	return field;
}

void RadiationSimulation::RadiationSimulator::add_callback_every_n_particles(std::function<void(std::shared_ptr<radfiled3d::IRadiationField>, size_t)> callback, size_t n_particles)
{
	RadiationSimulator::handler->add_callback_every_n_particles(callback, n_particles);
}

void RadiationSimulation::RadiationSimulator::display_gui()
{
	RadiationSimulator::handler->finalize();
	RadiationSimulator::handler->display_gui();
}

void RadiationSimulation::RadiationSimulator::add_geometry(const std::vector<std::shared_ptr<Geometry::Mesh>>& meshes)
{
	RadiationSimulator::handler->add_geometry(meshes);
	World::Get()->set_geometries(meshes);
}

void RadiationSimulation::RadiationSimulator::add_geometry(std::shared_ptr<Geometry::Mesh> mesh)
{
	RadiationSimulator::add_geometry(std::vector<std::shared_ptr<Geometry::Mesh>>({ mesh }));
}

void RadiationSimulation::RadiationSimulator::add_radiation_source(std::shared_ptr<RadiationSource> source)
{
	World::Get()->radiation_source = source;
}

void RadiationSimulation::RadiationSimulator::set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution, float statistical_error_threshold, float statistical_error_enforcement_ratio, const glm::uvec2& angular_resolution, uint32_t directional_lobes)
{
	RadiationSimulator::handler->set_radiation_field_resolution(radiation_field_dimensions, radiation_field_voxel_dimensions, radiation_field_max_energy, energy_resolution, statistical_error_threshold, statistical_error_enforcement_ratio, angular_resolution, directional_lobes);
}

void RadiationSimulation::RadiationSimulator::set_world_info(std::unique_ptr<RadiationSimulation::WorldInfo> info)
{
	World::world_info = std::move(info);
}

void RadiationSimulation::RadiationSimulator::set_random_seed(uint64_t seed)
{
	RadiationSimulator::handler->set_random_seed(seed);
}
