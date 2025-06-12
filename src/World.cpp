#include "World.hpp"
#include <stdexcept>
#include "Geometry.hpp"

using namespace RadiationSimulation;

std::shared_ptr<World>		World::instance = std::shared_ptr<World>(NULL);
std::unique_ptr<WorldInfo>	World::world_info = std::unique_ptr<WorldInfo>();
	
World::World()
{
}

void RadiationSimulation::World::initialize()
{
	if (World::instance.get() != NULL)
		throw std::runtime_error("Can't reinitialize World!");
	World::instance = std::make_shared<World>();
}

void RadiationSimulation::World::set_geometries(const std::vector<std::shared_ptr<Mesh>>& geoms)
{
	this->geometries = geoms;
	for (auto mesh : this->geometries)
		if (mesh->isPatient())
			this->patient = mesh;
}
