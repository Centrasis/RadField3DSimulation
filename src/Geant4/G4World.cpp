#include "Geant4/G4World.hpp"


using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Geant4;

RadiationSimulation::Geant4::World::World(std::shared_ptr<RadiationSimulation::World> raw_world)
	: RadiationSimulation::World(),
	  raw_world(raw_world)
{
}

void RadiationSimulation::Geant4::World::initialize(std::shared_ptr<G4Box> box, std::shared_ptr<G4Material> material, std::shared_ptr<G4LogicalVolume> volume)
{
	RadiationSimulation::World::instance = std::make_shared<Geant4::World>(RadiationSimulation::World::instance);
	static_cast<Geant4::World*>(RadiationSimulation::World::instance.get())->box = box;
	static_cast<Geant4::World*>(RadiationSimulation::World::instance.get())->material = material;
	static_cast<Geant4::World*>(RadiationSimulation::World::instance.get())->volume = volume;
}

std::shared_ptr<Geant4::World> RadiationSimulation::Geant4::World::Get()
{
	return std::dynamic_pointer_cast<Geant4::World>(RadiationSimulation::World::Get());
}
