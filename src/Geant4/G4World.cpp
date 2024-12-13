#include "Geant4/G4World.hpp"


using namespace RadiationSimulation;

RadiationSimulation::G4World::G4World(std::shared_ptr<World> raw_world)
	: World(),
	  raw_world(raw_world)
{
}

void RadiationSimulation::G4World::initialize(std::shared_ptr<G4Box> box, std::shared_ptr<G4Material> material, std::shared_ptr<G4LogicalVolume> volume)
{
	World::instance = std::make_shared<G4World>(World::instance);
	static_cast<G4World*>(World::instance.get())->box = box;
	static_cast<G4World*>(World::instance.get())->material = material;
	static_cast<G4World*>(World::instance.get())->volume = volume;
}

std::shared_ptr<G4World> RadiationSimulation::G4World::Get()
{
	return std::dynamic_pointer_cast<G4World>(World::Get());
}
