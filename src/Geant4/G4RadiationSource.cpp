#include "Geant4//G4RadiationSource.hpp"
#include <G4ParticleTable.hh>
#include "RadiationSource.hpp"
#include "G4SystemOfUnits.hh"
#include <G4ios.hh>
#include <G4Event.hh>
#include <Randomize.hh>

using namespace RadiationSimulation;
using namespace RadiationSimulation::Geant4;


Geant4::RadiationSource::RadiationSource(std::shared_ptr<RadiationSimulation::RadiationSource> source, int fluence_per_run)
	: source(source),
	  particle_gun(fluence_per_run)
{
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	G4ParticleDefinition* particle = particleTable->FindParticle(source->getParticleName());

	this->particle_gun.SetParticleDefinition(particle);
}

void Geant4::RadiationSource::GeneratePrimaries(G4Event* evt)
{
	this->particle_gun.SetParticlePosition(G4ThreeVector(source->getLocation().x * m, source->getLocation().y * m, source->getLocation().z * m));
	// The worker's engine is re-seeded by the run manager for every event, so the primary, like the rest of the event,
	// is determined by the run's seed and the event number.
	const RadiationSimulation::UniformRandom uniform = [] { return G4UniformRand(); };
	glm::vec3 direction = this->source->drawRayDirection(uniform);
	this->particle_gun.SetParticleMomentumDirection(G4ThreeVector(direction.x, direction.y, direction.z));
	double energy = this->source->drawEnergy_eV(uniform);
	this->particle_gun.SetParticleEnergy(energy * eV);
	this->particle_gun.GeneratePrimaryVertex(evt);
}
