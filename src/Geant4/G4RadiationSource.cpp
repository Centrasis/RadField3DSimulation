#include "Geant4//G4RadiationSource.hpp"
#include <G4ParticleTable.hh>
#include "RadiationSource.hpp"
#include "G4SystemOfUnits.hh"
#include <G4ios.hh>
#include <G4Event.hh>

using namespace RadiationSimulation;


G4RadiationSource::G4RadiationSource(std::shared_ptr<RadiationSource> source, int fluence_per_run)
	: source(source),
	  particle_gun(fluence_per_run)
{
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	G4ParticleDefinition* particle = particleTable->FindParticle(source->getParticleName());

	this->particle_gun.SetParticleDefinition(particle);
}

void G4RadiationSource::GeneratePrimaries(G4Event* evt)
{
	this->particle_gun.SetParticlePosition(G4ThreeVector(source->getLocation().x * m, source->getLocation().y * m, source->getLocation().z * m));
	glm::vec3 direction = this->source->drawRayDirection();
	this->particle_gun.SetParticleMomentumDirection(G4ThreeVector(direction.x, direction.y, direction.z));
	double energy = this->source->drawEnergy_eV();
	this->particle_gun.SetParticleEnergy(energy * eV);
	this->particle_gun.GeneratePrimaryVertex(evt);
}
