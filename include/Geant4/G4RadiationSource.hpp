#pragma once
#define _USE_MATH_DEFINES
#include <G4VUserPrimaryGeneratorAction.hh>
#include <memory>
#include <G4VUserActionInitialization.hh>
#include <G4ParticleGun.hh>


namespace RadiationSimulation {
	class RadiationSource;

	class G4RadiationSource : public G4VUserPrimaryGeneratorAction {
	private:
		G4ParticleGun particle_gun;
		const std::shared_ptr<RadiationSource> source;
	public:
		G4RadiationSource(std::shared_ptr<RadiationSource> source, int fluence_per_run = 1);
		void GeneratePrimaries(G4Event* evt);
	};
}