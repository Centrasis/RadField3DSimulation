#pragma once
#define _USE_MATH_DEFINES
#include <G4VUserPrimaryGeneratorAction.hh>
#include <memory>
#include <G4VUserActionInitialization.hh>
#include <G4ParticleGun.hh>


namespace RadiationSimulation {
	class RadiationSource;
}

namespace RadiationSimulation::Geant4 {
	class RadiationSource : public G4VUserPrimaryGeneratorAction {
	private:
		G4ParticleGun particle_gun;
		const std::shared_ptr<RadiationSimulation::RadiationSource> source;
	public:
		RadiationSource(std::shared_ptr<RadiationSimulation::RadiationSource> source, int fluence_per_run = 1);
		void GeneratePrimaries(G4Event* evt);

		virtual ~RadiationSource() {
			G4cout << "RadiationSource destroyed" << G4endl;
		}
	};
}