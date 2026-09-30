#pragma once
#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmParameters.hh"
#include "G4Version.hh"
#include "G4SystemOfUnits.hh"
#include <string>


// Electromagnetic-only physics for keV-range X-ray transport (photons + secondary electrons). Option4 is the
// most accurate low-energy EM constructor (Livermore/Penelope-based models). No hadronic/neutron physics:
// diagnostic-energy photons stay far below the ~MeV photonuclear threshold, so hadronic constructors never
// fire — they would only add startup cost and load neutron cross-section data that is never used.
namespace RadiationSimulation::Geant4 {
    class MedicalPhysicsList : public G4VModularPhysicsList {
    public:
        MedicalPhysicsList() {
            RegisterPhysics(new G4EmStandardPhysics_option4());
            // One combined photon process (photoelectric, Compton, Rayleigh, conversion): same models, one cross-section
            // lookup per photon step. Set after option4, whose constructor resets the EM parameters; 11.1+ enables it
            // by default, 11.0 does not.
            G4EmParameters::Instance()->SetGeneralProcessActive(true);
            // Fluorescence regardless of the production cuts: at 0.1 mm the gamma threshold is 29 keV in lead and 36 keV
            // in tungsten, which would drop their L lines (8-15 keV) that shields emit into the room.
            G4EmParameters::Instance()->SetDeexcitationIgnoreCut(true);
            // Electrons and positrons are stopped where they are created and deposit their energy there. At diagnostic
            // energies they travel < 0.3 mm in tissue and only photons are scored; their bremsstrahlung and impact
            // fluorescence are ~0.1 % of the energy. Full electron transport took ~92 % of the run time (Angio scene,
            // 11 M photons: 921 s vs 78 s, scored field within 0.004 %).
            G4EmParameters::Instance()->SetLowestElectronEnergy(1. * MeV);
        }

        /** Keeps the default cut for electrons and positrons, but produces secondary photons (bremsstrahlung) down to
        * the 990 eV lower edge in every material: 1 µm gives at most 1.7 keV even in tungsten and lead, where 0.1 mm
        * would drop photons below 36 and 29 keV. No photon above 5 keV is cut.
        */
        void SetCuts() override {
            G4VUserPhysicsList::SetCuts();
            SetCutValue(photon_range_cut, "gamma");
        }

        static constexpr double photon_range_cut = 0.001 * mm;

        /** Name of the physics for the field metadata: the Geant4 version the module is built against and the enabled
        * photon process, e.g. "G4-11.4.2:MedicalGeneralPhotonProcess". Runs with a different name are never appended to a field.
        */
        static std::string getName() {
            // G4VERSION_NUMBER is major * 100 + minor * 10 + patch, e.g. 1142
            std::string name = "G4-" + std::to_string(G4VERSION_NUMBER / 100) + "." + std::to_string((G4VERSION_NUMBER / 10) % 10) + "." + std::to_string(G4VERSION_NUMBER % 10);
            if (G4EmParameters::Instance()->GeneralProcessActive())
                name += ":MedicalGeneralPhotonProcess";
            return name;
        }
    };
}
