#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"


// Electromagnetic-only physics for keV-range X-ray transport (photons + secondary electrons). Option4 is the
// most accurate low-energy EM constructor (Livermore/Penelope-based models). No hadronic/neutron physics:
// diagnostic-energy photons stay far below the ~MeV photonuclear threshold, so hadronic constructors never
// fire — they would only add startup cost and load neutron cross-section data that is never used.
class MedicalPhysicsList : public G4VModularPhysicsList {
public:
    MedicalPhysicsList() {
        RegisterPhysics(new G4EmStandardPhysics_option4());
    }
};
