#include "G4VModularPhysicsList.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4RayleighScattering.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "QGSP_BIC.hh"


class MedicalPhysicsList : public QGSP_BIC_HP { // QGSP_BIC
protected:
    G4EmStandardPhysics_option4* emPhysicsList;
public:
    MedicalPhysicsList()
        :   emPhysicsList(new G4EmStandardPhysics_option4()),
            QGSP_BIC_HP() // QGSP_BIC
    {
        size_t i = 0;
        bool found_g4_phys_list = false;
		bool found_extra_phys = false;
		bool found_low_ep_phys = false;
        while(true) {
            try {
                auto p = this->GetPhysics(i);
                if (p == NULL || p == nullptr) break;
                if (std::string("G4EmStandard_opt4") == p->GetPhysicsName()) {
                    found_g4_phys_list = true;
                }
				if (std::string("G4EmExtraPhysics") == p->GetPhysicsName()) {
					found_extra_phys = true;
				}
                if (std::string("G4EmLowEPPhysics") == p->GetPhysicsName()) {
                    found_low_ep_phys = true;
                }
				if (found_g4_phys_list && found_extra_phys && found_low_ep_phys) break;
                i++;
            } catch (std::exception& e) {
				break;
			}
		}
        // Register electromagnetic physics
        if (!found_g4_phys_list) {
            RegisterPhysics(this->emPhysicsList);
        }
    }
};