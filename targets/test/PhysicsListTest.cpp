#include "Geant4/G4PhysicsList.hpp"
#include <gtest/gtest.h>
#include <G4Box.hh>
#include <G4EmParameters.hh>
#include <G4Gamma.hh>
#include <G4LogicalVolume.hh>
#include <G4NistManager.hh>
#include <G4PVPlacement.hh>
#include <G4ProcessManager.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4ProductionCutsTable.hh>
#include <G4MaterialCutsCouple.hh>
#include <G4ProductionCuts.hh>
#include <G4Version.hh>
#include <set>
#include <sstream>
#include <string>
#include <vector>


namespace {
	class AirWorld : public G4VUserDetectorConstruction {
	public:
		G4VPhysicalVolume* Construct() override
		{
			auto* box = new G4Box("World", 1. * m, 1. * m, 1. * m);
			auto* volume = new G4LogicalVolume(box, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "World");
			// the heaviest materials of the scenes: shields and the detector
			int i = 0;
			for (const char* name : { "G4_WATER", "G4_BONE_CORTICAL_ICRP", "G4_I", "G4_W", "G4_Pb" }) {
				auto* part = new G4LogicalVolume(new G4Box(name, 5. * cm, 5. * cm, 5. * cm), G4NistManager::Instance()->FindOrBuildMaterial(name), name);
				new G4PVPlacement(nullptr, G4ThreeVector(-0.8 * m + (i++) * 0.3 * m, 0., 0.), part, name, volume, false, 0);
			}
			return new G4PVPlacement(nullptr, G4ThreeVector(), volume, "World", nullptr, false, 0);
		}
	};

	// "x.y.z" from the release tag in G4Version, e.g. "$Name: geant4-11-04-patch-02 $" (no patch part in .0 releases)
	std::string version_from_release_tag()
	{
		const std::string tag = G4Version.substr(G4Version.find("geant4-") + 7);
		std::vector<std::string> parts;
		std::stringstream stream(tag.substr(0, tag.find_first_of(" $")));
		for (std::string part; std::getline(stream, part, '-');)
			parts.push_back(part);
		const std::string patch = (parts.size() >= 4 && parts[2] == "patch") ? parts[3] : "0";
		return std::to_string(std::stoi(parts[0])) + "." + std::to_string(std::stoi(parts[1])) + "." + std::to_string(std::stoi(patch));
	}
}

TEST(MedicalPhysicsList, NameHoldsTheGeant4VersionAndTheEnabledPhotonProcess) {
	// runs before any run manager exists: the EM parameters lock once a run was set up
	G4EmParameters::Instance()->SetGeneralProcessActive(true);
	EXPECT_EQ(RadiationSimulation::Geant4::MedicalPhysicsList::getName(), "G4-" + version_from_release_tag() + ":MedicalGeneralPhotonProcess");

	G4EmParameters::Instance()->SetGeneralProcessActive(false);
	EXPECT_EQ(RadiationSimulation::Geant4::MedicalPhysicsList::getName(), "G4-" + version_from_release_tag());
}

TEST(MedicalPhysicsList, TransportsEveryPhotonAboveFiveKeVWithTheGeneralGammaProcess) {
	// the state 11.0 starts in; option4 of 11.4 would switch it on by itself
	G4EmParameters::Instance()->SetGeneralProcessActive(false);

	auto* run_manager = new G4RunManager();
	run_manager->SetUserInitialization(new AirWorld());
	auto* physics = new RadiationSimulation::Geant4::MedicalPhysicsList();
	run_manager->SetUserInitialization(physics);
	EXPECT_TRUE(G4EmParameters::Instance()->GeneralProcessActive());
	// shields' L lines (8-15 keV) lie below lead's and tungsten's gamma threshold at the 0.1 mm cut
	EXPECT_TRUE(G4EmParameters::Instance()->DeexcitationIgnoreCut());
	EXPECT_TRUE(G4EmParameters::Instance()->Fluo());
	// electrons are stopped where they are created (fast, photons unaffected)
	EXPECT_GE(G4EmParameters::Instance()->LowestElectronEnergy(), 1. * MeV);
	// the simulator's default cut, as RadiationSimulationHandler sets it
	physics->SetDefaultCutValue(0.1 * mm);
	run_manager->Initialize();
	run_manager->RunInitialization();

	// no secondary photon above 5 keV is cut in any material; electrons keep the 0.1 mm cut
	const G4ProductionCutsTable* cuts = G4ProductionCutsTable::GetProductionCutsTable();
	for (size_t i = 0; i < cuts->GetTableSize(); i++) {
		const G4String& material = cuts->GetMaterialCutsCouple(i)->GetMaterial()->GetName();
		EXPECT_LE((*cuts->GetEnergyCutsVector(idxG4GammaCut))[i], 5. * keV) << material;
		EXPECT_DOUBLE_EQ(cuts->GetMaterialCutsCouple(i)->GetProductionCuts()->GetProductionCut("e-"), 0.1 * mm) << material;
	}

	std::set<std::string> processes;
	const G4ProcessVector* list = G4Gamma::Definition()->GetProcessManager()->GetProcessList();
	for (size_t i = 0; i < list->size(); i++)
		processes.insert((*list)[i]->GetProcessName());

	EXPECT_EQ(processes.count("GammaGeneralProc"), 1u);
	// the four photon processes are merged into it instead of being registered separately
	for (const char* separate : { "phot", "compt", "Rayl", "conv" })
		EXPECT_EQ(processes.count(separate), 0u) << separate;
	delete run_manager;
}
