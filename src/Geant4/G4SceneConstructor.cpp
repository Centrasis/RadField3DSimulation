#include "Geant4/G4SceneConstructor.hpp"
#include "Geometry.hpp"
#include <G4NistManager.hh>
#include <G4ThreeVector.hh>
#include <G4PVPlacement.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include "Geant4/G4World.hpp"
#include "RadiationFieldDetector.hpp"
#include "Geant4/G4RadiationFieldDetector.hpp"
#include <stdexcept>
#include "G4Tubs.hh"


using namespace RadiationSimulation;


G4NistManager* MaterialSolver::nist_man = NULL;

void RadiationSimulation::G4SceneConstructor::place_mesh(std::shared_ptr<G4Mesh> mesh, G4LogicalVolume* parent)
{
	auto m_name = mesh->getMesh()->getMaterialName();
	if (m_name.size() > 0)
		mesh->setMaterial(MaterialSolver::get_material(m_name));
	else
		mesh->setMaterial(this->world_material);
	mesh->place(parent);

	for (auto& child : mesh->getChildren()) {
		this->place_mesh(child, mesh->getVolume().get());
	}
}

G4SceneConstructor::G4SceneConstructor(const std::vector<std::shared_ptr<Mesh>>& meshes)
	: G4VUserDetectorConstruction()
{
	this->world_dim = G4ThreeVector(0, 0, 0);
	for (auto& m : meshes) {
		this->g4meshes.push_back(
			std::make_shared<G4Mesh>(m, length_unit_in_meshes)
		);
	}

	World::Get()->set_geometries(meshes);
}

G4VPhysicalVolume* G4SceneConstructor::Construct()
{
	//auto world_info = World::get_world_info();

	if (this->max_world_extend.x < World::get_world_info()->dimensions.x * m || this->max_world_extend.y < World::get_world_info()->dimensions.y * m || this->max_world_extend.z < World::get_world_info()->dimensions.z * m) {
		throw std::runtime_error("Tracking volume is too big for the max world size!");
	}
	
	this->world_dim = G4ThreeVector(World::get_world_info()->dimensions.x * m, World::get_world_info()->dimensions.y * m, World::get_world_info()->dimensions.z * m);
	G4Box* worldBox = new G4Box("World", this->max_world_extend.x / 2.0, this->max_world_extend.y / 2.0, this->max_world_extend.z / 2.0);
	this->world_material = MaterialSolver::get_material(World::get_world_info()->material);
	if (this->world_material == NULL)
		throw std::runtime_error("World Material could not be loaded!");

	G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, world_material, "World");
	G4Box* tracker_box = new G4Box("Tracker", world_dim.x() / 2.0, world_dim.y() / 2.0, world_dim.z() / 2.0);
	G4LogicalVolume* trackerLog = new G4LogicalVolume(tracker_box, world_material, "Tracker");
	G4VPhysicalVolume* trackerPhys = new G4PVPlacement(
		0,
		G4ThreeVector(0, 0, 0),
		trackerLog,
		"Tracker",
		worldLog,
		false,
		0
	);

	G4World::initialize(
		std::shared_ptr<G4Box>(worldBox),
		std::shared_ptr<G4Material>(world_material),
		std::shared_ptr<G4LogicalVolume>(trackerLog)
	);

	for (auto& m : g4meshes) {
		this->place_mesh(m, trackerLog);
	}

	if (World::Get()->get_radiation_field_detector()) {
		World::Get()->get_radiation_field_detector()->SetUp();
	}

	return trackerPhys;
}

void G4SceneConstructor::ConstructSDandField()
{
	
}

std::map<G4String, G4Material*> RadiationSimulation::MaterialSolver::custom_materials;

void RadiationSimulation::MaterialSolver::init_custom_materials()
{
	G4Element* H = MaterialSolver::nist_man->FindOrBuildElement("H");
	G4Element* C = MaterialSolver::nist_man->FindOrBuildElement("C");
	G4Element* N = MaterialSolver::nist_man->FindOrBuildElement("N");
	G4Element* O = MaterialSolver::nist_man->FindOrBuildElement("O");
	G4Element* P = MaterialSolver::nist_man->FindOrBuildElement("P");
	G4Element* Ca = MaterialSolver::nist_man->FindOrBuildElement("Ca");
	G4Element* Mg = MaterialSolver::nist_man->FindOrBuildElement("Mg");
	G4Element* Na = MaterialSolver::nist_man->FindOrBuildElement("Na");
	G4Element* K = MaterialSolver::nist_man->FindOrBuildElement("K");

	G4Material* Polyamide = new G4Material("Polyamide", 1.14 * g / cm3, 4);
	Polyamide->AddElement(H, 11);
	Polyamide->AddElement(C, 6);
	Polyamide->AddElement(N, 1);
	Polyamide->AddElement(O, 1);

	MaterialSolver::custom_materials["Polyamide"] = Polyamide;

	auto nist = G4NistManager::Instance();
	auto carbonFiber = nist->FindMaterial("G4_CARBON_FIBER");
	if (carbonFiber == nullptr) {
		carbonFiber = G4Material::GetMaterial("G4_CARBON_FIBER", false);

		if (carbonFiber == nullptr) {
			G4cout << "Creating CarbonFiber material as it is not present in this Geant4 version" << G4endl;
			carbonFiber = new G4Material("CarbonFiber", 1.60 * g / cm3, 2);
			G4Material* GraphiteCarbon = new G4Material("GraphiteCarbon", 1.80 * g / cm3, 1);
			GraphiteCarbon->AddElement(C, 1.0);
			G4Material* Epoxy = new G4Material("Epoxy", 1.20 * g / cm3, 3);
			Epoxy->AddElement(C, 0.60);
			Epoxy->AddElement(H, 0.08);
			Epoxy->AddElement(O, 0.32);
			carbonFiber->AddMaterial(GraphiteCarbon, 0.60);
			carbonFiber->AddMaterial(Epoxy, 0.40);

			MaterialSolver::custom_materials["CarbonFiber"] = carbonFiber;
		}
	}


	// Aproximate composition from ICRU Report 46, ICRP Publ. 110
	// from: https://gitlab.cern.ch/geant4/geant4/-/blob/master/examples/advanced/ICRP110_HumanPhantoms/src/ICRP110PhantomMaterial_Male.cc
	G4cout << "Creating Spongiosa material as it is not present in this Geant4 version" << G4endl;
	auto spongiosa = new G4Material("ThoraticSpongiosa", 1.074 * g / cm3, 9);
	spongiosa->AddElement(H, 0.055);
	spongiosa->AddElement(C, 0.185);
	spongiosa->AddElement(N, 0.035);
	spongiosa->AddElement(O, 0.420);
	spongiosa->AddElement(P, 0.065);
	spongiosa->AddElement(Ca, 0.205);
	spongiosa->AddElement(Mg, 0.010);
	spongiosa->AddElement(Na, 0.010);
	spongiosa->AddElement(K, 0.015);

	MaterialSolver::custom_materials["ThoraticSpongiosa"] = spongiosa;
}

G4Material* MaterialSolver::get_material(const G4String& name)
{
	G4String material_name = name;
	material_name.toUpper();
	if (material_name.substr(0, 2) != "G4") {
		material_name = G4String("G4_") + material_name;
	}
	G4Material* mat = nullptr;
	try {
		mat = MaterialSolver::get_nist_man()->FindOrBuildMaterial(material_name);
	}
	catch (std::runtime_error&) {
	}
	if (mat == nullptr) {
		try {
			mat = MaterialSolver::get_nist_man()->FindOrBuildMaterial(name);
		}
		catch (std::runtime_error&) {
		}
	}
	if (mat == nullptr) {
		mat = G4Material::GetMaterial(material_name, false);
	}
	if (mat == nullptr) {
		mat = G4Material::GetMaterial(name, false);
	}

	if (mat == nullptr) {
		auto found = MaterialSolver::custom_materials.find(name);
		if (found == MaterialSolver::custom_materials.end())
			throw std::runtime_error("Material '" + name + "' could not be found!");
		mat = found->second;
	}
	return mat;
}

G4NistManager* MaterialSolver::get_nist_man()
{
	if (MaterialSolver::nist_man == NULL) {
		MaterialSolver::nist_man = G4NistManager::Instance();
		MaterialSolver::init_custom_materials();
	}
	return MaterialSolver::nist_man;
}
