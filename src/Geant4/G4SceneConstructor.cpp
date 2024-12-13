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
		throw new std::runtime_error("Tracking volume is too big for the max world size!");
	}
	
	this->world_dim = G4ThreeVector(World::get_world_info()->dimensions.x * m, World::get_world_info()->dimensions.y * m, World::get_world_info()->dimensions.z * m);
	G4Box* worldBox = new G4Box("World", this->max_world_extend.x / 2.0, this->max_world_extend.y / 2.0, this->max_world_extend.z / 2.0);
	this->world_material = MaterialSolver::get_material(World::get_world_info()->material);
	if (this->world_material == NULL)
		throw new std::runtime_error("World Material could not be loaded!");

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

	G4Material* Polyamide = new G4Material("Polyamide", 1.14 * g / cm3, 4);
	Polyamide->AddElement(H, 11);
	Polyamide->AddElement(C, 6);
	Polyamide->AddElement(N, 1);
	Polyamide->AddElement(O, 1);

	MaterialSolver::custom_materials["Polyamide"] = Polyamide;
}

G4Material* MaterialSolver::get_material(const G4String& name)
{
	G4String material_name = name;
	material_name.toUpper();
	if (material_name.substr(0, 2) != "G4") {
		material_name = G4String("G4_") + material_name;
	}
	G4Material* mat = NULL;
	try {
		mat = MaterialSolver::get_nist_man()->FindOrBuildMaterial(material_name);
	}
	catch (std::runtime_error&) {
	}
	if (mat == NULL) {
		try {
			mat = MaterialSolver::get_nist_man()->FindOrBuildMaterial(name);
		}
		catch (std::runtime_error&) {
		}
	}
	if (mat == NULL) {
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
