#include "Geant4/G4SceneConstructor.hpp"
#include "Geometry.hpp"
#include "RadiationSource.hpp"
#include <G4NistManager.hh>
#include <G4ThreeVector.hh>
#include <G4PVPlacement.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include "Geant4/G4World.hpp"
#include "Geant4/G4RadiationFieldDetector.hpp"
#include <stdexcept>
#include <functional>
#include <glm/glm.hpp>
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
		// NON-OWNING: G4Mesh is a G4TessellatedSolid, owned by G4SolidStore — see the note in ::Construct.
		this->g4meshes.push_back(
			std::shared_ptr<G4Mesh>(new G4Mesh(m, length_unit_in_meshes), [](G4Mesh*) {})
		);
	}

	World::Get()->set_geometries(meshes);
}

G4VPhysicalVolume* G4SceneConstructor::Construct()
{
	//auto world_info = World::get_world_info();

	// Size the world so it at least matches the configured field (World::dimensions) but grows to enclose
	// the geometry and the source when either exceeds it: the geometry may be larger than the scored field
	// and the source can sit outside it, yet every primary vertex and all geometry must lie inside the
	// world. A snug fit also avoids tracking through an oversized air volume.
	this->world_dim = G4ThreeVector(World::get_world_info()->dimensions.x * m, World::get_world_info()->dimensions.y * m, World::get_world_info()->dimensions.z * m);

	// Start from the configured field size (WorldDim is the full extent, so its half is the floor). The
	// geometry and source reaches carry a small margin so neither sits exactly on the world boundary, but
	// the floor itself is not inflated: when nothing exceeds WorldDim, the world equals WorldDim.
	glm::vec3 half_extent = glm::vec3(this->world_dim.x(), this->world_dim.y(), this->world_dim.z()) / 2.f;

	const float margin = 1.05f;
	std::function<void(const std::shared_ptr<G4Mesh>&, const glm::vec3&)> grow_to_mesh =
		[&](const std::shared_ptr<G4Mesh>& gm, const glm::vec3& parent_pos) {
			const glm::vec3 pos = parent_pos + gm->getMesh()->getPosition();   // world position, G4 units
			const auto& bb = gm->getBoundingBox();                             // local AABB, G4 units
			half_extent = glm::max(half_extent, margin * glm::abs(pos + bb.first));
			half_extent = glm::max(half_extent, margin * glm::abs(pos + bb.second));
			for (const auto& child : gm->getChildren())
				grow_to_mesh(child, pos);
		};
	for (const auto& gm : this->g4meshes)
		grow_to_mesh(gm, glm::vec3(0.f));

	if (World::Get()->get_radiation_source() != nullptr) {
		const glm::vec3 source_pos = World::Get()->get_radiation_source()->getLocation() * static_cast<float>(m);
		half_extent = glm::max(half_extent, margin * glm::abs(source_pos));
	}

	this->max_world_extend = half_extent * 2.f;

	G4cout << "World volume full extent: "
	       << max_world_extend.x / m << " x " << max_world_extend.y / m << " x " << max_world_extend.z / m
	       << " m (recorded field: "
	       << world_dim.x() / m << " x " << world_dim.y() / m << " x " << world_dim.z() / m << " m)" << G4endl;

	G4Box* worldBox = new G4Box("World", this->max_world_extend.x / 2.0, this->max_world_extend.y / 2.0, this->max_world_extend.z / 2.0);
	this->world_material = MaterialSolver::get_material(World::get_world_info()->material);
	if (this->world_material == NULL)
		throw std::runtime_error("World Material could not be loaded!");

	G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, world_material, "World");
	// Geometry is placed directly in the (large) world, not in a box the size of the scored field. Flux is
	// scored per step from the voxel position (bounds-checked in the detector), so the recorded field
	// (World::dimensions) may be SMALLER than the geometry: the mesh simply extends beyond the field into
	// the surrounding world. A physical tracker box the size of the field would instead force the field to
	// enclose the whole mesh (the mesh would otherwise protrude its mother volume).
	G4VPhysicalVolume* worldPhys = new G4PVPlacement(
		0,
		G4ThreeVector(0, 0, 0),
		worldLog,
		"World",
		nullptr,
		false,
		0
	);

	// Geant4's stores own all geometry (G4SolidStore / G4LogicalVolumeStore / the material table delete these
	// at teardown). Hold them as NON-OWNING shared_ptr (no-op deleter) so only Geant4 deletes them.
	G4World::initialize(
		std::shared_ptr<G4Box>(worldBox, [](G4Box*) {}),
		std::shared_ptr<G4Material>(world_material, [](G4Material*) {}),
		std::shared_ptr<G4LogicalVolume>(worldLog, [](G4LogicalVolume*) {})
	);

	for (auto& m : g4meshes) {
		this->place_mesh(m, worldLog);
	}

	if (World::Get()->get_radiation_field_detector()) {
		World::Get()->get_radiation_field_detector()->SetUp();
	}

	return worldPhys;
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
	G4Element* S = MaterialSolver::nist_man->FindOrBuildElement("S");
	G4Element* Cl = MaterialSolver::nist_man->FindOrBuildElement("Cl");
	G4Element* Fe = MaterialSolver::nist_man->FindOrBuildElement("Fe");
	G4Element* Cs = MaterialSolver::nist_man->FindOrBuildElement("Cs");
	G4Element* I = MaterialSolver::nist_man->FindOrBuildElement("I");

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
		}
	}
	MaterialSolver::custom_materials["CarbonFiber"] = carbonFiber;


	// Aproximate composition from ICRU Report 46, ICRP Publ. 110
	// from: https://gitlab.cern.ch/geant4/geant4/-/blob/master/examples/advanced/ICRP110_HumanPhantoms/src/ICRP110PhantomMaterial_Male.cc
	G4cout << "Creating Spongiosa material as it is not present in this Geant4 version" << G4endl;
	auto spongiosa = new G4Material("ThoraticSpongiosa", 1.074 * g / cm3, 11);
	spongiosa->AddElement(H, 0.099);
	spongiosa->AddElement(C, 0.376);
	spongiosa->AddElement(N, 0.027);
	spongiosa->AddElement(O, 0.459);
	spongiosa->AddElement(Na, 0.001);
	spongiosa->AddElement(P, 0.012);
	spongiosa->AddElement(S, 0.002);
	spongiosa->AddElement(Cl, 0.002);
	spongiosa->AddElement(K, 0.001);
	spongiosa->AddElement(Ca, 0.020);
	spongiosa->AddElement(Fe, 0.001);
	MaterialSolver::custom_materials["ThoraticSpongiosa"] = spongiosa;

	// From https://apc.u-paris.fr/~franco/g4doxy4.10/html/class_materials.html
	G4Material* CsI = new G4Material("CsI", 4.51 * g / cm3, 2);
	CsI->AddElement(Cs, 0.5);
	CsI->AddElement(I, 0.5);
	MaterialSolver::custom_materials["CsI"] = CsI;

	G4Material* PU_foam = new G4Material("PU_foam", 0.05 * g / cm3, 4);
	PU_foam->AddElement(C, 0.48);
	PU_foam->AddElement(H, 0.36);
	PU_foam->AddElement(O, 0.10);
	PU_foam->AddElement(N, 0.06);
	MaterialSolver::custom_materials["PU_foam"] = PU_foam;
}

G4Material* MaterialSolver::get_material(const G4String& name)
{
	G4String material_name = name;
	G4StrUtil::to_upper(material_name);   // G4String::toUpper() is deprecated
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
