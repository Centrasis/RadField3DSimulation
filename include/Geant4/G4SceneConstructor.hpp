#pragma once
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "Geant4/G4Geometry.hpp"
#include <memory>
#include <vector>
#include <G4SystemOfUnits.hh>


class G4NistManager;
class G4Material;

namespace RadiationSimulation {
	class G4SceneConstructor : public G4VUserDetectorConstruction
	{
	protected:
		std::vector<std::shared_ptr<G4Mesh>> g4meshes;
		double length_unit_in_meshes = m;
		double max_step_width = 1e-7 * mm;
		glm::vec3 max_world_extend = glm::vec3(1e+3 * m);
		G4Material* world_material = NULL;
		G4ThreeVector world_dim;
		void place_mesh(std::shared_ptr<G4Mesh> mesh, G4LogicalVolume* parent);
	public:
		G4SceneConstructor(const std::vector<std::shared_ptr<Mesh>>& meshes);
		virtual ~G4SceneConstructor() {
			G4cout << "G4SceneConstructor destroyed" << G4endl;
		}
		G4VPhysicalVolume* Construct();
		virtual void ConstructSDandField() override;
	};

	class MaterialSolver {
	protected:
		static G4NistManager* nist_man;
		static void init_custom_materials();
		static std::map<G4String, G4Material*> custom_materials;
	public:
		static G4Material* get_material(const G4String& name);
		static G4NistManager* get_nist_man();
	};
}