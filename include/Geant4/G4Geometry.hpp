#pragma once
#include <G4TessellatedSolid.hh>
#include "Geometry.hpp"
#include <memory>
#include <G4PolyhedronArbitrary.hh>

class G4PVPlacement;
class G4LogicalVolume;
class G4Material;


namespace RadiationSimulation {
	class G4Mesh : public G4TessellatedSolid {
	protected:
		std::shared_ptr<Mesh> mesh;
		std::vector<std::shared_ptr<G4Mesh>> children;
		const double length_unit;
		G4ThreeVector position = G4ThreeVector(0);
		G4RotationMatrix rotation;
		std::shared_ptr<G4LogicalVolume> volume;
		std::shared_ptr<G4PVPlacement> physical;
	public:
		G4Mesh(std::shared_ptr<Mesh> mesh, double length_unit = 1.0);
		~G4Mesh() {
			G4cout << "G4Mesh destroyed" << G4endl;
		}
		std::shared_ptr<Mesh> getMesh() const { return this->mesh; }
		void setMaterial(G4Material* material);
		void place(G4LogicalVolume* parent);
		std::shared_ptr<G4LogicalVolume> getVolume();
		const std::pair<glm::vec3, glm::vec3>& getBoundingBox() const;
		const G4RotationMatrix& getRotation() const { return this->rotation; }
		double getLengthUnit() const { return this->length_unit; }
		inline const std::vector<std::shared_ptr<G4Mesh>>& getChildren() const { return this->children; }
	};
}