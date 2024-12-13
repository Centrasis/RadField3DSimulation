#include "Geant4/G4Geometry.hpp"
#include <G4QuadrangularFacet.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <stdexcept>


using namespace RadiationSimulation;


G4Mesh::G4Mesh(std::shared_ptr<Mesh> mesh, double length_unit)
	:	G4TessellatedSolid(mesh->getName()),
	    length_unit(length_unit),
		mesh(mesh)
{
	const glm::vec3 scale = mesh->getScale();
	for (size_t i = 0; i < mesh->vertices.size(); i++) {
		mesh->vertices[i].x *= length_unit * scale.x;
		mesh->vertices[i].y *= length_unit * scale.y;
		mesh->vertices[i].z *= length_unit * scale.z;
	}

	mesh->position *= length_unit;

	for (auto f: mesh->getFaces()) {
		G4VFacet* face = NULL;
		const TriFace* tf = NULL;
		const QuadFace* qf = NULL;

		switch (f->getType())
		{
		case FaceType::Tri:
			tf = static_cast<const TriFace*>(f);
			face = new G4TriangularFacet(
				G4ThreeVector(
					mesh->vertices[tf->getIndices().r].x,
					mesh->vertices[tf->getIndices().r].y,
					mesh->vertices[tf->getIndices().r].z
				),
				G4ThreeVector(
					mesh->vertices[tf->getIndices().g].x,
					mesh->vertices[tf->getIndices().g].y,
					mesh->vertices[tf->getIndices().g].z
				),
				G4ThreeVector(
					mesh->vertices[tf->getIndices().b].x,
					mesh->vertices[tf->getIndices().b].y,
					mesh->vertices[tf->getIndices().b].z
				),
				G4FacetVertexType::ABSOLUTE
			);
			break;
		case FaceType::Quad:
			qf = static_cast<const QuadFace*>(f);
			face = new G4QuadrangularFacet(
				G4ThreeVector(
					mesh->vertices[qf->getIndices().r].x,
					mesh->vertices[qf->getIndices().r].y,
					mesh->vertices[qf->getIndices().r].z
				),
				G4ThreeVector(
					mesh->vertices[qf->getIndices().g].x,
					mesh->vertices[qf->getIndices().g].y,
					mesh->vertices[qf->getIndices().g].z
				),
				G4ThreeVector(
					mesh->vertices[qf->getIndices().b].x,
					mesh->vertices[qf->getIndices().b].y,
					mesh->vertices[qf->getIndices().b].z
				),
				G4ThreeVector(
					mesh->vertices[qf->getIndices().a].x,
					mesh->vertices[qf->getIndices().a].y,
					mesh->vertices[qf->getIndices().a].z
				),
				G4FacetVertexType::ABSOLUTE
			);
			break;
		default:
			throw std::runtime_error("Unknown Face type!");
			break;
		}
		
		this->AddFacet(face);
	}
	this->SetSolidClosed(true);

	G4ThreeVector min;
	G4ThreeVector max;
	this->BoundingLimits(min, max);
	this->mesh->bounding_box = {
		glm::vec3(min.getX(), min.getY(), min.getZ()),
		glm::vec3(max.getX(), max.getY(), max.getZ())
	};

	this->rotation.rotateX(mesh->getRotation().x);
	this->rotation.rotateY(mesh->getRotation().y);
	this->rotation.rotateZ(mesh->getRotation().z);

	this->position = G4ThreeVector(mesh->position.x, mesh->position.y, mesh->position.z);

	for (auto& child : mesh->children) {
		this->children.push_back(std::make_shared<G4Mesh>(child, length_unit));
	}
}

void G4Mesh::place(G4LogicalVolume* parent)
{
	this->physical = std::make_shared<G4PVPlacement>(
		&this->rotation,
		this->position,
		this->getVolume().get(),
		this->GetName(),
		parent,
		false,
		0,
		true
	);
}

void G4Mesh::setMaterial(G4Material* material)
{
	G4String type = "Tracker";
	this->volume = std::make_shared<G4LogicalVolume>(this, material, type);
}

std::shared_ptr<G4LogicalVolume> G4Mesh::getVolume() {
	if (!this->volume.get()) {
		G4String type = "Tracker";
		this->volume = std::make_shared<G4LogicalVolume>(this, (G4Material*)NULL, type);
	}
	return this->volume;
}

const std::pair<glm::vec3, glm::vec3>& G4Mesh::getBoundingBox() const
{
	return this->mesh->bounding_box;
}

