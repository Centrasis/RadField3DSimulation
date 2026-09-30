#include "Geant4/G4Geometry.hpp"
#include <G4QuadrangularFacet.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4SystemOfUnits.hh>
#include <algorithm>
#include <chrono>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>


using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Geant4;


Geant4::Mesh::Mesh(std::shared_ptr<Geometry::Mesh> mesh, double length_unit)
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

	glm::vec3 rot_angles = glm::eulerAngles(mesh->getRotation());

	this->rotation.rotateX(rot_angles.x);
	this->rotation.rotateY(rot_angles.y);
	this->rotation.rotateZ(rot_angles.z);

	this->position = G4ThreeVector(mesh->position.x, mesh->position.y, mesh->position.z);

	for (auto& child : mesh->children) {
		// NON-OWNING: Geant4::Mesh is a G4TessellatedSolid (G4SolidStore-owned) — an owning shared_ptr double-frees.
		this->children.push_back(std::shared_ptr<Geant4::Mesh>(new Geant4::Mesh(child, length_unit), [](Geant4::Mesh*) {}));
	}
}

void Geant4::Mesh::place(G4LogicalVolume* parent)
{
	// NON-OWNING: G4PhysicalVolumeStore owns the placement — an owning shared_ptr double-frees at teardown.
	this->physical = std::shared_ptr<G4PVPlacement>(
		new G4PVPlacement(
			&this->rotation,
			this->position,
			this->getVolume().get(),
			this->GetName(),
			parent,
			false,
			0,
			true
		),
		[](G4PVPlacement*) {}
	);
}

std::vector<const G4Material*> Geant4::Mesh::getMaterials() const
{
	std::vector<const G4Material*> materials;
	if (this->volume && this->volume->GetMaterial() != nullptr)
		materials.push_back(this->volume->GetMaterial());
	for (const auto& child : this->children)
		for (const G4Material* material : child->getMaterials())
			if (std::find(materials.begin(), materials.end(), material) == materials.end())
				materials.push_back(material);
	return materials;
}

void Geant4::Mesh::setMaterial(G4Material* material)
{
	G4String type = "Tracker";
	// NON-OWNING: G4LogicalVolumeStore owns the volume — an owning shared_ptr double-frees at teardown.
	this->volume = std::shared_ptr<G4LogicalVolume>(new G4LogicalVolume(this, material, type), [](G4LogicalVolume*) {});
}

std::shared_ptr<G4LogicalVolume> Geant4::Mesh::getVolume() {
	if (!this->volume.get()) {
		G4String type = "Tracker";
		// NON-OWNING: G4LogicalVolumeStore owns the volume (see setMaterial).
		this->volume = std::shared_ptr<G4LogicalVolume>(new G4LogicalVolume(this, (G4Material*)NULL, type), [](G4LogicalVolume*) {});
	}
	return this->volume;
}

const std::pair<glm::vec3, glm::vec3>& Geant4::Mesh::getBoundingBox() const
{
	return this->mesh->bounding_box;
}

Geant4::PlacedMeshSurface::PlacedMeshSurface(const G4TessellatedSolid& solid, const G4RotationMatrix& rotation, const G4ThreeVector& translation)
	: solid(solid),
	  rotation(rotation),
	  translation(translation)
{
}

size_t Geant4::PlacedMeshSurface::polygon_count() const
{
	return static_cast<size_t>(this->solid.GetNumberOfFacets());
}

size_t Geant4::PlacedMeshSurface::vertex_count(size_t polygon) const
{
	return static_cast<size_t>(this->solid.GetFacet(static_cast<G4int>(polygon))->GetNumberOfVertices());
}

glm::dvec3 Geant4::PlacedMeshSurface::vertex(size_t polygon, size_t vertex) const
{
	const G4ThreeVector world = this->rotation * this->solid.GetFacet(static_cast<G4int>(polygon))->GetVertex(static_cast<G4int>(vertex)) + this->translation;
	return glm::dvec3(world.x(), world.y(), world.z()) / m;
}

namespace {
	struct PlacedMesh {
		const Geant4::Mesh* mesh;
		G4RotationMatrix rotation;
		G4ThreeVector translation;
	};

	void collect_placed_meshes(const G4LogicalVolume& volume, const G4RotationMatrix& rotation, const G4ThreeVector& translation, std::vector<PlacedMesh>& placed)
	{
		for (size_t i = 0; i < volume.GetNoDaughters(); i++) {
			const G4VPhysicalVolume* daughter = volume.GetDaughter(i);
			const G4RotationMatrix daughter_rotation = rotation * daughter->GetObjectRotationValue();
			const G4ThreeVector daughter_translation = rotation * daughter->GetObjectTranslation() + translation;
			const G4LogicalVolume* daughter_volume = daughter->GetLogicalVolume();
			if (const auto* mesh = dynamic_cast<const Geant4::Mesh*>(daughter_volume->GetSolid()))
				placed.push_back({ mesh, daughter_rotation, daughter_translation });
			collect_placed_meshes(*daughter_volume, daughter_rotation, daughter_translation, placed);
		}
	}
}

void Geant4::add_geometry_channel(radfiled3d::CartesianRadiationField& field, const G4LogicalVolume& world_volume, int max_threads)
{
	std::vector<PlacedMesh> placed;
	collect_placed_meshes(world_volume, G4RotationMatrix(), G4ThreeVector(), placed);
	if (placed.empty())
		return;

	const auto start = std::chrono::steady_clock::now();
	auto channel = std::static_pointer_cast<radfiled3d::VoxelGridBuffer>(field.add_channel("geometry"));
	const glm::uvec3 counts = channel->get_voxel_counts();
	const double voxel_size = channel->get_voxel_dimensions().x;
	const Voxelization::VoxelGrid grid{
		counts,
		voxel_size,
		// world origin at the grid centre, like the detector's voxel mapping
		-0.5 * glm::dvec3(counts) * voxel_size
	};

	std::map<std::string, size_t> overlapping_per_type;
	for (const PlacedMesh& placed_mesh : placed) {
		const std::string& type = placed_mesh.mesh->getMesh()->getType();
		if (!channel->has_layer(type))
			channel->add_layer<uint8_t>(type, 0, "occupancy");
		const Geant4::PlacedMeshSurface surface(*placed_mesh.mesh, placed_mesh.rotation, placed_mesh.translation);
		overlapping_per_type[type] += Voxelization::mark_overlapping_voxels(surface, grid, channel->get_layer<uint8_t>(type), max_threads);
	}

	const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
	G4cout << "Voxelized " << placed.size() << " geometries into " << overlapping_per_type.size() << " type layer(s) in " << seconds << " s:" << G4endl;
	for (const auto& [type, voxels] : overlapping_per_type)
		G4cout << "  " << type << ": " << voxels << " overlapping voxel marks" << G4endl;
}
