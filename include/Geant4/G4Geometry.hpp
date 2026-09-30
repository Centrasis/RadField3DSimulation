#pragma once
#include <G4TessellatedSolid.hh>
#include "Geometry.hpp"
#include "Voxelization.hpp"
#include <radfiled3d/radiation_field.hpp>
#include <memory>
#include <G4PolyhedronArbitrary.hh>

class G4PVPlacement;
class G4LogicalVolume;
class G4Material;


namespace RadiationSimulation::Geant4 {
	class Mesh : public G4TessellatedSolid {
	protected:
		std::shared_ptr<Geometry::Mesh> mesh;
		std::vector<std::shared_ptr<Mesh>> children;
		const double length_unit;
		G4ThreeVector position = G4ThreeVector(0);
		G4RotationMatrix rotation;
		std::shared_ptr<G4LogicalVolume> volume;
		std::shared_ptr<G4PVPlacement> physical;
	public:
		Mesh(std::shared_ptr<Geometry::Mesh> mesh, double length_unit = 1.0);
		~Mesh() {
			G4cout << "Mesh destroyed" << G4endl;
		}
		std::shared_ptr<Geometry::Mesh> getMesh() const { return this->mesh; }
		void setMaterial(G4Material* material);
		void place(G4LogicalVolume* parent);
		std::shared_ptr<G4LogicalVolume> getVolume();
		const std::pair<glm::vec3, glm::vec3>& getBoundingBox() const;
		const G4RotationMatrix& getRotation() const { return this->rotation; }
		/** Turns the placement with the C-arm: the mesh as placed so far (its own transform) is rotated by `c_arm` about
		* the origin and then moved by `pivot` (G4 length units): the isocentre (zero) for an image detector, the focal
		* spot for a tube. */
		void turnWithCArm(const glm::quat& c_arm, const G4ThreeVector& pivot);
		/** Moves an image detector, placed in the base pose (beam along +Y, before turnWithCArm), along the beam axis:
		* its entrance face (lowest Y) goes as far from the isocentre as the beam still fits between its X and Z edges,
		* clamped to [min_distance, max_distance]. `half_tangents` are the beam's (see ISourceShape::getHalfTangents),
		* lengths in G4 units.
		* @return The distance of the entrance face from the isocentre.
		* @throws std::runtime_error if the beam axis misses the detector. */
		double fitToBeam(double source_distance, const glm::vec2& half_tangents, double min_distance, double max_distance);
		/** The eight corners of the mesh's bounding box where its placement puts them (G4 length units). */
		std::vector<G4ThreeVector> placedBoundingBoxCorners() const;
		double getLengthUnit() const { return this->length_unit; }
		inline const std::vector<std::shared_ptr<Mesh>>& getChildren() const { return this->children; }
		/** Materials of this volume and of all its children, each once, as set by the SceneConstructor. */
		std::vector<const G4Material*> getMaterials() const;
	};

	/** Surface of a placed Geant4::Mesh volume: its tessellated facets moved into the world by all placements from the
	* volume up to the world volume, in metres.
	*/
	class PlacedMeshSurface : public Geometry::Voxelization::ISurface {
	protected:
		const G4TessellatedSolid& solid;
		const G4RotationMatrix rotation;
		const G4ThreeVector translation;
	public:
		PlacedMeshSurface(const G4TessellatedSolid& solid, const G4RotationMatrix& rotation, const G4ThreeVector& translation);
		virtual size_t polygon_count() const override;
		virtual size_t vertex_count(size_t polygon) const override;
		virtual glm::dvec3 vertex(size_t polygon, size_t vertex) const override;
	};

	/** Voxelizes the geometry placed below `world_volume` into the channel "geometry" of `field`, using the field's voxel
	* grid. Every mesh Type gets its own 8-bit layer named like the type, holding 255 where a voxel overlaps any mesh of
	* that type and 0 elsewhere (see Geometry::Voxelization::mark_overlapping_voxels). Meshes are processed one at a
	* time, so the additional memory is the layers plus the distance field of the largest mesh.
	* @param field The field to add the channel to.
	* @param world_volume The logical world volume of the constructed Geant4 scene.
	* @param max_threads Maximum number of worker threads. -1 uses all available cores.
	*/
	void add_geometry_channel(radfiled3d::CartesianRadiationField& field, const G4LogicalVolume& world_volume, int max_threads = -1);
}
