#pragma once
#include "Geometry.hpp"
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/quaternion.hpp>
#include <optional>
#include <stdexcept>


namespace RadiationSimulation::Geometry::CArm {
	/** Thrown when the C-arm structure (meshes of Type CArm) collides with other geometry in the pose of a run: such a
	* pose cannot be taken by a real C-arm, and Geant4 cannot transport through overlapping volumes. The run is to be
	* rejected, not simulated. */
	class PoseCollision : public std::runtime_error {
	public:
		using std::runtime_error::runtime_error;
	};

	/// Exit code of RadField3D for a rejected pose (PoseCollision); the DatasetGenerator draws a new one.
	constexpr int POSE_COLLISION_EXIT_CODE = 3;

	/** The pose of the C-arm in one run, all in metres. */
	struct Pose {
		glm::quat rotation;       ///< RadiationSource::getCArmRotation()
		glm::vec3 focal_spot;     ///< RadiationSource::getLocation()
		glm::vec2 half_tangents;  ///< ISourceShape::getHalfTangents() of the beam
		glm::quat structure_rotation; ///< RadiationSource::getCArmStructureRotation(), used for meshes of Type CArm
	};

	/** Places a root mesh of type ImageDetector, XRayTube or CArm for `pose` (anything else is left alone). The mesh is
	* modelled in the base pose of the tube (below the table, beam along +Y): its own transform is turned by the
	* C-arm rotation about the isocentre (image detector, C-arm) or about the focal spot (tube, moved to it). An image
	* detector with an IsocenterDistance range is first moved along the beam so that its entrance face (lowest Y)
	* lies as far out as the beam still fits between its X and Z edges, within the range.
	*
	* This is the Geant4-free twin of Geant4::Mesh::fitToBeam and Geant4::Mesh::turnWithCArm (used by the analytical
	* simulator); CArmPlacementTest checks that both place a mesh identically.
	*
	* @return The distance (m) of the fitted entrance face from the isocentre, if the mesh was fitted.
	* @throws std::runtime_error if the beam axis does not cross a mesh that is to be fitted in the base pose.
	*/
	std::optional<float> place(Mesh& mesh, const Pose& pose);

	/** The eight corners (m) of the mesh's bounding box where its own transform puts them, in the order
	* (x, y, z) = lo/hi with x varying fastest. */
	std::vector<glm::vec3> placed_bounding_box_corners(const Mesh& mesh);
}
