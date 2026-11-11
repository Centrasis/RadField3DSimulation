#include "CArm.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <glm/gtc/matrix_transform.hpp>

using namespace RadiationSimulation::Geometry;


std::vector<glm::vec3> CArm::placed_bounding_box_corners(const Mesh& mesh)
{
	glm::vec3 lo(std::numeric_limits<float>::infinity()), hi(-std::numeric_limits<float>::infinity());
	for (const glm::vec3& v : mesh.getVertices()) {
		const glm::vec3 scaled = v * mesh.getScale();
		lo = glm::min(lo, scaled);
		hi = glm::max(hi, scaled);
	}
	const glm::mat3 rotation = glm::mat3_cast(mesh.getRotation());
	std::vector<glm::vec3> corners;
	for (int i = 0; i < 8; i++) {
		const glm::vec3 local((i & 1) ? hi.x : lo.x, (i & 2) ? hi.y : lo.y, (i & 4) ? hi.z : lo.z);
		corners.push_back(rotation * local + mesh.getPosition());
	}
	return corners;
}

namespace {
	float fit_to_beam(Mesh& mesh, const CArm::Pose& pose, float min_distance, float max_distance)
	{
		float entrance = std::numeric_limits<float>::infinity();
		float x_lo = entrance, x_hi = -entrance, z_lo = entrance, z_hi = -entrance;
		for (const glm::vec3& c : CArm::placed_bounding_box_corners(mesh)) {
			entrance = std::min(entrance, c.y);
			x_lo = std::min(x_lo, c.x);
			x_hi = std::max(x_hi, c.x);
			z_lo = std::min(z_lo, c.z);
			z_hi = std::max(z_hi, c.z);
		}
		// the beam axis is the Y axis; the edge nearest to it limits the beam
		const float reach_x = std::min(x_hi, -x_lo);
		const float reach_z = std::min(z_hi, -z_lo);
		if (!(reach_x > 0.f && reach_z > 0.f))
			throw std::runtime_error("Image detector \"" + mesh.getName() + "\": the beam axis does not cross it in the base pose, so it cannot be fitted to the beam.");
		// the base pose maps the source frame's x onto X and its y onto Z, so the collimator's edges run along the detector's
		const float source_distance = glm::length(pose.focal_spot);
		const float fit = std::min(reach_x / pose.half_tangents.x, reach_z / pose.half_tangents.y) - source_distance;
		const float distance = std::clamp(fit, min_distance, max_distance);
		mesh.setPosition(mesh.getPosition() + glm::vec3(0.f, distance - entrance, 0.f));
		return distance;
	}
}

std::optional<float> CArm::place(Mesh& mesh, const Pose& pose)
{
	if (!mesh.turnsWithCArm())
		return std::nullopt;
	std::optional<float> distance;
	if (mesh.isImageDetector() && mesh.getIsocenterDistanceRange().has_value()) {
		const glm::vec2 range = mesh.getIsocenterDistanceRange().value();
		distance = fit_to_beam(mesh, pose, range.x, range.y);
	}
	const glm::quat turn = glm::normalize(mesh.isCArm() ? pose.structure_rotation : pose.rotation);
	mesh.setRotation(turn * mesh.getRotation());
	mesh.setPosition((mesh.isXRayTube() ? pose.focal_spot : glm::vec3(0.f)) + turn * mesh.getPosition());
	return distance;
}
