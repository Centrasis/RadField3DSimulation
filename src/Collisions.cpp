#include "Collisions.h"
#define GLM_ENABLE_EXPERIMENTAL
#include "glm/glm.hpp"
#include "glm/gtx/intersect.hpp"
#include <array>
#include <algorithm>


glm::vec3 Collisions::Line::closest_point(const glm::vec3& p) const
{
	glm::vec3 diff = this->p2 - this->p1;
	float t = glm::dot(p - this->p1, diff) / glm::dot(diff, diff);
	return this->p1 + std::min(std::max(t, 0.f), 1.f) * diff;
}

float Collisions::Line::length() const
{
	return glm::length(this->p1 - this->p2);
}

glm::vec3 Collisions::Line::direction() const
{
	return glm::normalize(this->p2 - this->p1);
}

Collisions::Box::Box(const glm::vec3& min, const glm::vec3& max) 
	: min(min),
	  max(max),
	  planes()
	  /*planes(
		  {
			  // 1st Quad
			  Plane(
				min, glm::vec3(min.x, max.y, min.z), glm::vec3(max.x, max.y, min.z), glm::vec3(max.x, min.y, max.z)
			  ),

			  // 2nd Quad
			  Plane(
				min, glm::vec3(max.x, min.y, min.z), glm::vec3(max.x, min.y, max.z), glm::vec3(min.x, min.y, max.z)
			  ),

			  // 3rd Quad
			  Plane(
				glm::vec3(min.x, min.y, max.z), glm::vec3(max.x, min.y, max.z), max, glm::vec3(min.x, max.y, max.z)
			  ),

			  // 4th Quad
			  Plane(
				min, glm::vec3(min.x, max.y, min.z), glm::vec3(min.x, max.y, max.z), glm::vec3(min.x, min.y, max.z)
			  ),

			  // 5th Quad
			  Plane(
				glm::vec3(max.x, min.y, min.z), glm::vec3(max.x, max.y, min.z), max, glm::vec3(max.x, min.y, max.z)
			  ),

			  // 6th Quad
			  Plane(
				  glm::vec3(min.x, max.y, min.z), glm::vec3(max.x, max.y, min.z), max, glm::vec3(min.x, max.y, max.z)
			  )
		  }
	)*/
{
}

bool Collisions::Box::intersect(const Line& line) const
{
	for (const auto& p : this->planes) {
		if (p.intersect(line))
			return true;
	}
	return false;
}

Collisions::Triangle::Triangle(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3)
	: p1(p1),
	  p2(p2),
	  p3(p3),
	  normal(glm::normalize(glm::cross(p2 - p1, p3 - p1)))
{
}

bool Collisions::Triangle::intersect(const Line& line) const
{
	const glm::vec3 p1p2 = line.p2 - line.p1;
	const float line_length = glm::length(p1p2);
	glm::vec3 plane_intersection;
	
	if (line_length <= this->POINT_THRESHOLD_LENGTH) {
		// if the line segment is essentially a point and has no real extend
		// Treat line like a point and check that
		plane_intersection = line.p1;
	} else {
		// if line segment has an extend treat it as a line
		const glm::vec3 plane_normal = glm::cross(p2 - p1, p3 - p1);
		const glm::vec3 ray_dir = glm::normalize(p1p2);
		if (glm::dot(plane_normal, ray_dir) == 0.f) // check if ray is parallel to the triangle plane
			return false;

		const float alpha = (glm::dot(plane_normal, p1) - glm::dot(plane_normal, line.p1)) / glm::dot(plane_normal, ray_dir);
		plane_intersection = line.p1 + alpha * ray_dir;

		if (line_length != glm::length(plane_intersection - line.p1) + glm::length(plane_intersection - line.p2)) // check if plane intersection is part of line segment
			return false;
	}

	const glm::vec3 np1 = this->p1 - plane_intersection;
	const glm::vec3 np2 = this->p2 - plane_intersection;
	const glm::vec3 np3 = this->p3 - plane_intersection;

	const glm::vec3 u = glm::cross(np2, np3);
	const glm::vec3 v = glm::cross(np3, np1);
	const glm::vec3 w = glm::cross(np1, np2);

	return !(glm::dot(u, v) < 0.f || glm::dot(u, w) < 0.f);
}

Collisions::Plane::Plane(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3)
	: base_point(p1),
	  d1(glm::normalize(p2 - p1)),
	  d2(glm::normalize(p3 - p1)),
	  normal(glm::normalize(glm::cross(d1, d2)))	  
{
}

bool Collisions::Plane::intersect(const Line& line) const
{
	float distance = 0.f;
	if (glm::intersectRayPlane(line.p1, line.direction(), this->base_point, this->normal, distance) && distance >= 0.f) {
		const glm::vec3 intersection = line.p1 + distance * line.direction();
		return line.length() == (glm::length(intersection - line.p1) + glm::length(intersection - line.p2));
	}
	return false;
}

bool Collisions::Plane::intersect_directed_ray(const Line& line, glm::vec3& intersection) const
{
	float distance = 0.f;
	if (glm::intersectRayPlane(line.p1, line.direction(), this->base_point, this->normal, distance) && distance >= 0.f) {
		intersection = line.p1 + distance * line.direction();
		return true;
	}
	return false;
}

Collisions::Capsule::Capsule(const glm::vec3& p1, const glm::vec3& p2, float radius)
	: base_line(p1, p2),
	  radius(radius)
{
}

bool Collisions::Capsule::intersect(const glm::vec3& p) const
{
	const glm::vec3 cp = this->base_line.closest_point(p);
	const float dist = glm::length(cp - p);

	return dist <= radius;
}
