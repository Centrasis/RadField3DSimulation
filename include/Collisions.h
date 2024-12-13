#pragma once
#include <glm/vec3.hpp>
#include <vector>
#include <array>


namespace Collisions {
	struct Line {
		const glm::vec3 p1;
		const glm::vec3 p2;

		Line(const glm::vec3& p1, const glm::vec3& p2)
			: p1(p1), p2(p2)
		{}

		glm::vec3 closest_point(const glm::vec3& p) const;

		float length() const;

		glm::vec3 direction() const;
	};

	struct Triangle {
		// Lines with a extend shorter than this threshold will be treated like points
		const float POINT_THRESHOLD_LENGTH = 0.000001f;
		const glm::vec3 p1;
		const glm::vec3 p2;
		const glm::vec3 p3;
		const glm::vec3 normal;

		Triangle(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);

		bool intersect(const Line& line) const;
	};

	struct Plane {
		const glm::vec3 base_point;
		const glm::vec3 d1;
		const glm::vec3 d2;
		const glm::vec3 normal;
		
		Plane(const glm::vec3& p1, const glm::vec3& p2, const glm::vec3& p3);
		bool intersect(const Line& line) const;
		bool intersect_directed_ray(const Line& line, glm::vec3& intersection) const;
	};

	struct Box {
		const glm::vec3 min;
		const glm::vec3 max;
		const std::vector<Plane> planes;

		Box(const glm::vec3& min, const glm::vec3& max);

		bool intersect(const Line& line) const;
	};

	struct Capsule {
		const float radius;
		const Line base_line;

		Capsule(const glm::vec3& p1, const glm::vec3& p2, float radius);
		bool intersect(const glm::vec3& p) const;
	};
}