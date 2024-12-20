#pragma once
#include <vector>
#include <array>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <string>
#include <memory>
#include "Collisions.h"


namespace RadiationSimulation {
	enum class FaceType {
		Tri,
		Quad
	};

	class Face {
	public:
		virtual FaceType getType() const = 0;
	};

	struct OrientedBoundingBox {
		const glm::vec3 centroid;
		const std::array<glm::vec3, 3> axes;
		const glm::vec3 min_coeffs;
		const glm::vec3 max_coeffs;

		OrientedBoundingBox(const glm::vec3& centroid, const std::array<glm::vec3, 3>& axes, const glm::vec3& min_coeffs, const glm::vec3& max_coeffs);
	};

	template<class T, FaceType type>
	class FaceT : public Face {
	protected:
		const T indices;
	public:
		FaceT(const T& indices) : indices(indices) {};

		const T& getIndices() const { return this->indices; }
		virtual FaceType getType() const override { return type; }
	};

	typedef FaceT<glm::uvec3, FaceType::Tri> TriFace;
	typedef FaceT<glm::uvec4, FaceType::Quad> QuadFace;

	class Mesh {
		friend class G4Mesh;
	protected:
		std::vector<glm::vec3> vertices;
		const std::vector<Face*> faces;
		const std::string name;
		std::string material_name = "";
		glm::vec3 rotation = glm::vec3(0.f);
		glm::vec3 position = glm::vec3(0.f);
		glm::vec3 scale    = glm::vec3(1.f);
		bool bis_patient = false;
		std::pair<glm::vec3, glm::vec3> bounding_box = { glm::vec3(0.f), glm::vec3(0.f) };
		std::vector<std::shared_ptr<Mesh>> children;

	public:
		Mesh(const std::vector<glm::vec3>& vertices, const std::vector<Face*>& faces, const std::string& name);

		void attachMaterialName(const std::string name);

		const std::string& getMaterialName() const;

		const std::vector<glm::vec3>& getVertices() const;

		const std::vector<Face*>& getFaces() const;

		inline const bool is_patient() const { return this->bis_patient; };

		void mark_patient() { this->bis_patient = true; };

		const std::string& getName() const;

		size_t vertexCount() const;

		size_t faceCount() const;

		const std::pair<glm::vec3, glm::vec3>& getAxisAlignedBoundingBox() const { return this->bounding_box; }

		std::shared_ptr<OrientedBoundingBox> getOrientedBoundingBox() const;

		std::shared_ptr<Collisions::Capsule> getBoundingCapsule() const;

		inline const glm::vec3& getRotation() const {
			return this->rotation;
		};

		inline void setRotation(const glm::vec3& rotation) {
			this->rotation = rotation;
		}

		inline const glm::vec3& getPosition() const {
			return this->position;
		};

		inline void setPosition(const glm::vec3& pos) {
			this->position = pos;
		}

		inline const glm::vec3& getScale() const {
			return this->scale;
		};

		inline void setScale(const glm::vec3& scale) {
			this->scale = scale;
		}

		inline const std::vector<std::shared_ptr<Mesh>>& getChildren() const {
			return this->children;
		}

		inline void addChild(std::shared_ptr<Mesh> child) {
			this->children.push_back(child);
		}
	};
}