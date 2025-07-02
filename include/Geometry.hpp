#pragma once
#include <vector>
#include <array>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/quaternion.hpp>
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
		glm::quat rotation = glm::quat_cast(glm::mat4(1.f));
		glm::vec3 position = glm::vec3(0.f);
		glm::vec3 scale    = glm::vec3(1.f);
		bool bis_patient = false;
		struct {
			bool bis_source = false;
			float concentric_distance = 0.f;
			glm::vec2 rotation_offset_radians = glm::vec2(0.f);
		} source_info;
		std::pair<glm::vec3, glm::vec3> bounding_box = { glm::vec3(0.f), glm::vec3(0.f) };
		std::vector<std::shared_ptr<Mesh>> children;

	public:
		Mesh(const std::vector<glm::vec3>& vertices, const std::vector<Face*>& faces, const std::string& name);

		void attachMaterialName(const std::string name);

		const std::string& getMaterialName() const;

		const std::vector<glm::vec3>& getVertices() const;

		const std::vector<Face*>& getFaces() const;

		inline const bool isPatient() const { return this->bis_patient; };

		inline const bool isSource() const { return this->source_info.bis_source; };
		inline const float getSourceConcentricDistance() const { return this->source_info.concentric_distance; };
		const glm::vec2& getSourceRotationOffset() const { return this->source_info.rotation_offset_radians; };
		inline void markAsSource(float concentric_distance, const glm::vec2& rotation_offset_radians) {
			this->source_info.bis_source = true;
			this->source_info.concentric_distance = concentric_distance;
			this->source_info.rotation_offset_radians = rotation_offset_radians;
		};

		inline void markAsPatient() { this->bis_patient = true; };

		const std::string& getName() const;

		size_t vertexCount() const;

		size_t faceCount() const;

		const std::pair<glm::vec3, glm::vec3>& getAxisAlignedBoundingBox() const { return this->bounding_box; }

		std::shared_ptr<OrientedBoundingBox> getOrientedBoundingBox() const;

		std::shared_ptr<Collisions::Capsule> getBoundingCapsule() const;

		inline const glm::quat& getRotation() const {
			return this->rotation;
		};

		inline void setRotation(const glm::quat& rotation) {
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