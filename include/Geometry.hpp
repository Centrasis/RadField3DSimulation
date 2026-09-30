#pragma once
#include <vector>
#include <array>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/quaternion.hpp>
#include <string>
#include <memory>
#include <optional>
#include "Collisions.h"


namespace RadiationSimulation::Geant4 {
	class Mesh;
}

namespace RadiationSimulation::Geometry {
	enum class FaceType {
		Tri,
		Quad
	};

	class Face {
	public:
		virtual ~Face() = default;
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
		friend class RadiationSimulation::Geant4::Mesh;
	public:
		/// Type of a mesh without an explicit "Type" in its geometry description.
		static constexpr const char* DEFAULT_TYPE = "unknown";
		/// Type marking the patient. At most one mesh of a scene may have it.
		static constexpr const char* PATIENT_TYPE = "patient";
		/** Type of meshes that turn with the C-arm like an image detector: modelled for the base pose of the tube (below
		* the isocentre, beam along +Y), they are rotated about the isocentre with the tube, so they stay on the side
		* opposite the tube at their modelled distance, facing the beam. */
		static constexpr const char* IMAGE_DETECTOR_TYPE = "imagedetector";
		/** Type of meshes that move with the X-ray tube: modelled around the focal spot at the origin in the base pose
		* (beam along +Y), they are rotated like the tube and moved to its focal spot. */
		static constexpr const char* XRAY_TUBE_TYPE = "xraytube";
		/// Types name layers of the stored geometry channel, whose names hold at most this many characters.
		static constexpr size_t MAX_TYPE_LENGTH = 63;
	protected:
		std::vector<glm::vec3> vertices;
		const std::vector<Face*> faces;
		const std::string name;
		std::string material_name = "";
		glm::quat rotation = glm::quat_cast(glm::mat4(1.f));
		glm::vec3 position = glm::vec3(0.f);
		glm::vec3 scale    = glm::vec3(1.f);
		std::string type = Mesh::DEFAULT_TYPE;
		std::pair<glm::vec3, glm::vec3> bounding_box = { glm::vec3(0.f), glm::vec3(0.f) };
		std::vector<std::shared_ptr<Mesh>> children;
		std::optional<glm::vec2> isocenter_distance_range;

	public:
		Mesh(const std::vector<glm::vec3>& vertices, const std::vector<Face*>& faces, const std::string& name);

		// Mesh owns the heap-allocated Face* in `faces` (created by the geometry loader) and frees them here.
		~Mesh();

		// The faces are raw owning pointers: copying a Mesh would alias them and double-free on teardown.
		// Meshes are always handled via shared_ptr, so forbid copies rather than deep-copy the faces.
		Mesh(const Mesh&) = delete;
		Mesh& operator=(const Mesh&) = delete;

		void attachMaterialName(const std::string name);

		const std::string& getMaterialName() const;

		const std::vector<glm::vec3>& getVertices() const;

		const std::vector<Face*>& getFaces() const;

		/// The mesh's type, always in lower case.
		inline const std::string& getType() const { return this->type; };
		/** Sets the type, matched case-insensitively: it is stored in lower case, so e.g. "Shield" and "SHIELD" are one type.
		* @throws std::invalid_argument if the type is empty or longer than MAX_TYPE_LENGTH.
		*/
		void setType(const std::string& type);

		inline const bool isPatient() const { return this->type == Mesh::PATIENT_TYPE; };

		inline const bool isImageDetector() const { return this->type == Mesh::IMAGE_DETECTOR_TYPE; };
		inline const bool isXRayTube() const { return this->type == Mesh::XRAY_TUBE_TYPE; };

		/** Distance range (min, max) in metres of an image detector's entrance face from the isocentre. When set, the
		* detector is not kept at its modelled distance but moved as far out as the beam still fits on it (never
		* spilling past its edges), within this range. Unset: the modelled distance is kept. */
		inline const std::optional<glm::vec2>& getIsocenterDistanceRange() const { return this->isocenter_distance_range; }
		inline void setIsocenterDistanceRange(const glm::vec2& range) { this->isocenter_distance_range = range; }

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

		inline void setRotation(const glm::vec3& rotation) {
			this->rotation = glm::quat_cast(glm::mat4(1.f) * glm::rotate(glm::mat4(1.f), rotation.x, glm::vec3(1.f, 0.f, 0.f)) * glm::rotate(glm::mat4(1.f), rotation.y, glm::vec3(0.f, 1.f, 0.f)) * glm::rotate(glm::mat4(1.f), rotation.z, glm::vec3(0.f, 0.f, 1.f)));
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