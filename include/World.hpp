#pragma once
#include <memory>
#include <vector>
#include <string>
#include <glm/glm.hpp>

namespace RadiationSimulation {
	class RadiationSource;
	class Mesh;
	class RadiationFieldDetector;

	/**
	 * @brief Structure to hold information about the world.
	 */
	struct WorldInfo {
		const std::string material; ///< Material of the world
		const glm::vec3 dimensions; ///< Dimensions of the world

		/**
		 * @brief Constructor for WorldInfo.
		 * @param material Material of the world.
		 * @param dimensions Dimensions of the world.
		 */
		WorldInfo(const std::string& material, const glm::vec3& dimensions)
			: material(material), dimensions(dimensions)
		{
		}
	};

	/**
	 * @brief Class representing the simulation world.
	 */
	class World {
		friend class RadiationSimulator;
		friend class G4RadiationSimulationHandler;
		friend class G4World;
	protected:
		std::shared_ptr<RadiationSource> radiation_source; ///< Radiation source in the world
		std::shared_ptr<RadiationFieldDetector> field_detector; ///< Detector for radiation field scoring
		std::vector<std::shared_ptr<Mesh>> geometries; ///< Geometries in the world
		std::shared_ptr<Mesh> patient = std::shared_ptr<Mesh>(NULL); ///< Patient mesh
		static std::shared_ptr<World> instance; ///< Singleton instance of the world
		static std::unique_ptr<WorldInfo> world_info; ///< Information about the world
	public:
		/**
		 * @brief Constructor for World.
		 */
		World();

		/**
		 * @brief Initialize the world.
		 */
		static void initialize();

		/**
		 * @brief Set world information.
		 * @param info Unique pointer to WorldInfo.
		 */
		inline static void set_world_info(std::unique_ptr<WorldInfo>& info) { World::world_info = std::move(info); }

		/**
		 * @brief Get world information.
		 * @return Unique pointer to WorldInfo.
		 */
		inline static const std::unique_ptr<WorldInfo>& get_world_info() { return World::world_info; }

		/**
		 * @brief Get the singleton instance of the world.
		 * @return Shared pointer to the World instance.
		 */
		inline static std::shared_ptr<World> Get() { return World::instance; }

		/**
		 * @brief Set the radiation field detector.
		 * @param field_detector Shared pointer to RadiationFieldDetector.
		 */
		virtual void set_radiation_field_detector(std::shared_ptr<RadiationFieldDetector> field_detector) { this->field_detector = field_detector; }

		/**
		 * @brief Get the radiation field detector.
		 * @return Shared pointer to RadiationFieldDetector.
		 */
		virtual std::shared_ptr<RadiationFieldDetector> get_radiation_field_detector() const { return this->field_detector; }

		/**
		 * @brief Set the radiation source.
		 * @param source Shared pointer to RadiationSource.
		 */
		virtual void set_radiation_source(std::shared_ptr<RadiationSource> source) { this->radiation_source = source; }

		/**
		 * @brief Set the geometries in the world.
		 * @param geoms Vector of shared pointers to Mesh.
		 */
		virtual void set_geometries(const std::vector<std::shared_ptr<Mesh>>& geoms);

		/**
		 * @brief Get the radiation source.
		 * @return Shared pointer to RadiationSource.
		 */
		virtual std::shared_ptr<RadiationSource> get_radiation_source() const { return this->radiation_source; }

		/**
		 * @brief Get the geometries.
		 * @return Vector of shared pointers to Mesh.
		 */
		virtual const std::vector<std::shared_ptr<Mesh>>& get_geometries() const { return this->geometries; }

		/**
		 * @brief Get the patient mesh.
		 * @return Shared pointer to Mesh.
		 */
		virtual std::shared_ptr<Mesh> get_patient() const { return this->patient; }
	};
}
