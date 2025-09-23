#pragma once
#include <memory>
#include "RadiationSimulationHandler.hpp"
#include "Geometry.hpp"
#include <RadFiled3D/storage/RadiationFieldStore.hpp>
#include <RadFiled3D/GridTracer.hpp>

namespace RadiationSimulation {
	class World;
	struct WorldInfo;

	/**
	 * @brief Enum class for different types of radiation handlers.
	 */
	enum class RadiationHandlerType {
		Geant4MedicalXRay ///< Handler for Geant4 Medical X-Ray simulations.
	};

	/**
	 * @brief Class for managing radiation simulations.
	 */
	class RadiationSimulator {
	protected:
		/**
		 * @brief Shared pointer to the radiation simulation handler.
		 */
		static std::shared_ptr<RadiationSimulationHandler> handler;

		/**
		 * @brief Flag indicating if the simulator is busy.
		 */
		static bool bIsBusy;

	public:
		/**
		 * @brief Initializes the radiation simulation handler.
		 * @param handler_type The type of radiation handler to initialize.
		 * @return Shared pointer to the initialized radiation simulation handler.
		 */
		static std::shared_ptr<RadiationSimulationHandler> initialize(const RadiationHandlerType handler_type, const int cpu_count = -1);

		/**
		 * @brief Simulates the radiation field. Shall support multi-threading, but is not required to allow a non-blocking call.
		 * @param n_particles Number of particles to simulate.
		 * @param tracing_algorithm Algorithm to use for grid tracing.
		 * @return Promise of a shared pointer to the simulated radiation field.
		 */
		static std::promise<std::shared_ptr<RadFiled3D::IRadiationField>> simulate_radiation_field(size_t n_particles = 1e+6, RadFiled3D::GridTracerAlgorithm tracing_algorithm = RadFiled3D::GridTracerAlgorithm::SAMPLING);

		/**
		 * @brief Adds a callback function to be called every n particles.
		 * @param callback The callback function.
		 * @param n_particles Number of particles after which the callback is called.
		 */
		static void add_callback_every_n_particles(std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)> callback, size_t n_particles);

		/**
		 * @brief Displays the graphical user interface.
		 */
		static void display_gui();

		/**
		 * @brief Adds geometry to the simulation.
		 * @param meshes Vector of shared pointers to the meshes to add.
		 */
		static void add_geometry(const std::vector<std::shared_ptr<Mesh>>& meshes);

		/**
		 * @brief Adds a single mesh to the simulation.
		 * @param mesh Shared pointer to the mesh to add.
		 */
		static void add_geometry(std::shared_ptr<Mesh> mesh);

		/**
		 * @brief Adds a radiation source to the simulation.
		 * @param source Shared pointer to the radiation source to add.
		 */
		static void add_radiation_source(std::shared_ptr<RadiationSource> source);

		/**
		 * @brief Sets the resolution of the radiation field.
		 * @param radiation_field_dimensions Dimensions of the radiation field.
		 * @param radiation_field_voxel_dimensions Dimensions of the voxels in the radiation field.
		 * @param radiation_field_max_energy Maximum energy of the radiation field.
		 * @param energy_resolution Resolution of the energy in the radiation field.
		 * @param statistical_error_threshold Statistical error threshold that needs to be fullfilled by a certain amount of voxels for early stopping of the simulation.
		 * @param statistical_error_enforcement_ratio Ratio of voxels that need to fullfill the statistical error threshold. The enforcement is done on the top x% of voxels sorted by their errors.
		 */
		static void set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution, float statistical_error_threshold, float statistical_error_enforcement_ratio);

		/**
		 * @brief Deinitializes the radiation simulator.
		 */
		static void deinitialize();

		/**
		 * @brief Checks if the simulator is busy.
		 * @return True if the simulator is busy, false otherwise.
		 */
		inline static bool is_busy() { return RadiationSimulator::bIsBusy; };

		/**
		 * @brief Sets the world information for the simulation.
		 * @param info Unique pointer to the world information.
		 */
		static void set_world_info(std::unique_ptr<WorldInfo> info);
	};
}