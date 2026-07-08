#pragma once
#include <memory>
#include "RadFiled3D/RadiationField.hpp"
#include "RadiationSource.hpp"
#include <CLHEP/Random/MTwistEngine.h>
#include <RadFiled3D/GridTracer.hpp>
#include <G4RunManager.hh>

#ifdef WITH_GEANT4_UIVIS
#include <G4UImanager.hh>
#include <G4VisExecutive.hh>
#endif
class G4VModularPhysicsList;

namespace RadiationSimulation {
	class Mesh;
	class G4RadiationFieldDetector;

	/**
	 * @brief Handler for Geant4-based radiation simulation.
	 */
	class G4RadiationSimulationHandler {
	protected:
		struct {
			glm::vec3 radiation_field_dimensions; ///< Dimensions of the radiation field.
			glm::vec3 radiation_field_voxel_dimensions; ///< Voxel dimensions of the radiation field.
			float radiation_field_max_energy; ///< Maximum energy of the radiation field.
			float energy_resolution; ///< Energy resolution of the radiation field.
			struct {
				float threshold = 0.1f; ///< Statistical error threshold for the simulation.
				float enforcement_ratio = 0.9f; ///< Statistical error enforcement over the top particles inside the ratio.
				float enforcement_resolution = 1.f; ///< Statistical error resolution for enforcement e.g. check for every n-th voxel with 1.f -> every voxel and 0.5f -> every second, ....
			} statistical_error;
			glm::uvec2 angular_resolution = glm::uvec2(0); //< Number of segments (phi, theta) for angular distribution per voxel. 0 = disabled
		} radiation_field_resolution;

		std::vector<std::shared_ptr<Mesh>> meshes;
		int cpu_count; ///< Number of CPU cores to use.
		bool run_mgr_initialized = false; ///< Flag indicating if the run manager is initialized.
		G4VModularPhysicsList* physics; ///< Pointer to the physics list.
		const int base_particle_count = 1e+6; ///< Base particle count for the simulation.
		bool has_ui = false; ///< Flag indicating if the UI is available.
		CLHEP::MTwistEngine random_generator; ///< Random number generator.
		std::unique_ptr<G4RunManager> G4mgr; ///< Unique pointer to the Geant4 run manager.
#ifdef WITH_GEANT4_UIVIS
		std::shared_ptr<G4UImanager> G4UIManager; ///< Shared pointer to the Geant4 UI manager.
		std::unique_ptr<G4VisExecutive> G4VisManager; ///< Unique pointer to the Geant4 visualization manager.
#endif
		std::shared_ptr<G4RadiationFieldDetector> field_detector; ///< Shared pointer to the radiation field detector.
		std::vector<std::pair<size_t, std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)>>> callbacks; ///< Vector of callbacks.

	public:
		G4RadiationSimulationHandler(const int cpu_count = -1);

		/**
		 * @brief Set the resolution of the radiation field.
		 */
		void set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution, float statistical_error_threshold, float statistical_error_enforcement_ratio, glm::uvec2 angular_resolution = glm::uvec2(0));

		/**
		 * @brief Initialize the Geant4 simulation handler.
		 * @return True if initialization is successful, false otherwise.
		 */
		virtual bool initialize();

		/**
		 * @brief Finalize the Geant4 simulation handler.
		 */
		virtual void finalize();

		/**
		 * @brief Display the graphical user interface.
		 */
		virtual void display_gui();

		/**
		 * @brief Update the graphical user interface.
		 */
		virtual void update_gui();

		/**
		 * @brief Add geometry to the simulation.
		 * @param meshes Vector of shared pointers to Mesh objects.
		 */
		virtual void add_geometry(const std::vector<std::shared_ptr<Mesh>>& meshes);

		/**
		 * @brief Simulate the radiation field.
		 * @param n_particles Number of particles to simulate.
		 * @param tracing_algorithm Algorithm to use for grid tracing.
		 * @return Shared pointer to the simulated radiation field.
		 */
		virtual std::shared_ptr<RadFiled3D::IRadiationField> simulate_radiation_field(size_t n_particles = 1e+6, RadFiled3D::GridTracerAlgorithm tracing_algorithm = RadFiled3D::GridTracerAlgorithm::SAMPLING);

		/**
		 * @brief Deinitialize the Geant4 simulation handler.
		 */
		virtual void deinitialize();

		/**
		 * @brief Get the used field detector.
		 * @return Shared pointer to the used field detector.
		 */
		inline std::shared_ptr<G4RadiationFieldDetector> get_used_field_detector() const { return this->field_detector; }

		/**
		 * @brief Add a callback to be executed every n particles.
		 * @param callback Function to be called.
		 * @param n_particles Number of particles after which the callback is executed.
		 */
		virtual void add_callback_every_n_particles(std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)> callback, size_t n_particles);
	};
}