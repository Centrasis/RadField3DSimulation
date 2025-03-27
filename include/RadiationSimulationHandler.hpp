#pragma once
#include <memory>
#include <future>
#include "RadFiled3D/RadiationField.hpp"
#include "RadiationSource.hpp"
#include <CLHEP/Random/MTwistEngine.h>
#include <RadFiled3D/GridTracer.hpp>
#include <G4RunManager.hh>

#ifdef WITH_GEANT4_UIVIS
class G4UImanager;
#include <G4VisExecutive.hh>
#endif
class G4VModularPhysicsList;

namespace RadiationSimulation {
	class Mesh;
	class G4RadiationFieldDetector;

	/**
	 * @brief Handler for a radiation simulation process.
	 */
	class RadiationSimulationHandler {
	protected:
		struct {
			glm::vec3 radiation_field_dimensions; ///< Dimensions of the radiation field.
			glm::vec3 radiation_field_voxel_dimensions; ///< Voxel dimensions of the radiation field.
			float radiation_field_max_energy; ///< Maximum energy of the radiation field.
			float energy_resolution; ///< Energy resolution of the radiation field.
		} radiation_field_resolution;
	public:
		/**
		 * @brief Initialize the simulation handler.
		 * @return True if initialization is successful, false otherwise.
		 */
		virtual bool initialize() = 0;

		/**
		 * @brief Finalize the simulation handler.
		 */
		virtual void finalize() = 0;

		/**
		 * @brief Deinitialize the simulation handler.
		 */
		virtual void deinitialize() = 0;

		/**
		 * @brief Set the resolution of the radiation field.
		 * @param radiation_field_dimensions Dimensions of the radiation field.
		 * @param radiation_field_voxel_dimensions Voxel dimensions of the radiation field.
		 * @param radiation_field_max_energy Maximum energy of the radiation field. Defines the histogram bins of the spectra layer.
		 * @param energy_resolution Energy resolution of the radiation field. Defines the histogram bins of the spectra layer.
		 */
		void set_radiation_field_resolution(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, float radiation_field_max_energy, float energy_resolution);

		/**
		 * @brief Add geometry to the simulation.
		 * @param meshes Vector of shared pointers to Mesh objects.
		 */
		virtual void add_geometry(const std::vector<std::shared_ptr<Mesh>>& meshes) = 0;

		/**
		 * @brief Simulate the radiation field.
		 * @param n_particles Maximum number of particles to simulate.
		 * @param tracing_algorithm Algorithm to use for grid tracing.
		 * @return Shared pointer to the simulated radiation field.
		 */
		virtual std::shared_ptr<RadFiled3D::IRadiationField> simulate_radiation_field(size_t n_particles = 1e+6, RadFiled3D::GridTracerAlgorithm tracing_algorithm = RadFiled3D::GridTracerAlgorithm::SAMPLING) = 0;

		/**
		 * @brief Display the graphical user interface. Works only if the CMake option WITH_GEANT4_UIVIS is enabled.
		 */
		virtual void display_gui() {};

		/**
		 * @brief Update the graphical user interface. Works only if the CMake option WITH_GEANT4_UIVIS is enabled.
		 */
		virtual void update_gui() {};

		/**
		 * @brief Destructor for the simulation handler.
		 */
		virtual ~RadiationSimulationHandler() {};

		/**
		 * @brief Add a callback to be executed every n particles.
		 * @param callback Function to be called.
		 * @param n_particles Number of particles after which the callback is executed.
		 */
		virtual void add_callback_every_n_particles(std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)> callback, size_t n_particles) {
			throw std::runtime_error("Not implemented");
		};
	};

	/**
	 * @brief Handler for Geant4-based radiation simulation.
	 */
	class G4RadiationSimulationHandler : public RadiationSimulationHandler {
	protected:
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
		 * @brief Initialize the Geant4 simulation handler.
		 * @return True if initialization is successful, false otherwise.
		 */
		virtual bool initialize() override;

		/**
		 * @brief Finalize the Geant4 simulation handler.
		 */
		virtual void finalize() override;

		/**
		 * @brief Display the graphical user interface.
		 */
		virtual void display_gui() override;

		/**
		 * @brief Update the graphical user interface.
		 */
		virtual void update_gui() override;

		/**
		 * @brief Add geometry to the simulation.
		 * @param meshes Vector of shared pointers to Mesh objects.
		 */
		virtual void add_geometry(const std::vector<std::shared_ptr<Mesh>>& meshes) override;

		/**
		 * @brief Simulate the radiation field.
		 * @param n_particles Number of particles to simulate.
		 * @param tracing_algorithm Algorithm to use for grid tracing.
		 * @return Shared pointer to the simulated radiation field.
		 */
		virtual std::shared_ptr<RadFiled3D::IRadiationField> simulate_radiation_field(size_t n_particles = 1e+6, RadFiled3D::GridTracerAlgorithm tracing_algorithm = RadFiled3D::GridTracerAlgorithm::SAMPLING) override;

		/**
		 * @brief Deinitialize the Geant4 simulation handler.
		 */
		virtual void deinitialize() override;

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
		virtual void add_callback_every_n_particles(std::function<void(std::shared_ptr<RadFiled3D::IRadiationField>, size_t)> callback, size_t n_particles) override;
	};
}