#pragma once
#include <memory>
#include <vector>
#include <radfiled3d/radiation_field.hpp>
#include <G4VSensitiveDetector.hh>
#include <glm/glm.hpp>
#include <algorithm>
#include <G4UserSteppingAction.hh>
#include <G4VUserActionInitialization.hh>
#include "Collisions.h"
#include <unordered_set>
#include <map>
#include <G4SystemOfUnits.hh>
#include <utils/RollingBuffer.hpp>
#include <radfiled3d/grid_tracer.hpp>
#include <G4EmCalculator.hh>
#include <shared_mutex>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include "Statistics.hpp"
#include "VMFTrainer.hpp"
#include <algorithm>
#include <numbers>


class G4RunManager;
#ifdef WITH_GEANT4_UIVIS
class G4UImanager;
class G4VisExecutive;
#endif
class G4LogicalVolume;

namespace RadiationSimulation {
	class RadiationSource;
}

namespace RadiationSimulation::Geant4 {
	class RadiationSource;

	enum class TrackStage : char {
		BEAM,
		SCATTER,
		PATIENT
	};

	/** Lets any number of threads score concurrently and one thread pause all of them, e.g. to take a consistent
	* snapshot of the field. Pausing waits for the steps in progress, then keeps new steps out until resumed; it takes
	* precedence over new steps, so a snapshot cannot be starved by continuously scoring workers.
	*/
	class ScoringGate {
		std::atomic<size_t> active{ 0 };
		std::atomic<bool> paused{ false };
		std::mutex mutex;
		std::condition_variable changed;
	public:
		void enter() {
			while (true) {
				this->active.fetch_add(1);
				if (!this->paused.load())
					return;
				this->leave();
				std::unique_lock lock(this->mutex);
				this->changed.wait(lock, [this] { return !this->paused.load(); });
			}
		}
		void leave() {
			if (this->active.fetch_sub(1) == 1 && this->paused.load()) {
				std::lock_guard lock(this->mutex);
				this->changed.notify_all();
			}
		}
		void pause() {
			std::unique_lock lock(this->mutex);
			this->changed.wait(lock, [this] { return !this->paused.load(); });
			this->paused.store(true);
			this->changed.wait(lock, [this] { return this->active.load() == 0; });
		}
		void resume() {
			{
				std::lock_guard lock(this->mutex);
				this->paused.store(false);
			}
			this->changed.notify_all();
		}
	};

	// A single app-owned object shared across all MT workers: one accumulated field, guarded by striped
	// per-voxel locks. It is a G4UserSteppingAction but is not registered with Geant4 directly (a per-worker
	// RadiationFieldSteppingAction forwarder is registered instead and calls UserSteppingAction() on it), so
	// Geant4 owns the per-worker forwarders while the app remains the sole owner of this shared detector.
	class RadiationFieldDetector: public G4UserSteppingAction {
	protected:
		// The accumulated scoring field and its spectrum layout. Declared first so it is constructed before
		// `buffers`, whose initializer calls field->add_channel(...).
		std::shared_ptr<radfiled3d::CartesianRadiationField> field;
		const size_t spectra_bins;
		const double spectra_bin_width;

		class ChannelBuffers {
		protected:
			// Striped voxel locks: a voxel always maps to the same mutex (idx % pool), locks are
			// taken one voxel at a time (never nested), so contention semantics match the old
			// one-mutex-per-voxel layout at a fraction of the memory.
			static constexpr size_t MUTEX_POOL_SIZE = 4096;
			Statistics::VoxelSpectraVariance spectra_variance;
			float total_energy = 0.f;
			std::vector<std::shared_mutex> mutexes;
			mutable std::shared_mutex buffer_mutex;
			const uint32_t directional_lobes;

			std::unique_ptr<VMFTrainer> make_trainer() const {
				if (this->directional_lobes == 0)
					return nullptr;
				// the first lobe of every voxel starts pointing away from the isocentre (the field centre)
				const glm::uvec3 counts = this->buffer.get_voxel_counts();
				const glm::vec3 size = this->buffer.get_voxel_dimensions();
				const glm::vec3 half = this->half_field_dim / static_cast<float>(m);
				return std::make_unique<VMFTrainer>(this->buffer.get_voxel_count(), this->directional_lobes, [counts, size, half](size_t voxel) {
					const glm::vec3 index(voxel % counts.x, (voxel / counts.x) % counts.y, voxel / (static_cast<size_t>(counts.x) * counts.y));
					return (index + 0.5f) * size - half;
				});
			}
		public:
			radfiled3d::VoxelGridBuffer& buffer;
			const glm::vec3 half_field_dim;
			float statistical_error_resolution = 0.5f;
			// Directions of travel per voxel as a vMF mixture, learned during the run (nullptr if disabled).
			std::unique_ptr<VMFTrainer> vmf;

			void reset() {
				std::unique_lock lock(this->buffer_mutex);
				this->total_energy = 0.f;
				this->buffer.clear_layer<double>("flux", 0.0);
				this->buffer.clear_layer<double>("spectrum", 0.0, this->buffer.get_voxel_flat<radfiled3d::HistogramVoxel<double>>("spectrum", 0).get_bins());
				if (this->buffer.has_layer("angular_flux"))
					this->buffer.clear_layer<double>("angular_flux", 0.0, this->buffer.get_voxel_flat<radfiled3d::AngularResolvedVoxel<double>>("angular_flux", 0).get_total_segments());
				this->vmf = this->make_trainer();
			}

			inline float get_total_energy() const { return this->total_energy; }

			inline int get_voxel_idx(const glm::vec3& position) {
				const glm::vec3 positive_position = (position + this->half_field_dim) / glm::vec3(m);

				if (positive_position.x <= 0.f || positive_position.y <= 0.f || positive_position.z <= 0.f)
					return -2;
				
				size_t idx = this->buffer.get_voxel_idx_by_coord(positive_position.x, positive_position.y, positive_position.z);
				if (idx >= this->buffer.get_voxel_count())
					return -1;
				return static_cast<int>(idx);
			}

			float get_overall_statistical_error_estimate(size_t primary_particle_count, float statistical_error_enforcement_ratio = 1.f) {
				assert(statistical_error_enforcement_ratio >= 0.f && statistical_error_enforcement_ratio <= 1.f);

				std::unique_lock lock(this->buffer_mutex);

				size_t step_width = static_cast<size_t>(1.f / this->statistical_error_resolution);
				if (step_width == 0)
					step_width = 1;

				glm::uvec3 steps_per_dim = glm::uvec3(
					static_cast<size_t>(this->buffer.get_voxel_counts().x) / step_width,
					static_cast<size_t>(this->buffer.get_voxel_counts().y) / step_width,
					static_cast<size_t>(this->buffer.get_voxel_counts().z) / step_width
				);

				std::vector<float> errors;

				for (size_t step_x = 0; step_x < steps_per_dim.x; step_x++) {
					const size_t x = step_x * step_width;
					for (size_t step_y = 0; step_y < steps_per_dim.y; step_y++) {
						const size_t y = step_y * step_width;
						for (size_t step_z = 0; step_z < steps_per_dim.z; step_z++) {
							const size_t z = step_z * step_width;
							errors.push_back(this->get_statistical_error(primary_particle_count, x, y, z));
						}
					}
				}

				std::sort(errors.begin(), errors.end());
				float stat_error = errors[std::max<size_t>(static_cast<size_t>(errors.size() * statistical_error_enforcement_ratio), 1) - 1];
				
				return stat_error;
			}

			/** Scores one step p1 -> p2 (grid frame, metres) into the voxels it enters, each weighted by
			* TracedVoxel::path_fraction (the path length there with the line tracer, 1 otherwise): flux, spectrum, angular
			* bins and vMF lobes all hold the incident radiation. The voxel the step starts in is left out: after an
			* interaction, the radiation there leaves the voxel instead of entering it.
			*/
			inline void score(float energy, const glm::vec3& p1, const glm::vec3& p2, const std::vector<radfiled3d::TracedVoxel>& voxels) {
				const glm::vec3 direction = p2 - p1;
				const float r = glm::length(direction);
				double weights = 0.0;
				for (const radfiled3d::TracedVoxel& traced : voxels)
					if (!traced.starts_segment)
						weights += traced.path_fraction;
				{
					std::unique_lock lock(this->buffer_mutex);
					this->total_energy += static_cast<float>(energy * weights);
				}

				for (const radfiled3d::TracedVoxel& traced : voxels) {
					if (traced.starts_segment || !(traced.path_fraction > 0.f))
						continue;
					const size_t voxel_idx = traced.index;
					const double weight = traced.path_fraction;
					auto& hist_voxel = buffer.get_voxel_flat<radfiled3d::HistogramVoxel<double>>("spectrum", voxel_idx);
					size_t index = static_cast<size_t>(energy / hist_voxel.get_histogram_bin_width());
					if (index >= hist_voxel.get_bins()) {
						index = hist_voxel.get_bins() - 1;
						G4cout << "WARNING: Energy value exceeds histogram energy range. Energy: " << energy << " MeV" << G4endl;
					}

					// The pool is sized min(voxel_count, MUTEX_POOL_SIZE), so index by the actual size, not the
					// cap — else a field with < MUTEX_POOL_SIZE voxels indexes past the vector (OOB -> crash).
					std::unique_lock lock(this->mutexes[voxel_idx % this->mutexes.size()]);
					this->buffer.get_voxel_flat<radfiled3d::ScalarVoxel<double>>("flux", voxel_idx) += weight;

					(&hist_voxel.get_data())[index] += weight;
					this->spectra_variance.add(voxel_idx, hist_voxel);

					if (r > 0.f) {
						if (this->buffer.has_layer("angular_flux")) {
							float theta = std::acos(glm::clamp(direction.z / r, -1.f, 1.f));
							float phi = std::atan2(direction.y, direction.x);
							if (phi < 0.f) phi += 2.f * std::numbers::pi_v<float>;
							this->buffer.get_voxel_flat<radfiled3d::AngularResolvedVoxel<double>>("angular_flux", voxel_idx).add_value(phi, theta, weight);
						}
						if (this->vmf)
							this->vmf->add(voxel_idx, direction / r, weight);
					}
				}
			}

			inline float get_statistical_error(size_t primary_particle_count, size_t x, size_t y, size_t z) {
				return this->spectra_variance.get_relative_error(this->buffer.get_voxel_idx(x, y, z));
			}

			ChannelBuffers(radfiled3d::VoxelGridBuffer& buffer, float spectra_bin_width, size_t spectra_bins, glm::uvec2 angular_resolution = glm::uvec2(0), uint32_t directional_lobes = 0)
				: directional_lobes(directional_lobes),
				  buffer(buffer),
				  half_field_dim(
					  glm::vec3(
						  static_cast<float>(buffer.get_voxel_counts().x * buffer.get_voxel_dimensions().x * m) / 2.f,
					      static_cast<float>(buffer.get_voxel_counts().y * buffer.get_voxel_dimensions().y * m) / 2.f,
						  static_cast<float>(buffer.get_voxel_counts().z * buffer.get_voxel_dimensions().z * m) / 2.f
					  )
				  ),
				  spectra_variance(buffer.get_voxel_count(), spectra_bins, 50),
				  mutexes(std::min(buffer.get_voxel_count(), MUTEX_POOL_SIZE))
			{
				// Scoring accumulates in DOUBLE: float += 1 saturates at 2^24 counts, which the
				// beam-entry voxels reach at ~1e8 primaries. The stored field is converted to
				// fp32 after normalization (get_normalized_field_copy).
				buffer.add_layer<float>("error", 1.f, "Variance");
				buffer.add_layer<double>("flux", 0.0, "counts / primary_particles");
				buffer.add_custom_layer<radfiled3d::HistogramVoxel<double>>("spectrum", radfiled3d::HistogramVoxel<double>(spectra_bins, static_cast<double>(spectra_bin_width), nullptr), 0.0, "eV");
				if (angular_resolution.x > 0 && angular_resolution.y > 0)
					buffer.add_custom_layer<radfiled3d::AngularResolvedVoxel<double>>("angular_flux", radfiled3d::AngularResolvedVoxel<double>(angular_resolution, nullptr), 0.0, "counts / primary_particles");
				this->vmf = this->make_trainer();
			}

			// Copying would re-run the main ctor on the SAME VoxelGridBuffer: the duplicate
			// add_layer calls allocate arrays that map::insert then silently drops (a hard leak),
			// and the per-voxel statistics would exist twice. Channels constructs in place.
			ChannelBuffers(const ChannelBuffers&) = delete;
			ChannelBuffers& operator=(const ChannelBuffers&) = delete;
		};

		size_t primary_particle_count = 0;
		mutable std::shared_mutex global_detector_mutex;

		struct Channels {
			ChannelBuffers scatter_field;
			ChannelBuffers xray_beam;

			Channels(radfiled3d::VoxelGridBuffer& scatter_buffer, radfiled3d::VoxelGridBuffer& xray_buffer, float spectra_bin_width, size_t spectra_bins, const glm::uvec2& angular_resolution, uint32_t directional_lobes)
				: scatter_field(scatter_buffer, spectra_bin_width, spectra_bins, angular_resolution, directional_lobes),
				  xray_beam(xray_buffer, spectra_bin_width, spectra_bins, angular_resolution) {}

			void reset() {
				this->scatter_field.reset();
				this->xray_beam.reset();
			}
		} buffers;

		struct EventContext {
			size_t event_id;
			std::map<size_t, TrackStage> track_stage;

			EventContext(size_t event_id) : event_id(event_id) {}
		};

		std::map<size_t, EventContext> thread_contexts;
		std::atomic<size_t> tracked_events_counter;
		// vMF training passes end at 1e6, 2e6, 4e6, ... primaries (doubling passes as in path guiding)
		static constexpr size_t VMF_FIRST_PASS = 1000000;
		std::atomic<size_t> next_vmf_pass{ VMF_FIRST_PASS };
		void train_directional_lobes();
		
		// Materials of the patient geometry: nothing is scored inside them. A handful of entries, so a linear scan is fastest.
		std::vector<const G4Material*> patient_materials;
		bool is_patient(const G4Material* material) const {
			for (const G4Material* m : this->patient_materials)
				if (m == material)
					return true;
			return false;
		}
		const float statistical_error_threshold;
		const float statistical_error_enforcement_ratio;
		bool is_tracking = true;
		const float simulation_energy_lower_threshold = 1 * keV;
		void score_step_for(const G4Step* step, const glm::vec3& p1, const glm::vec3& p2, const std::vector<radfiled3d::TracedVoxel>& voxels, TrackStage stage);
		void evaluate_field();
		std::shared_ptr<radfiled3d::GridTracer> tracer;
		// Every voxel a step enters weighted by the path length inside it (track-length estimate, line tracer only);
		// otherwise each voxel the tracer returns counts once.
		bool path_length_weighting = false;
		std::vector< std::function<void(size_t, const G4Step*)>> new_particle_callbacks;
		// Paused while a normalized copy is taken, so the copy and its primary count describe the same events.
		ScoringGate scoring_gate;
		// Voxelized scene geometry on the scoring grid, computed once per run.
		std::shared_ptr<radfiled3d::CartesianRadiationField> geometry;
	public:
		RadiationFieldDetector(
			const glm::vec3& radiation_field_dimensions,
			const glm::vec3& radiation_field_voxel_dimensions,
			size_t spectra_bins,
			double spectra_bin_width,
			float statistical_error_threshold = 0.1f,
			float statistical_error_enforcement_ratio = 0.9f,
			float statistical_error_enforcement_resolution = 0.5f,
			const glm::uvec2& angular_resolution = glm::uvec2(0),
			uint32_t directional_lobes = 0
		);
		virtual ~RadiationFieldDetector() {
			G4cout << "RadiationFieldDetector destroyed" << G4endl;
		}
		/** Sets the materials of the patient geometry (its own and its children's): nothing is scored inside them. */
		void SetUp(const std::vector<const G4Material*>& patient_materials);
		const std::vector<const G4Material*>& get_patient_materials() const { return this->patient_materials; }
		virtual void finalize(size_t particle_count);

		// Runs the field evaluation and returns the normalized fp32 copy.
		std::shared_ptr<radfiled3d::IRadiationField> evaluate();

		/** @param path_length_weighting With the line tracer, weight every voxel a step enters by the path length inside
		* it (track-length estimate) instead of counting it once; ignored by the other tracers. */
		template<class T>
		inline void define_grid_tracer(bool path_length_weighting = true) {
			static_assert(std::is_base_of<radfiled3d::GridTracer, T>::value, "T must be derived from radfiled3d::GridTracer");
			this->tracer = std::make_shared<T>(this->buffers.scatter_field.buffer);
			this->path_length_weighting = path_length_weighting && std::is_same_v<T, radfiled3d::LinetracingGridTracer>;
		}

		size_t get_number_of_tracked_particles() const;
		size_t get_primary_particle_count() const;
		
		// Scores one step into the field; invoked for every step by the per-worker forwarder.
		virtual void UserSteppingAction(const G4Step* step) override;
		std::shared_ptr<radfiled3d::IRadiationField> get_normalized_field_copy();
		/** Voxelizes the scene geometry below `world_volume` onto the scoring grid (see Geant4::add_geometry_channel). */
		void voxelize_geometry(const G4LogicalVolume& world_volume, int max_threads = -1);
		/** Adds a copy of the voxelized geometry as channel "geometry" to `field`, which must share the scoring grid.
		* Does nothing if no geometry was voxelized.
		*/
		void add_geometry_channel_to(radfiled3d::CartesianRadiationField& field) const;

		float get_statistical_error(size_t primary_particle_count = 0);
		void register_on_new_particle(std::function<void(size_t, const G4Step*)> callback);
	};

	// Lightweight per-worker stepping action owned by Geant4 (one per worker thread, created in Build()). It
	// owns nothing and routes each step into the single app-owned detector shared by all workers, so every
	// worker scores into one field.
	class RadiationFieldSteppingAction : public G4UserSteppingAction {
		RadiationFieldDetector* detector;   // non-owning: the app owns the shared detector
	public:
		explicit RadiationFieldSteppingAction(RadiationFieldDetector* detector) : detector(detector) {}
		virtual void UserSteppingAction(const G4Step* step) override { this->detector->UserSteppingAction(step); }
	};

	class RadiationFieldAction : public G4VUserActionInitialization {
	protected:
		std::shared_ptr<RadiationFieldDetector> det;
		// The physics source is shared read-only; Build() (called per worker thread) constructs a fresh
		// RadiationSource per worker from it — a single shared generator across MT workers corrupts the
		// gun's thread-local allocations (non-deterministic mid-run segfault). See the .cpp.
		std::shared_ptr<RadiationSimulation::RadiationSource> rad_source;
		int fluence_per_run;
	public:
		RadiationFieldAction(std::shared_ptr<RadiationFieldDetector> det, std::shared_ptr<RadiationSimulation::RadiationSource> rad_source, int fluence_per_run = 1) : det(det), rad_source(rad_source), fluence_per_run(fluence_per_run) {};
		void Build() const;
		virtual ~RadiationFieldAction() {
			G4cout << "RadiationFieldAction destroyed" << G4endl;
		}
	};
}