#pragma once
#include <memory>
#include <vector>
#include <RadiationFieldDetector.hpp>
#include <RadFiled3D/RadiationField.hpp>
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
#include <RadFiled3D/GridTracer.hpp>
#include <G4EmCalculator.hh>
#include <shared_mutex>
#include "Statistics.hpp"
#include <algorithm>
#include <numbers>


class G4RunManager;
#ifdef WITH_GEANT4_UIVIS
class G4UImanager;
class G4VisExecutive;
#endif
class G4LogicalVolume;

namespace RadiationSimulation {
	class Mesh;
	class G4RadiationSource;
	class RadiationSource;

	enum class TrackStage : char {
		BEAM,
		SCATTER,
		PATIENT
	};

	// A single app-owned object shared across all MT workers: one accumulated field, guarded by striped
	// per-voxel locks. It is a G4UserSteppingAction but is not registered with Geant4 directly (a per-worker
	// G4RadiationFieldSteppingAction forwarder is registered instead and calls UserSteppingAction() on it), so
	// Geant4 owns the per-worker forwarders while the app remains the sole owner of this shared detector.
	class G4RadiationFieldDetector: public G4UserSteppingAction, public RadiationFieldDetector {
	protected:
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
		public:
			RadFiled3D::VoxelGridBuffer& buffer;
			const glm::vec3 half_field_dim;
			float statistical_error_resolution = 0.5f;

			void reset() {
				std::unique_lock lock(this->buffer_mutex);
				this->total_energy = 0.f;
				this->buffer.clear_layer<double>("flux", 0.0);
				this->buffer.clear_layer<double>("spectrum", 0.0, this->buffer.get_voxel_flat<RadFiled3D::HistogramVoxel<double>>("spectrum", 0).get_bins());
				if (this->buffer.has_layer("angular_flux"))
					this->buffer.clear_layer<double>("angular_flux", 0.0, this->buffer.get_voxel_flat<RadFiled3D::AngularResolvedVoxel<double>>("angular_flux", 0).get_total_segments());
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

			inline void score(float energy, const glm::vec3& direction, const std::vector<size_t>& voxel_indices) {
				{
					std::unique_lock lock(this->buffer_mutex);
					this->total_energy += energy * voxel_indices.size();
				}

				for (size_t voxel_idx : voxel_indices) {
					auto& hist_voxel = buffer.get_voxel_flat<RadFiled3D::HistogramVoxel<double>>("spectrum", voxel_idx);
					size_t index = static_cast<size_t>(energy / hist_voxel.get_histogram_bin_width());
					if (index >= hist_voxel.get_bins()) {
						index = hist_voxel.get_bins() - 1;
						G4cout << "WARNING: Energy value exceeds histogram energy range. Energy: " << energy << " MeV" << G4endl;
					}

					// The pool is sized min(voxel_count, MUTEX_POOL_SIZE), so index by the actual size, not the
				// cap — else a field with < MUTEX_POOL_SIZE voxels indexes past the vector (OOB -> crash).
				std::unique_lock lock(this->mutexes[voxel_idx % this->mutexes.size()]);
					this->buffer.get_voxel_flat<RadFiled3D::ScalarVoxel<double>>("flux", voxel_idx) += 1.0;

					(&hist_voxel.get_data())[index] += 1.0;
					this->spectra_variance.add(voxel_idx, hist_voxel);

					if (this->buffer.has_layer("angular_flux")) {
						float r = glm::length(direction);
						if (r > 0.f) {
							float theta = std::acos(glm::clamp(direction.z / r, -1.f, 1.f));
							float phi = std::atan2(direction.y, direction.x);
							if (phi < 0.f) phi += 2.f * std::numbers::pi_v<float>;
							this->buffer.get_voxel_flat<RadFiled3D::AngularResolvedVoxel<double>>("angular_flux", voxel_idx).add_value(phi, theta, 1.0);
						}
					}
				}
			}

			inline float get_statistical_error(size_t primary_particle_count, size_t x, size_t y, size_t z) {
				return this->spectra_variance.get_relative_error(this->buffer.get_voxel_idx(x, y, z));
			}

			ChannelBuffers(RadFiled3D::VoxelGridBuffer& buffer, float spectra_bin_width, size_t spectra_bins, glm::uvec2 angular_resolution = glm::uvec2(0))
				: buffer(buffer),
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
				buffer.add_custom_layer<RadFiled3D::HistogramVoxel<double>>("spectrum", RadFiled3D::HistogramVoxel<double>(spectra_bins, static_cast<double>(spectra_bin_width), nullptr), 0.0, "eV");
				if (angular_resolution.x > 0 && angular_resolution.y > 0)
					buffer.add_custom_layer<RadFiled3D::AngularResolvedVoxel<double>>("angular_flux", RadFiled3D::AngularResolvedVoxel<double>(angular_resolution, nullptr), 0.0, "counts / primary_particles");
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

			Channels(RadFiled3D::VoxelGridBuffer& scatter_buffer, RadFiled3D::VoxelGridBuffer& xray_buffer, float spectra_bin_width, size_t spectra_bins, const glm::uvec2& angular_resolution)
				: scatter_field(scatter_buffer, spectra_bin_width, spectra_bins, angular_resolution),
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
		
		G4Material* air_material = NULL;
		const float statistical_error_threshold;
		const float statistical_error_enforcement_ratio;
		bool is_tracking = true;
		const float simulation_energy_lower_threshold = 1 * keV;
		void score_step_for(const G4Step* step, const std::vector<size_t>& voxel_indices, TrackStage stage);
		virtual void evaluate_field() override;
		std::shared_ptr<RadFiled3D::GridTracer> tracer;
		std::vector< std::function<void(size_t, const G4Step*)>> new_particle_callbacks;
	public:
		G4RadiationFieldDetector(
			const glm::vec3& radiation_field_dimensions,
			const glm::vec3& radiation_field_voxel_dimensions,
			size_t spectra_bins,
			double spectra_bin_width,
			float statistical_error_threshold = 0.1f,
			float statistical_error_enforcement_ratio = 0.9f,
			float statistical_error_enforcement_resolution = 0.5f,
			const glm::uvec2& angular_resolution = glm::uvec2(0)
		);
		virtual ~G4RadiationFieldDetector() {
			G4cout << "G4RadiationFieldDetector destroyed" << G4endl;
		}
		virtual void SetUp() override;
		virtual void finalize(size_t particle_count);

		template<class T>
		inline void define_grid_tracer() {
			static_assert(std::is_base_of<RadFiled3D::GridTracer, T>::value, "T must be derived from RadFiled3D::GridTracer");
			this->tracer = std::make_shared<T>(this->buffers.scatter_field.buffer);
		}

		virtual size_t get_number_of_tracked_particles() const override;
		virtual size_t get_primary_particle_count() const override;
		
		// Scores one step into the field; invoked for every step by the per-worker forwarder.
		virtual void UserSteppingAction(const G4Step* step) override;
		virtual std::shared_ptr<RadFiled3D::IRadiationField> get_normalized_field_copy() override;

		virtual float get_statistical_error(size_t primary_particle_count = 0) override;
		void register_on_new_particle(std::function<void(size_t, const G4Step*)> callback);
	};

	// Lightweight per-worker stepping action owned by Geant4 (one per worker thread, created in Build()). It
	// owns nothing and routes each step into the single app-owned detector shared by all workers, so every
	// worker scores into one field.
	class G4RadiationFieldSteppingAction : public G4UserSteppingAction {
		G4RadiationFieldDetector* detector;   // non-owning: the app owns the shared detector
	public:
		explicit G4RadiationFieldSteppingAction(G4RadiationFieldDetector* detector) : detector(detector) {}
		virtual void UserSteppingAction(const G4Step* step) override { this->detector->UserSteppingAction(step); }
	};

	class G4RadiationFieldAction : public G4VUserActionInitialization {
	protected:
		std::shared_ptr<G4RadiationFieldDetector> det;
		// The physics source is shared read-only; Build() (called per worker thread) constructs a fresh
		// G4RadiationSource per worker from it — a single shared generator across MT workers corrupts the
		// gun's thread-local allocations (non-deterministic mid-run segfault). See the .cpp.
		std::shared_ptr<RadiationSource> rad_source;
		int fluence_per_run;
	public:
		G4RadiationFieldAction(std::shared_ptr<G4RadiationFieldDetector> det, std::shared_ptr<RadiationSource> rad_source, int fluence_per_run = 1) : det(det), rad_source(rad_source), fluence_per_run(fluence_per_run) {};
		void Build() const;
		virtual ~G4RadiationFieldAction() {
			G4cout << "G4RadiationFieldAction destroyed" << G4endl;
		}
	};
}