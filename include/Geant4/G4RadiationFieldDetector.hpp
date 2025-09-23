#pragma once
#include <memory>
#include <vector>
#include <RadiationFieldDetector.hpp>
#include <RadFiled3D/RadiationField.hpp>
#include <G4VSensitiveDetector.hh>
#include <glm/vec3.hpp>
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


class G4RunManager;
#ifdef WITH_GEANT4_UIVIS
class G4UImanager;
class G4VisExecutive;
#endif
class G4LogicalVolume;

namespace RadiationSimulation {
	class Mesh;
	class G4RadiationSource;

	enum class TrackStage : char {
		BEAM,
		SCATTER,
		PATIENT
	};

	class G4RadiationFieldDetector: public G4UserSteppingAction, public RadiationFieldDetector {
	protected:
		class ChannelBuffers {
		protected:
			std::vector<OnlinePCA> pca;
			std::vector<Statistics::HistogramDistributionVariance> spectra_variance;
			float total_energy = 0.f;
			std::vector<std::shared_mutex> mutexes;
			mutable std::shared_mutex buffer_mutex;
			float statistical_error_resolution = 0.5f;
		public:
			RadFiled3D::VoxelGridBuffer& buffer;
			const glm::vec3 half_field_dim;

			void reset() {
				std::unique_lock lock(this->buffer_mutex);
				this->total_energy = 0.f;
				for (auto& pca : this->pca)
					pca.reset();
				/*for (auto& var : this->spectra_variance)
					var.reset();*/
				this->buffer.clear_layer<float>("energy", 0.f);
				this->buffer.clear_layer<float>("hits", 0.f);
				this->buffer.clear_layer<glm::vec3>("direction", glm::vec3(0.f));
				this->buffer.clear_layer<float>("spectrum", 0.f, this->buffer.get_voxel_flat<RadFiled3D::HistogramVoxel>("spectrum", 0).get_bins());
			}

			inline float get_total_energy() const { return this->total_energy; }
			inline float get_energy(size_t x, size_t y, size_t z) const { return this->buffer.get_voxel<RadFiled3D::ScalarVoxel<float>>("energy", x, y, z).get_data(); }

			inline int get_voxel_idx(const glm::vec3& position) {
				const glm::vec3 positive_position = (position + this->half_field_dim) / glm::vec3(m);

				if (positive_position.x <= 0.f || positive_position.y <= 0.f || positive_position.z <= 0.f)
					return -2;
				
				size_t idx = this->buffer.get_voxel_idx_by_coord(positive_position.x, positive_position.y, positive_position.z);
				if (idx >= pca.size())
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

			inline const std::vector<OnlinePCA>& get_pcas() const { return this->pca; }

			inline void score(float energy, const glm::vec3& direction, const std::vector<size_t>& voxel_indices) {
				{
					std::unique_lock lock(this->buffer_mutex);
					this->total_energy += energy * voxel_indices.size();
				}

				for (size_t voxel_idx : voxel_indices) {
					auto& hist_voxel = buffer.get_voxel_flat<RadFiled3D::HistogramVoxel>("spectrum", voxel_idx);
					size_t index = static_cast<size_t>(energy / hist_voxel.get_histogram_bin_width());
					if (index >= hist_voxel.get_bins()) {
						index = hist_voxel.get_bins() - 1;
						G4cout << "WARNING: Energy value exceeds histogram energy range. Energy: " << energy << " MeV" << G4endl;
					}

					std::unique_lock lock(this->mutexes[voxel_idx]);
					this->buffer.get_voxel_flat<RadFiled3D::ScalarVoxel<float>>("energy", voxel_idx) += energy;
					this->buffer.get_voxel_flat<RadFiled3D::ScalarVoxel<float>>("hits", voxel_idx) += 1.f;

					this->pca[voxel_idx].addVector(direction);

					(&hist_voxel.get_data())[index] += 1.f;
					this->spectra_variance[voxel_idx].add(hist_voxel);
				}
			}

			inline float get_statistical_error(size_t primary_particle_count, size_t x, size_t y, size_t z) {
				return this->spectra_variance[this->buffer.get_voxel_idx(x, y, z)].get_relative_error();
			}

			ChannelBuffers(RadFiled3D::VoxelGridBuffer& buffer, float spectra_bin_width, size_t spectra_bins) 
				: buffer(buffer),
				  half_field_dim(
					  glm::vec3(
						  static_cast<float>(buffer.get_voxel_counts().x * buffer.get_voxel_dimensions().x * m) / 2.f,
					      static_cast<float>(buffer.get_voxel_counts().y * buffer.get_voxel_dimensions().y * m) / 2.f,
						  static_cast<float>(buffer.get_voxel_counts().z * buffer.get_voxel_dimensions().z * m) / 2.f
					  )
				  ),
				  pca(buffer.get_voxel_count()),
				  spectra_variance(buffer.get_voxel_count(), Statistics::HistogramDistributionVariance(spectra_bins, 50)),
				  mutexes(buffer.get_voxel_count())
			{
				buffer.add_layer<float>("energy", 0.f, "eV");
				buffer.add_layer<float>("error", 1.f, "Variance");
				buffer.add_layer<float>("hits", 0.f, "counts");
				buffer.add_layer<glm::vec3>("direction", glm::vec3(0.f), "direction vector");
				buffer.add_custom_layer<RadFiled3D::HistogramVoxel>("spectrum", RadFiled3D::HistogramVoxel(spectra_bins, spectra_bin_width, nullptr), 0.f, "eV");
			}

			ChannelBuffers(const ChannelBuffers& other) 
				: ChannelBuffers(
					other.buffer,
					other.buffer.get_voxel_flat<RadFiled3D::HistogramVoxel>("spectrum", 0).get_histogram_bin_width(),
					other.buffer.get_voxel_flat<RadFiled3D::HistogramVoxel>("spectrum", 0).get_bins()
				) {}
		};

		size_t primary_particle_count = 0;
		mutable std::shared_mutex global_detector_mutex;

		struct Channels {
			ChannelBuffers scatter_field;
			ChannelBuffers xray_beam;

			Channels(const ChannelBuffers& scatter_field, const ChannelBuffers& xray_beam): scatter_field(scatter_field), xray_beam(xray_beam) {}
			
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
		G4RadiationFieldDetector(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, size_t spectra_bins, double spectra_bin_width, float statistical_error_threshold = 0.1f, float statistical_error_enforcement_ratio = 0.9f);
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
		
		virtual void UserSteppingAction(const G4Step*) override;
		virtual std::shared_ptr<RadFiled3D::IRadiationField> get_normalized_field_copy() override;

		virtual float get_statistical_error(size_t primary_particle_count = 0) override;
		void register_on_new_particle(std::function<void(size_t, const G4Step*)> callback);
	};

	class G4RadiationFieldAction : public G4VUserActionInitialization {
	protected:
		std::shared_ptr<G4RadiationFieldDetector> det;
		std::shared_ptr<G4RadiationSource> source;
	public:
		G4RadiationFieldAction(std::shared_ptr<G4RadiationFieldDetector> det, std::shared_ptr<G4RadiationSource> source) : det(det), source(source) {};
		void Build() const;
		virtual ~G4RadiationFieldAction() {
			G4cout << "G4RadiationFieldAction destroyed" << G4endl;
		}
	};
}