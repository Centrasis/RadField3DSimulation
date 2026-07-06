#include "Geant4/G4RadiationFieldDetector.hpp"
#include "RadFiled3D/RadiationField.hpp"
#include "RadiationFieldDetector.hpp"
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4SystemOfUnits.hh>
#include <G4PVPlacement.hh>
#include "Geant4/G4World.hpp"
#include <G4SDManager.hh>
#include <G4ios.hh>
#include <G4SteppingManager.hh>
#include "G4Track.hh"
#include "G4Step.hh"
#include <Geant4/G4RadiationSource.hpp>
#include <G4RunManager.hh>
#include "World.hpp"
#include "Geometry.hpp"
#include "G4NistManager.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Threading.hh"

#ifdef WITH_GEANT4_UIVIS
// visualization stuff
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Point3D.hh"
#endif


using namespace RadiationSimulation;

static std::string AIR_NAME = "G4_AIR";

void RadiationSimulation::G4RadiationFieldDetector::evaluate_field()
{
	// Diagnostics only: the double accumulators stay untouched; normalization and the fp32
	// conversion happen in get_normalized_field_copy(), which every store path uses.
	const size_t primary_particles = this->tracked_events_counter;
	std::shared_ptr<RadFiled3D::VoxelGridBuffer> scatter_buffer = this->field->get_channel("scatter_field");
	// iterate the BUFFER's voxel count: the field's count can exceed it by one per axis for
	// dimension ratios just below an integer (FLT_EPSILON asymmetry in RadFiled3D) — indexing
	// by the field count would read/write past the layer allocations.
	const size_t n = scatter_buffer->get_voxel_count();
	double total_hits = 0.0;
	double max_hits = 0.0;
	for (size_t i = 0; i < n; i++) {
		const double hits = scatter_buffer->get_voxel_flat<RadFiled3D::ScalarVoxel<double>>("flux", i).get_data();
		total_hits += hits;
		if (hits > max_hits)
			max_hits = hits;
	}

	if (max_hits == 0.0)
		G4cout << "WARNING: No voxel was hit by any particle. This is an indication of an unmatching tracking volume size or errorneous scene definitions." << G4endl;

	const double accumulated_hits_per_particle = (primary_particles > 0) ? total_hits / static_cast<double>(primary_particles) : 0.0;
	if (accumulated_hits_per_particle < 1.0)
		G4cout << "WARNING: On average there wasn't at least one voxel hit per particle. This is an indication of an unmatching tracking volume size or errorneous scene definitions. Average hits per voxel was: " << accumulated_hits_per_particle << G4endl;
}

std::shared_ptr<RadFiled3D::IRadiationField> RadiationSimulation::G4RadiationFieldDetector::get_normalized_field_copy()
{
	// Builds the STORED field: per-primary normalization of the double accumulators, written
	// into fresh fp32 layers. The scoring field itself stays untouched.
	const size_t primary_particles = this->tracked_events_counter;
	const double norm = (primary_particles > 0) ? 1.0 / static_cast<double>(primary_particles) : 0.0;
	auto copy = std::make_shared<RadFiled3D::CartesianRadiationField>(this->field->get_field_dimensions(), this->field->get_voxel_dimensions());

	struct ChannelPair { const char* name; ChannelBuffers& src; };
	ChannelPair channels[2] = { {"scatter_field", this->buffers.scatter_field}, {"direct_beam", this->buffers.xray_beam} };

	for (auto& [name, src] : channels) {
		RadFiled3D::VoxelGridBuffer& in = src.buffer;
		// iterate the BUFFER's voxel count, not the field's (FLT_EPSILON asymmetry — see evaluate_field)
		const size_t n = in.get_voxel_count();
		auto& hist0 = in.get_voxel_flat<RadFiled3D::HistogramVoxel<double>>("spectrum", 0);
		const size_t bins = hist0.get_bins();
		const bool has_angular = in.has_layer("angular_flux");

		RadFiled3D::VoxelGridBuffer* out = static_cast<RadFiled3D::VoxelGridBuffer*>(copy->add_channel(name).get());
		out->add_layer<float>("error", 1.f, "Variance");
		out->add_layer<float>("flux", 0.f, "counts / primary_particles");
		out->add_custom_layer<RadFiled3D::HistogramVoxel<float>>("spectrum", RadFiled3D::HistogramVoxel<float>(bins, static_cast<float>(hist0.get_histogram_bin_width()), nullptr), 0.f, "eV");
		if (has_angular)
			out->add_custom_layer<RadFiled3D::AngularResolvedVoxel<float>>("angular_flux", RadFiled3D::AngularResolvedVoxel<float>(in.get_voxel_flat<RadFiled3D::AngularResolvedVoxel<double>>("angular_flux", 0).get_segments(), nullptr), 0.f, "counts / primary_particles");

		for (size_t i = 0; i < n; i++) {
			out->get_voxel_flat<RadFiled3D::ScalarVoxel<float>>("flux", i) = static_cast<float>(in.get_voxel_flat<RadFiled3D::ScalarVoxel<double>>("flux", i).get_data() * norm);

			const double* h_in = &in.get_voxel_flat<RadFiled3D::HistogramVoxel<double>>("spectrum", i).get_data();
			float* h_out = &out->get_voxel_flat<RadFiled3D::HistogramVoxel<float>>("spectrum", i).get_data();
			double sum = 0.0;
			for (size_t b = 0; b < bins; b++)
				sum += h_in[b];
			if (sum > 0.0)
				for (size_t b = 0; b < bins; b++)
					h_out[b] = static_cast<float>(h_in[b] / sum);

			if (has_angular) {
				auto& a_in = in.get_voxel_flat<RadFiled3D::AngularResolvedVoxel<double>>("angular_flux", i);
				const size_t segments = a_in.get_total_segments();
				const double* s_in = &a_in.get_data();
				float* s_out = &out->get_voxel_flat<RadFiled3D::AngularResolvedVoxel<float>>("angular_flux", i).get_data();
				for (size_t s = 0; s < segments; s++)
					s_out[s] = static_cast<float>(s_in[s] * norm);
			}
		}

		for (size_t x = 0; x < out->get_voxel_counts().x; x++)
			for (size_t y = 0; y < out->get_voxel_counts().y; y++)
				for (size_t z = 0; z < out->get_voxel_counts().z; z++)
					out->get_voxel<RadFiled3D::ScalarVoxel<float>>("error", x, y, z) = src.get_statistical_error(primary_particles, x, y, z);

		out->set_statistical_error("spectrum", src.get_overall_statistical_error_estimate(primary_particles));
	}

	return copy;
}

void RadiationSimulation::G4RadiationFieldDetector::score_step_for(const G4Step* step, const std::vector<size_t>& voxel_indices, TrackStage stage)
{
	const float energy = static_cast<float>(step->GetTrack()->GetTotalEnergy());
	auto p1 = step->GetPreStepPoint()->GetPosition();
	auto p2 = step->GetPostStepPoint()->GetPosition();
	auto track_line = Collisions::Line(
		glm::vec3(p1.getX(), p1.getY(), p1.getZ()),
		glm::vec3(p2.getX(), p2.getY(), p2.getZ())
	);
	const glm::vec3 direction = track_line.direction();

	if (stage == TrackStage::SCATTER) {
		this->buffers.scatter_field.score(energy, direction, voxel_indices);
	} else if (stage == TrackStage::BEAM) {
		this->buffers.xray_beam.score(energy, direction, voxel_indices);
	}
}

RadiationSimulation::G4RadiationFieldDetector::G4RadiationFieldDetector(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, size_t spectra_bins, double spectra_bin_width, float statistical_error_threshold, float statistical_error_enforcement_ratio, float statistical_error_enforcement_resolution, const glm::uvec2& angular_resolution)
	: RadiationFieldDetector(radiation_field_dimensions, radiation_field_voxel_dimensions, spectra_bins, spectra_bin_width),
	  statistical_error_threshold(statistical_error_threshold),
	  statistical_error_enforcement_ratio(statistical_error_enforcement_ratio),
	  buffers(
	      *static_cast<RadFiled3D::VoxelGridBuffer*>(field->add_channel("scatter_field").get()),
		  *static_cast<RadFiled3D::VoxelGridBuffer*>(field->add_channel("direct_beam").get()),
		  static_cast<float>(spectra_bin_width), spectra_bins, angular_resolution
	  ),
	  G4UserSteppingAction()
{
	this->define_grid_tracer<RadFiled3D::SamplingGridTracer>();
	this->tracked_events_counter = 0;
	this->buffers.scatter_field.statistical_error_resolution = statistical_error_enforcement_resolution;
}

void RadiationSimulation::G4RadiationFieldDetector::SetUp()
{
	std::unique_lock lock(this->global_detector_mutex);
	if (this->air_material == NULL) {
		G4NistManager* nist = G4NistManager::Instance();
		this->air_material = nist->FindOrBuildMaterial(AIR_NAME);
	}
}

void RadiationSimulation::G4RadiationFieldDetector::finalize(size_t particle_count)
{
	this->primary_particle_count = particle_count;
	if (this->thread_contexts.size() > 0)
		this->thread_contexts.clear();
	this->buffers.reset();
	this->tracked_events_counter = 0;
	this->is_tracking = true;
}

size_t RadiationSimulation::G4RadiationFieldDetector::get_primary_particle_count() const
{
	return this->primary_particle_count;
}

size_t RadiationSimulation::G4RadiationFieldDetector::get_number_of_tracked_particles() const
{
	return this->tracked_events_counter;
}

void RadiationSimulation::G4RadiationFieldDetector::UserSteppingAction(const G4Step* step)
{
	if (!this->is_tracking) {
		G4RunManager::GetRunManager()->AbortEvent();
		G4RunManager::GetRunManager()->AbortRun(false);
		return;
	}

	if (this->air_material == NULL) {
		std::unique_lock lock(this->global_detector_mutex);
		if (this->air_material == NULL) {
			this->air_material = G4NistManager::Instance()->FindOrBuildMaterial(AIR_NAME);
		}
	}

	const G4int thread_id = G4Threading::G4GetThreadId();
	const size_t event_id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

	// The read must be locked too: another worker's first-event insert rebalances the map while
	// this thread traverses it. std::map iterators stay valid across inserts, so holding the
	// found iterator after releasing the lock is safe.
	std::map<size_t, EventContext>::iterator thread_context_itr;
	{
		std::shared_lock read_lock(this->global_detector_mutex);
		thread_context_itr = this->thread_contexts.find(thread_id);
	}
	if (thread_context_itr == this->thread_contexts.end()) {
		std::unique_lock lock(this->global_detector_mutex);
		this->tracked_events_counter++;
		thread_context_itr = this->thread_contexts.find(thread_id);
		if (thread_context_itr == this->thread_contexts.end()) {
			this->thread_contexts.insert({ thread_id, EventContext(event_id) });
			thread_context_itr = this->thread_contexts.find(thread_id);
		}
	}

	if (thread_context_itr->second.event_id != event_id) {
		thread_context_itr->second = EventContext(event_id);
		this->tracked_events_counter++;
		
		try {
			for (auto& cb : this->new_particle_callbacks)
				cb(this->tracked_events_counter, step);
		}
		catch (const std::exception& e) {
			G4cout << "Error in new particle callback: " << e.what() << G4endl;
		}
	}

	EventContext& event_context = thread_context_itr->second;
	const size_t track_id = step->GetTrack()->GetTrackID();

	std::map<size_t, TrackStage>::iterator current_track_stage_itr = event_context.track_stage.find(track_id);
	if (current_track_stage_itr == event_context.track_stage.end()){
		if (step->GetTrack()->GetParentID() > 0) {
			event_context.track_stage.insert({ track_id, event_context.track_stage.find(step->GetTrack()->GetParentID())->second });
		} else {
			event_context.track_stage.insert({ track_id, TrackStage::BEAM });
		}
		current_track_stage_itr = event_context.track_stage.find(track_id);
	}
	TrackStage& current_track_stage = current_track_stage_itr->second;

	// keep track of the directoion history of this step
	auto pre_step_point = step->GetPreStepPoint();
	auto post_step_point = step->GetPostStepPoint();
	auto p1 = pre_step_point->GetPosition();
	auto p2 = post_step_point->GetPosition();
	const Collisions::Line track_line = Collisions::Line(
		glm::vec3(p1.getX(), p1.getY(), p1.getZ()),
		glm::vec3(p2.getX(), p2.getY(), p2.getZ())
	);
	auto pre_mat = pre_step_point->GetMaterial();
	auto post_mat = post_step_point->GetMaterial();

	if (pre_mat != NULL && post_mat != NULL) {
		if (pre_mat != post_mat && current_track_stage == TrackStage::SCATTER) {
			if (post_mat != this->air_material) {
				// if particle was just shortly leaving the model
				current_track_stage = TrackStage::PATIENT;
			}
		}
	}

	if (current_track_stage == TrackStage::PATIENT && (post_mat == NULL || post_mat == this->air_material || (pre_mat != NULL && pre_mat == this->air_material))) {
		current_track_stage = TrackStage::SCATTER;
	}

#ifdef WITH_GEANT4_UIVIS
	G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
	if (visManager)	{
		G4Polyline polyline;
		polyline.push_back(G4Point3D(track_line.p1.x, track_line.p1.y, track_line.p1.z));
		polyline.push_back(G4Point3D(track_line.p2.x, track_line.p2.y, track_line.p2.z));
		G4Colour color(0.0, 0.0, 0.0);

		switch (current_track_stage) {
		case TrackStage::BEAM:
			color.SetRed(1.0);
			break;
		case TrackStage::PATIENT:
			color.SetGreen(1.0);
			break;
		case TrackStage::SCATTER:
			color.SetBlue(1.0);
			break;
		}
		G4VisAttributes attributes(color);
		attributes.SetLineWidth(2.0);
		polyline.SetVisAttributes(attributes);
		visManager->Draw(polyline);
	}
#endif
	
	if (current_track_stage != TrackStage::PATIENT) {
		if (step->GetTrack()->GetParticleDefinition() == G4Gamma::Definition() && step->GetTrack()->GetTotalEnergy() > 0.0) {
			const std::vector<size_t> voxel_indices = this->tracer->trace(
				(track_line.p1 + this->buffers.scatter_field.half_field_dim) / glm::vec3(m),
				(track_line.p2 + this->buffers.scatter_field.half_field_dim) / glm::vec3(m)
			);

			this->score_step_for(step, voxel_indices, current_track_stage);
		}
	}
	
	if (pre_mat != NULL && post_mat != NULL) {
		if (pre_mat != post_mat) {
			if (pre_mat == this->air_material && current_track_stage != TrackStage::PATIENT) {
				// if particle is entering the model
				current_track_stage = TrackStage::PATIENT;
			}
			else if (post_mat == this->air_material && current_track_stage != TrackStage::SCATTER) {
				// if particle is leaving the model
				current_track_stage = TrackStage::SCATTER;
			}
		}
	}

	if (this->statistical_error_threshold > 0.f && this->tracked_events_counter % 50000 == 0 && current_track_stage == TrackStage::SCATTER) {
		// check if this was enough to complete the whole simulation
		const float stat_error = this->get_statistical_error(this->tracked_events_counter);

		if (stat_error < this->statistical_error_threshold) {
			this->is_tracking = false;
			G4cout << "Reached relative error threshold: " << stat_error << " < " << this->statistical_error_threshold << G4endl;
			G4cout << "Aborting!" << G4endl;
			G4RunManager::GetRunManager()->AbortEvent();
			G4RunManager::GetRunManager()->AbortRun(false);
		}
	}
}

float RadiationSimulation::G4RadiationFieldDetector::get_statistical_error(size_t primary_particle_count)
{
	return this->buffers.scatter_field.get_overall_statistical_error_estimate(primary_particle_count, this->statistical_error_enforcement_ratio);
}

void RadiationSimulation::G4RadiationFieldDetector::register_on_new_particle(std::function<void(size_t, const G4Step*)> callback)
{
	this->new_particle_callbacks.push_back(callback);
}

void RadiationSimulation::G4RadiationFieldAction::Build() const
{
	SetUserAction(this->source.get());
	SetUserAction(this->det.get());
}
