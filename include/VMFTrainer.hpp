#pragma once
#include <glm/vec3.hpp>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>


namespace RadiationSimulation {
	/** Learns, per voxel, a von Mises-Fisher mixture of the directions of travel while the simulation runs.
	* Training runs in passes as in path guiding (Müller et al. 2017; stepwise EM: Vorba et al. 2014): during a pass,
	* add() accumulates per voxel and lobe the sufficient statistics (responsibility, responsibility-weighted direction)
	* with the lobes of the previous pass; m_step() turns them into new lobes and starts the next pass.
	* The lobe count is a maximum: in the stored result, lobes that describe the same source (similar width, means
	* closer than the lobes' width) are merged and the freed lobe gets weight 0, so a voxel lit by one source ends with
	* a single lobe; a compact source and a wide background around it stay two lobes. Training itself keeps all lobes.
	* A voxel with too few scored directions to fit lobes is stored as the uniform distribution (one lobe, weight 1,
	* kappa 0): every voxel starts uniform and stays so until it has `min_samples` effective samples.
	* add() must be called under the caller's per-voxel lock, m_step() while nothing is added.
	*/
	class VMFTrainer {
	public:
		/** Values per lobe in the stored layout (radfiled3d::VMFMixtureVoxel): weight, mean x, y, z, kappa. */
		static constexpr size_t VALUES_PER_LOBE = 5;
		/** Effective samples (sum w)^2 / sum w^2 a voxel needs before its fitted lobes are stored instead of the uniform
		* distribution. Directions from n uniform samples have a mean resultant length of about 1 / sqrt(n), which a fit
		* reads as kappa of about 3 / sqrt(n): at 20 samples that spurious concentration stays below 0.7. */
		static constexpr double DEFAULT_MIN_SAMPLES = 20.0;

		/**
		* @param voxel_count Number of voxels.
		* @param lobes Lobes per voxel (at least 1).
		* @param initial_direction Unit vector for the first lobe of a voxel before any training (e.g. away from the
		*        isocentre); the other lobes start opposite and perpendicular to it, all wide. Training needs these distinct
		*        starts: identical (e.g. all uniform) lobes would take equal shares of every sample and never separate.
		* @param min_samples Effective samples a voxel needs before its fitted lobes are stored (see DEFAULT_MIN_SAMPLES).
		*/
		VMFTrainer(size_t voxel_count, uint32_t lobes, const std::function<glm::vec3(size_t)>& initial_direction, double min_samples = DEFAULT_MIN_SAMPLES);

		/** Adds one scored direction of travel (unit vector) to a voxel, with a weight (e.g. the path length in the voxel). */
		void add(size_t voxel_idx, const glm::vec3& direction, double weight = 1.0);

		/** Fits new lobes from the finished pass and starts the next one. */
		void m_step();

		/** Writes the voxel's lobes for storing (fitted from the previous and the current pass) to out[0 .. lobes * 5).
		* Always all `lobes` slots, so every field shares the same layout: sorted by weight (strongest first), unused
		* slots (weight 0) all zero. A voxel with fewer than `min_samples` effective samples in these two passes (none
		* included) gets the uniform distribution: weight 1, mean (0, 0, 1), kappa 0 in the first slot.
		*/
		void write_lobes(size_t voxel_idx, float* out) const;

		uint32_t get_lobes() const { return this->lobes; }
		size_t get_passes() const { return this->passes; }

	private:
		static constexpr size_t STATS_PER_LOBE = 4;   // sum r, sum r * direction

		const size_t voxel_count;
		const uint32_t lobes;
		const double min_samples;
		size_t passes = 0;
		std::vector<float> model;          // voxel * lobes * VALUES_PER_LOBE
		std::vector<float> log_norm;       // voxel * lobes: log(weight * C(kappa)) of the model
		std::vector<double> current;       // voxel * lobes * STATS_PER_LOBE
		std::vector<double> previous;
		std::vector<double> current_squared_weights;    // voxel: sum w^2, for the effective sample count
		std::vector<double> previous_squared_weights;

		void fit(const double* stats, const float* fallback, float* out, bool merge) const;
		void update_log_norm(size_t voxel_idx);
	};
}
