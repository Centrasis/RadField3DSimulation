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
	* add() must be called under the caller's per-voxel lock, m_step() while nothing is added.
	*/
	class VMFTrainer {
	public:
		/** Values per lobe in the stored layout (radfiled3d::VMFMixtureVoxel): weight, mean x, y, z, kappa. */
		static constexpr size_t VALUES_PER_LOBE = 5;

		/**
		* @param voxel_count Number of voxels.
		* @param lobes Lobes per voxel (1..8).
		* @param initial_direction Unit vector for the first lobe of a voxel before any training (e.g. away from the
		*        isocentre); the other lobes start opposite and perpendicular to it, all wide.
		*/
		VMFTrainer(size_t voxel_count, uint32_t lobes, const std::function<glm::vec3(size_t)>& initial_direction);

		/** Adds one scored direction of travel (unit vector) to a voxel, with a weight (e.g. the path length in the voxel). */
		void add(size_t voxel, const glm::vec3& direction, double weight = 1.0);

		/** Fits new lobes from the finished pass and starts the next one. */
		void m_step();

		/** Writes the voxel's lobes for storing (fitted from the previous and the current pass) to out[0 .. lobes * 5).
		* Always all `lobes` slots, so every field shares the same layout: sorted by weight (strongest first), unused
		* slots (weight 0) and voxels without any scored direction all zero.
		*/
		void write_lobes(size_t voxel, float* out) const;

		uint32_t get_lobes() const { return this->lobes; }
		size_t get_passes() const { return this->passes; }

	private:
		static constexpr size_t STATS_PER_LOBE = 4;   // sum r, sum r * direction

		const size_t voxel_count;
		const uint32_t lobes;
		size_t passes = 0;
		std::vector<float> model;          // voxel * lobes * VALUES_PER_LOBE
		std::vector<float> log_norm;       // voxel * lobes: log(weight * C(kappa)) of the model
		std::vector<double> current;       // voxel * lobes * STATS_PER_LOBE
		std::vector<double> previous;

		void fit(const double* stats, const float* fallback, float* out, bool merge) const;
		void update_log_norm(size_t voxel);
	};
}
