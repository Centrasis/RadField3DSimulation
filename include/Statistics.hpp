#pragma once
#include <vector>
#include <cstdint>
#include <atomic>
#include <unordered_map>
#include "radfiled3d/voxel.hpp"

namespace Statistics {
	/** @brief Class to calculate the variance of a set of values.
	* Variance is calculated incrementally to avoid storing all values.
	* The algorithm is based on shifted data to avoid numerical instability.
	* The algorithm is based on the following paper:
	*
	* Welford, B. P. (1962). Note on a method for calculating corrected sums of squares and products. Technometrics, 4(3), 419-420.
	*/
	class Variance {
	protected:
		float mean;
		float M2;
		size_t count = 0;
	public:
		Variance() {
			this->reset();
		};

		void add(float value);
		float get_variance() const;
		float get_mean() const;

		inline size_t get_count() const {
			return this->count;
		}

		virtual void reset() {
			this->mean = 0.f;
			this->M2 = 0.f;
			this->count = 0;
		}

		float get_relative_error() const;
		float get_standard_error() const;
	};

	/** @brief Class to calculate the mean, distribution and variance of a histogram.
	 * The histogram is assumed to be a probability distribution.
	 */
	class HistogramMeanDistributionVariance {
	protected:
		float bin_width = 0.f;
		size_t bins;
		size_t count = 0;
		size_t score_every_n = 100;
		size_t add_count = 0;
		std::vector<float> cumulative_histogram;
		float cumulative_error;
	public:
		HistogramMeanDistributionVariance(float bin_width, size_t bins, size_t score_every_n = 100);
		void add(const radfiled3d::HistogramVoxel<float>& vx);
		void reset();

		float get_relative_error() const;
	};

	/** @brief Class to calculate the variance of a histogram.
	 * The histogram is assumed to be a probability distribution.
	 */
	class HistogramDistributionVariance {
	protected:
		size_t bins;
		size_t count = 0;
		size_t score_every_n = 100;
		size_t add_count = 0;
		std::vector<Variance> variances;
	public:
		HistogramDistributionVariance(size_t bins, size_t score_every_n = 100);
		void add(const radfiled3d::HistogramVoxel<float>& vx);
		void reset();

		float get_variance() const;

		float get_relative_error() const;
	};

	/** @brief Flat-array equivalent of one HistogramDistributionVariance PER VOXEL.
	 * Same math (subsampled per-bin Welford over the normalized spectrum), but four contiguous
	 * allocations instead of millions of per-voxel objects with their own heap vectors:
	 * bins * 2 floats + 2 counters per voxel (~264 B/voxel at 32 bins vs ~840 B fragmented).
	 */
	class VoxelSpectraVariance {
	protected:
		size_t bins;
		size_t voxel_count;
		size_t score_every_n;
		std::vector<uint32_t> add_counts;
		std::vector<uint32_t> counts;
		std::vector<float> means;   // voxel-major: [voxel_idx * bins + bin]
		std::vector<float> m2s;
	public:
		VoxelSpectraVariance(size_t voxel_count, size_t bins, size_t score_every_n = 100);
		void add(size_t voxel_idx, const radfiled3d::HistogramVoxel<double>& vx);
		void reset();

		float get_relative_error(size_t voxel_idx) const;
	};

	/** @brief Per-voxel statistical error of a score per primary particle by the history-by-history method.
	* The primaries are the independent samples: x_i is everything one primary and all its secondaries scored in a
	* voxel. Squaring single steps instead would treat correlated contributions of one history (a track re-entering
	* the voxel, a step ending inside it, scattered photons coming back, secondaries) as independent and
	* underestimate the error. With S1 = sum x_i (the scored flux itself), S2 = sum x_i^2 and N primaries, the
	* standard error of the mean per primary is
	*   s = sqrt((S2 / N - (S1 / N)^2) / (N - 1)),
	* the relative error R = s / (S1 / N) = sqrt((N S2 / S1^2 - 1) / (N - 1)) ~ sqrt(S2 / S1^2 - 1 / N).
	* Only S2 is kept here: a history's per-voxel totals are collected in a History (one per worker thread) and
	* squared when it ends.
	*
	* Walters, Kawrakow, Rogers (2002). History by history statistical estimators in the BEAM code system.
	*   Med. Phys. 29(12), 2745-2752.
	* Chetty et al. (2007). Report of the AAPM Task Group No. 105, Med. Phys. 34(12), 4818-4853, eq. (3b).
	* X-5 Monte Carlo Team (2003). MCNP - A General Monte Carlo N-Particle Transport Code, Version 5, LA-UR-03-1987,
	*   Vol. I, ch. 2 (relative error R, R < 0.1 for a reliable tally).
	*/
	class VoxelHistoryVariance {
	protected:
		std::vector<std::atomic<double>> sum_of_squares;
	public:
		/** Per-voxel totals of one history that is still being simulated. */
		using History = std::unordered_map<size_t, double>;

		explicit VoxelHistoryVariance(size_t voxel_count);
		void reset();
		/** Adds the squared per-voxel totals of a finished history and clears it. Thread-safe. */
		void end_history(History& history);
		double get_sum_of_squares(size_t voxel_idx) const { return this->sum_of_squares[voxel_idx].load(std::memory_order_relaxed); }
		/** Relative standard error R of the mean score per primary from S1 = `sum`, S2 = `sum_of_squares` and
		* N = `histories`. 1 where nothing was scored or N < 2. */
		static double relative_error(double sum, double sum_of_squares, size_t histories);
	};
};