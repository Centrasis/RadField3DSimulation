#pragma once
#include <vector>
#include "RadFiled3D/Voxel.hpp"

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
		void add(const RadFiled3D::HistogramVoxel& vx);
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
		void add(const RadFiled3D::HistogramVoxel& vx);
		void reset();

		float get_variance() const;

		float get_relative_error() const;
	};
};