#pragma once
#include <random>
#include "StatisticsHelpers.hpp"
#include <functional>
#include <limits>


namespace Statistics {
	template<typename T>
	class ProbabilityDensityFunction : public XYSeries<T> {
		protected:
			std::uniform_real_distribution<T> distribution;
			// RNG is per worker thread (thread_local in draw_sample) — the PDF is shared (shared_ptr)
			// across MT workers via the radiation source, so a shared generator would be a data race.
			XYSeries<T> cdf;

		public:
			ProbabilityDensityFunction(const std::vector<std::pair<T, T>>& points)
				: XYSeries<T>(),
				  distribution(std::numeric_limits<T>::min(), 1.0 + std::numeric_limits<T>::min()),
				  cdf()
			{
				if (points.size() > 0) {
					// norm all points in pdf and cdf
					T sum = 0.0;
					for (auto& point : points) {
						sum += point.second;
					}

					if (sum == 0.0)
						throw std::runtime_error("ProbabilityDensityFunction: sum of all points is zero");

					T integral = 0.0;
					for (auto& point : points)
					{
						T y = point.second / sum;
						this->add(point.first, integral);
						this->cdf.add(integral, point.first);
						integral += y;
					}
					this->finalize();
					this->cdf.finalize();
				}
			};

			T draw_sample() {
				// Per-worker-thread RNG: this PDF is shared across MT workers, so a shared generator would
				// be a data race. thread_local gives each worker its own source (same std::random_device
				// behaviour as before, just no longer shared). The distribution member is read-only here.
				thread_local std::random_device generator;
				T random_probability = this->distribution(generator);
				return this->cdf.get_interpolated(random_probability, &Statistics::Interpolation::Linear<T>);
			};

			T min() const {
				return this->cdf.minY();
			};

			T max() const {
				return this->cdf.maxY();
			};
	};
}
