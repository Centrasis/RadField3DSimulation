#pragma once
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>


namespace Statistics {
	/** Probability distribution of a tabulated spectrum, sampled by inverting its cumulative distribution.
	* Every point (x, weight) is the centre of a bin that reaches halfway to the neighbouring points; the first and last
	* bins reach as far outwards as inwards (the first one never below 0). Values are uniformly distributed within their
	* bin. A single point is a delta distribution at its x.
	* This is the convention of the spectra RadField3D reads: SpekPy lists mid-bin energies, and the tube spectrum stored
	* in a field labels each bin with its centre.
	* The distribution holds no random number engine: the caller passes a uniform random number to sample().
	*/
	template<typename T>
	class ProbabilityDensityFunction {
	protected:
		std::vector<double> lower_edges;
		std::vector<double> upper_edges;
		// cumulative[i] is the probability of all bins before bin i; cumulative.back() == 1
		std::vector<double> cumulative;
		size_t first_bin = 0;
		size_t last_bin = 0;

	public:
		ProbabilityDensityFunction(std::vector<std::pair<T, T>> points)
		{
			if (points.empty())
				throw std::runtime_error("ProbabilityDensityFunction: no points given");
			std::sort(points.begin(), points.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

			const size_t n = points.size();
			double total = 0.0;
			for (const auto& point : points) {
				if (point.second < 0)
					throw std::runtime_error("ProbabilityDensityFunction: negative weight");
				total += static_cast<double>(point.second);
			}
			if (total <= 0.0)
				throw std::runtime_error("ProbabilityDensityFunction: sum of all points is zero");

			this->lower_edges.resize(n);
			this->upper_edges.resize(n);
			for (size_t i = 0; i < n; i++) {
				const double x = static_cast<double>(points[i].first);
				const double half_below = (i > 0) ? 0.5 * (x - static_cast<double>(points[i - 1].first)) : (n > 1 ? 0.5 * (static_cast<double>(points[1].first) - x) : 0.0);
				const double half_above = (i + 1 < n) ? 0.5 * (static_cast<double>(points[i + 1].first) - x) : half_below;
				this->lower_edges[i] = std::max(0.0, x - half_below);
				this->upper_edges[i] = x + half_above;
			}

			this->first_bin = n;
			for (size_t i = 0; i < n; i++) {
				if (points[i].second > 0) {
					this->first_bin = std::min(this->first_bin, i);
					this->last_bin = i;
				}
			}

			this->cumulative.resize(n + 1, 0.0);
			for (size_t i = 0; i < n; i++)
				this->cumulative[i + 1] = this->cumulative[i] + static_cast<double>(points[i].second) / total;
			// no rounding residue after the last bin with a weight
			for (size_t i = this->last_bin + 1; i <= n; i++)
				this->cumulative[i] = 1.0;
		}

		/** Probability of a value below x. */
		double cdf(double x) const
		{
			if (x <= this->lower_edges[this->first_bin])
				return 0.0;
			if (x >= this->upper_edges[this->last_bin])
				return 1.0;
			const size_t bin = static_cast<size_t>(std::upper_bound(this->lower_edges.begin(), this->lower_edges.end(), x) - this->lower_edges.begin()) - 1;
			const double width = this->upper_edges[bin] - this->lower_edges[bin];
			const double inside = (width > 0.0) ? std::clamp((x - this->lower_edges[bin]) / width, 0.0, 1.0) : 1.0;
			return this->cumulative[bin] + inside * (this->cumulative[bin + 1] - this->cumulative[bin]);
		}

		/** Probability of a value in [low, high). */
		double probability(double low, double high) const
		{
			return (high > low) ? this->cdf(high) - this->cdf(low) : 0.0;
		}

		/** The value below which a share p of the distribution lies (inverse of cdf), p in [0, 1]. */
		T quantile(double p) const
		{
			p = std::clamp(p, 0.0, 1.0);
			size_t bin = static_cast<size_t>(std::upper_bound(this->cumulative.begin(), this->cumulative.end(), p) - this->cumulative.begin());
			bin = std::clamp(bin, this->first_bin + 1, this->last_bin + 1) - 1;
			const double probability = this->cumulative[bin + 1] - this->cumulative[bin];
			const double inside = std::clamp((p - this->cumulative[bin]) / probability, 0.0, 1.0);
			return static_cast<T>(this->lower_edges[bin] + inside * (this->upper_edges[bin] - this->lower_edges[bin]));
		}

		/** Draws a value from the distribution.
		* @param uniform A random number uniformly distributed in [0, 1).
		*/
		T sample(double uniform) const
		{
			return this->quantile(uniform);
		}

		/** Smallest value the distribution can produce: the lower edge of the first bin with a weight. */
		T min() const
		{
			return static_cast<T>(this->lower_edges[this->first_bin]);
		}

		/** Largest value the distribution can produce: the upper edge of the last bin with a weight. */
		T max() const
		{
			return static_cast<T>(this->upper_edges[this->last_bin]);
		}

		size_t get_point_count() const
		{
			return this->lower_edges.size();
		}
	};
}
