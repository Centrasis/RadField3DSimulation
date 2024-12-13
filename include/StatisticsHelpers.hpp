#pragma once
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <iostream>


namespace Statistics {
	namespace Interpolation {
		template<typename T>
		static T Linear(const T& y0, const T& y1, const double& alpha) {
			return y0 + alpha * (y1 - y0);
		};
	}

	template<typename T>
	class XYSeries {
		using XYMap = std::map<T, T>;
		using Interpolation = std::function<T(T, T, double)>;
		using iterator = typename XYMap::iterator;
		using const_iterator = typename XYMap::const_iterator;
		private:
			XYMap xy_map;
			std::vector<T> x_values;
		public:
			XYSeries() {};

			void finalize() {
				this->x_values.clear();
				for (auto& xy : this->xy_map) {
					this->x_values.push_back(xy.first);
				}
				std::sort(this->x_values.begin(), this->x_values.end(), [](T x, T y) { return x <= y; });
			};

			virtual void add(const T& x, const T& y) {
				this->xy_map.insert({ x, y });
			};

			inline void add(const std::pair<T, T>& xy) { this->add(xy.first, xy.second); };

			virtual T get(const T& x) const {
				return this->xy_map.find(x)->second;
			};

			virtual double get_interpolated(const T& x, const Interpolation& interpolation) {
				if (this->x_values.size() != this->xy_map.size()) {
					this->finalize();
				}
				//const std::vector<double>::iterator lower = std::lower_bound<std::vector<double>::iterator, double>(this->x_values.begin(), this->x_values.end(), x);
				auto upper = std::upper_bound<typename std::vector<T>::iterator, T>(this->x_values.begin(), this->x_values.end(), x);
				auto lower = (upper != this->x_values.begin()) ? std::prev(upper) : upper;

				const T& lower_y = this->xy_map.find(*lower)->second;
				if (upper == this->x_values.end() || lower == upper) {
					return lower_y;
				}

				const T& upper_y = this->xy_map.find(*upper)->second;
				
				T distance_x = *upper - *lower;
				T alpha = (distance_x > 0) ? (x - *lower) / (distance_x) : 0.0;
				return interpolation(lower_y, upper_y, alpha);
			};

			virtual size_t get_point_count() const { return this->xy_map.size(); };

			virtual iterator begin() noexcept { return this->xy_map.begin(); };
			virtual iterator end() noexcept { return this->xy_map.end(); };
			virtual const_iterator begin() const noexcept { return this->xy_map.begin(); };
			virtual const_iterator end() const noexcept { return this->xy_map.end(); };

			virtual T minX() const {
				return this->xy_map.begin()->first;
			};

			virtual T maxX() const {
				return this->xy_map.rbegin()->first;
			};

			virtual T minY() const {
				return this->xy_map.begin()->second;
			};

			virtual T maxY() const {
				return this->xy_map.rbegin()->second;
			};
	};
}
