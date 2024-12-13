#pragma once
#include <memory>
#include "RadFiled3D/RadiationField.hpp"
#include <utils/PCA.hpp>
#include <glm/vec3.hpp>
#include <RadFiled3D/Voxel.hpp>


namespace RadiationSimulation {
	class Voxel {
	public:
		virtual inline const RadFiled3D::IVoxel* getVoxel() const = 0;
		virtual void normalize(size_t particle_count) = 0;
	};

	class FloatVoxel : public Voxel {
	protected:
		RadFiled3D::ScalarVoxel<float>& voxel;
	public:
		FloatVoxel(RadFiled3D::ScalarVoxel<float>& voxel) : voxel(voxel) {}
		inline const RadFiled3D::IVoxel* getVoxel() const override { return &this->voxel; }
		inline void add_value(float value) { this->voxel += value; }
		virtual void normalize(size_t particle_count) override {
			this->voxel /= static_cast<float>(particle_count);
			assert(this->value() == this->value() && "NAN in voxel");
		}

		float& value() const { return this->voxel.get_data(); }
	};

	class CounterVoxel : public Voxel {
	protected:
		RadFiled3D::ScalarVoxel<float>& voxel;
	public:
		CounterVoxel(RadFiled3D::ScalarVoxel<float>& voxel) : voxel(voxel) {}
		inline const RadFiled3D::IVoxel* getVoxel() const override { return &this->voxel; }
		inline void add_value() { this->voxel += 1.f; }
		virtual void normalize(size_t particle_count) override {
			this->voxel /= static_cast<float>(particle_count);
		}

		float& value() const { return this->voxel.get_data(); }
	};

	class DirectionVoxel : public Voxel {
	protected:
		RadFiled3D::ScalarVoxel<glm::vec3>& voxel;
		OnlinePCA pca;
	public:
		DirectionVoxel(RadFiled3D::ScalarVoxel<glm::vec3>& voxel) : voxel(voxel) {}
		inline const RadFiled3D::IVoxel* getVoxel() const override { return &this->voxel; }
		inline void add_value(const glm::vec3& dir) { pca.addVector(dir); }
		virtual void normalize(size_t particle_count) override {
			this->voxel = pca.getPrincipalDirection();
		}

		glm::vec3& value() const { return this->voxel.get_data(); }
	};

	class SpectrumVoxel : public Voxel {
	protected:
		RadFiled3D::HistogramVoxel& voxel;
	public:
		SpectrumVoxel(RadFiled3D::HistogramVoxel& voxel) : voxel(voxel) {}

		inline const RadFiled3D::IVoxel* getVoxel() const override { return &this->voxel; }

		inline void add_value(double value) {
			size_t index = static_cast<size_t>(value / this->voxel.get_histogram_bin_width());
			if (index < this->voxel.get_bins()) {
				(&this->voxel.get_data())[index] += 1.0;
			}
			else {
				(&this->voxel.get_data())[this->voxel.get_bins() - 1] += 1.0;
			}
		}

		virtual void normalize(size_t particle_count) override;

		std::span<float> value() const { return this->voxel.get_histogram(); }
	};

	class RadiationDetector {
	protected:
		FloatVoxel energy;
		CounterVoxel hits;
		DirectionVoxel direction;
		SpectrumVoxel spectrum;
		float energy_squared_sum = 0.f;
		
	public:
		RadiationDetector(size_t spectrum_bins, double spectrum_bin_width, std::shared_ptr<RadFiled3D::VoxelGridBuffer> voxels, const glm::uvec3& coord);
		void normalize(size_t particle_count);
	};

	class RadiationFieldDetector {
		friend class RadiationFieldDetector;
	protected:
		std::shared_ptr<RadFiled3D::CartesianRadiationField> field;
		const size_t spectra_bins;
		const double spectra_bin_width;
		virtual void evaluate_field() = 0;
	public:
		RadiationFieldDetector(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, size_t spectra_bins, double spectra_bin_width);
		virtual void SetUp() = 0;
		virtual std::shared_ptr<RadFiled3D::IRadiationField> get_normalized_field_copy() = 0;
		std::shared_ptr<RadFiled3D::IRadiationField> evaluate();
		virtual float get_statistical_error(size_t primary_particle_count = 0) = 0;
		virtual size_t get_number_of_tracked_particles() = 0;
	};
};