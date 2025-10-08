#pragma once
#include <string>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/quaternion.hpp>
#include <memory>
#include <random>
#include "ProbabilityFunctions.hpp"
#include "RadFiled3D/Voxel.hpp"
#include <mutex>
#include <shared_mutex>

namespace RadiationSimulation {

	/**
	 * @brief Interface for source shapes.
	 */
	class ISourceShape {
	public:
		/**
		 * @brief Draws a ray direction at random.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection() = 0;
	};

	/**
	 * @brief Cone-shaped radiation source.
	 */
	class ConeSourceShape : public ISourceShape {
	protected:
		std::mt19937 rand_engine; ///< Random engine for generating directions.
		std::uniform_real_distribution<float> rot_angle_distribution; ///< Distribution for rotation angles.
		float opening_angle_radians; ///< Opening angle in radians.
	public:
		/**
		 * @brief Constructor for ConeSourceShape.
		 * @param opening_angle_deg Opening angle in degrees.
		 */
		ConeSourceShape(float opening_angle_deg);

		/**
		 * @brief Draws a ray direction within the cone.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection() override;

		float getOpeningAngleDegrees() const { return glm::degrees(this->opening_angle_radians); }
	};

	/**
	 * @brief Rectangle-shaped radiation source.
	 */
	class RectangleSourceShape : public ISourceShape {
	protected:
		std::mt19937 rand_engine; ///< Random engine for generating directions.
		std::uniform_real_distribution<float> distribution_x; ///< Distribution for generating directions.
		std::uniform_real_distribution<float> distribution_y; ///< Distribution for generating directions.
		glm::vec2 size; ///< Size of the rectangle.
		float distance; ///< Distance from the source.
	public:
		/**
		 * @brief Constructor for RectangleSourceShape.
		 * @param size Size of the rectangle at the specified distance.
		 * @param distance Distance from the source at which the rectangle size shouldbe fulfilled.
		 */
		RectangleSourceShape(const glm::vec2& size, float distance);

		/**
		 * @brief Draws a ray direction within the rectangle.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection() override;

		glm::vec2 getFieldSizeMeters() const { return this->size; }
	};

	/**
	 * @brief Ellipsoid-shaped radiation source.
	 */
	class EllipsoidSourceShape : public ISourceShape {
	protected:
		glm::vec2 angles; ///< Angles defining the ellipsoid.
		std::mt19937 rand_engine;  ///< Random engine for generating directions.
		std::uniform_real_distribution<float> distribution; ///< Distribution for generating directions.
	public:
		/**
		 * @brief Constructor for EllipsoidSourceShape.
		 * @param angles Angles defining the ellipsoid.
		 */
		EllipsoidSourceShape(const glm::vec2& angles);

		/**
		 * @brief Draws a ray direction within the ellipsoid.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection() override;

		glm::vec2 getOpeningAnglesDegrees() const { return this->angles; }
	};

	/**
	 * @brief Base class for radiation sources.
	 */
	class RadiationSource {
	protected:
		float energy_eV; ///< Energy of the radiation in electron volts.
		const std::string particle_name; ///< Name of the particle.
		glm::quat rotation = glm::quat_cast(glm::mat3(1.f)); ///< Transformation matrix.
		glm::vec3 location = glm::vec3(0.f);
		std::unique_ptr<ISourceShape> shape; ///< Shape of the radiation source.
	public:
		/**
		 * @brief Constructor for RadiationSource.
		 * @param energy_eV Energy of the radiation in electron volts.
		 * @param particle_name Name of the particle.
		 * @param shape Shape of the radiation source.
		 */
		RadiationSource(float energy_eV, std::string particle_name, std::unique_ptr<ISourceShape> shape);

		/**
		 * @brief Gets the name of the particle.
		 * @return The name of the particle.
		 */
		inline const std::string& getParticleName() const { return this->particle_name; }

		/**
		 * @brief Draws the energy of the radiation.
		 * @return The energy of the radiation in electron volts.
		 */
		virtual float drawEnergy_eV() { return this->energy_eV; }

		/**
		 * @brief Gets the number of possibilities.
		 * @return The number of possibilities.
		 */
		virtual size_t getPossibilitiesCount() const { return 1; }

		/**
		 * @brief Gets the location of the radiation source.
		 * @return The location as a glm::vec3.
		 */
		inline glm::vec3 getLocation() const { return this->location; }

		/**
		 * @brief Gets the rotation quaternion of the radiation source.
		 * @return The transformation matrix.
		 */
		inline const glm::quat& getRotation() const { return this->rotation; }

		/**
		 * @brief Draws a ray direction.
		 * @return The direction of the ray as a glm::vec3.
		 */
		glm::vec3 drawRayDirection();

		/**
		 * @brief Sets the transformation matrix of the radiation source.
		 * @param location The location of the source.
		 * @param orientation The orientation of the source.
		 */
		void setTransform(const glm::vec3& location, const glm::vec3& orientation);

		const ISourceShape* getShape() const { return this->shape.get(); }
	};

	/**
	 * @brief X-ray radiation source.
	 */
	class XRaySource : public RadiationSource {
	public:
		/**
		 * @brief Constructor for XRaySource.
		 * @param energy_eV Energy of the radiation in electron volts.
		 * @param shape Shape of the radiation source.
		 */
		XRaySource(float energy_eV, std::unique_ptr<ISourceShape> shape);
	};

	/**
	 * @brief X-ray spectrum radiation source.
	 * @details This source generates radiation with a given spectrum. Thread-safe.
	 */
	class XRaySpectrumSource : public XRaySource {
	protected:
		std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> spectrum_probabilities; ///< Spectrum probabilities.
		const float energy_lower_cut_eV; ///< Lower cut-off energy in electron volts.
		RadFiled3D::HistogramVoxel hist; ///< Histogram voxel for the generated spectrum.
		std::vector<float> hist_buffer; ///< Buffer for the histogram.
		mutable std::shared_mutex hist_mutex; ///< Mutex for thread-safe access to the histogram.
	public:
		/**
		 * @brief Constructor for XRaySpectrumSource.
		 * @param spectrum_probabilities Spectrum probabilities.
		 * @param shape Shape of the radiation source.
		 * @param energy_lower_cut_eV Lower cut-off energy in electron volts.
		 */
		XRaySpectrumSource(std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> spectrum_probabilities, std::unique_ptr<ISourceShape> shape, float energy_lower_cut_eV = 3e+4, float max_energy_eV = 0.f);

		/**
		 * @brief Destructor for XRaySpectrumSource.
		 */
		~XRaySpectrumSource();

		/**
		 * @brief Gets the number of possibilities.
		 * @return The number of possibilities.
		 */
		virtual size_t getPossibilitiesCount() const override;

		/**
		 * @brief Draws the energy of the radiation.
		 * @return The energy of the radiation in electron volts.
		 */
		virtual float drawEnergy_eV() override;

		/**
		 * @brief Gets the generated spectrum.
		 * @return The generated spectrum as a HistogramVoxel.
		 */
		const RadFiled3D::HistogramVoxel& getGeneratedSpectrum() const { return this->hist; }
	};

	/**
	 * @brief Loader for radiation spectra from csv files.
	 */
	class SpectrumLoader {
	public:
		/**
		 * @brief Loads a spectrum from a csv file.
		 * @param filename The name of the file.
		 * @return A shared pointer to the loaded spectrum.
		 */
		static std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> LoadSpectrum(const std::string& filename);
	};
}