#pragma once
#include <string>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <glm/gtc/quaternion.hpp>
#include <functional>
#include <memory>
#include "ProbabilityFunctions.hpp"
#include "radfiled3d/voxel.hpp"
#include <atomic>
#include <cstdint>
#include <vector>

namespace RadiationSimulation {
	/** Source of random numbers uniformly distributed in [0, 1). The Geant4 primary generator passes the worker's
	* per-event engine (G4UniformRand), so every sampled primary is determined by the run's seed.
	*/
	using UniformRandom = std::function<double()>;

	/**
	 * @brief Interface for source shapes.
	 */
	class ISourceShape {
	public:
		virtual ~ISourceShape() = default;
		/**
		 * @brief Draws a ray direction at random.
		 * @param uniform Source of the random numbers.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection(const UniformRandom& uniform) = 0;

		/**
		 * @brief Half size of the beam at unit distance from the focal spot, along the source frame's x and y axes
		 * (the tangents of the half opening angles). Infinite for beams opening 90 degrees or more.
		 */
		virtual glm::vec2 getHalfTangents() const = 0;
	};

	/**
	 * @brief Cone-shaped radiation source.
	 */
	class ConeSourceShape : public ISourceShape {
	protected:
		float opening_angle_radians; ///< Opening angle in radians.
	public:
		/**
		 * @brief Constructor for ConeSourceShape.
		 * @param opening_angle_deg Opening angle in degrees.
		 */
		ConeSourceShape(float opening_angle_deg);

		/**
		 * @brief Draws a ray direction within the cone.
		 * @param uniform Source of the random numbers.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection(const UniformRandom& uniform) override;

		virtual glm::vec2 getHalfTangents() const override;

		float getOpeningAngleDegrees() const { return glm::degrees(this->opening_angle_radians); }
	};

	/**
	 * @brief Rectangle-shaped radiation source: a point source behind a rectangular collimator.
	 * Rays are uniform per solid angle inside the pyramid through the rectangle, so the photons per area of the
	 * rectangle's plane fall off towards its edges as cos³θ.
	 */
	class RectangleSourceShape : public ISourceShape {
	protected:
		glm::vec2 size; ///< Size of the rectangle.
		float distance; ///< Distance from the source.
	public:
		/**
		 * @brief Constructor for RectangleSourceShape.
		 * @param size Size of the rectangle at the specified distance.
		 * @param distance Distance from the source at which the rectangle size shouldbe fulfilled, must be positive.
		 */
		RectangleSourceShape(const glm::vec2& size, float distance);

		/**
		 * @brief Draws a ray direction within the rectangle.
		 * @param uniform Source of the random numbers.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection(const UniformRandom& uniform) override;

		virtual glm::vec2 getHalfTangents() const override { return this->size / (2.f * this->distance); }

		glm::vec2 getFieldSizeMeters() const { return this->size; }
	};

	/**
	 * @brief Ellipsoid-shaped radiation source.
	 */
	class EllipsoidSourceShape : public ISourceShape {
	protected:
		glm::vec2 half_angles_degrees; ///< Half opening angles along x and y in degrees.
		glm::dvec2 tan_half_angles; ///< Semi-axes of the ellipse the beam cuts from the plane at unit distance.
		double cos_enclosing_angle; ///< Cosine of the larger half opening angle.
	public:
		/**
		 * @brief Constructor for EllipsoidSourceShape: an elliptical cone around the beam axis, as behind an elliptical
		 * collimator aperture. At unit distance along the axis, the beam covers the ellipse with the semi-axes
		 * tan(half_angles.x) and tan(half_angles.y).
		 * @param half_angles Half opening angles along x and y in degrees, each in (0, 90).
		 * @throws std::invalid_argument if an angle lies outside (0, 90) degrees.
		 */
		EllipsoidSourceShape(const glm::vec2& half_angles);

		/**
		 * @brief Draws a ray direction uniformly per solid angle within the elliptical cone.
		 * @param uniform Source of the random numbers.
		 * @return The direction of the ray as a glm::vec3.
		 */
		virtual glm::vec3 drawRayDirection(const UniformRandom& uniform) override;

		virtual glm::vec2 getHalfTangents() const override { return glm::vec2(this->tan_half_angles); }

		glm::vec2 getOpeningAnglesDegrees() const { return this->half_angles_degrees; }
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
		 * @param uniform Source of the random numbers.
		 * @return The energy of the radiation in electron volts.
		 */
		virtual float drawEnergy_eV(const UniformRandom& uniform) { return this->energy_eV; }

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
		 * @brief Rotation of the tube (and with it the C-arm) from its base pose: below the isocentre, beam along +Y
		 * (the pose of --source-theta 90 --source-phi 0). Uses the rotation that orients the beam and its collimated
		 * field, so meshes turned with it stay aligned with the field.
		 */
		glm::quat getCArmRotation() const;

		/**
		 * @brief Draws a ray direction.
		 * @param uniform Source of the random numbers.
		 * @return The direction of the ray as a glm::vec3.
		 */
		glm::vec3 drawRayDirection(const UniformRandom& uniform);

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
		double lower_cut_probability; ///< Share of the spectrum below the lower cut, which is never sampled.
		/// Generated energies per 1 keV bin, bin i centred at i keV and the last bin also holding everything above.
		/// 64-bit atomic counts: exact beyond 2^24 per bin, and no lock between the workers.
		std::vector<std::atomic<uint64_t>> generated_counts;
		static constexpr float generated_bin_width_eV = 1e+3f;
	public:
		static constexpr float max_energy_lower_cut_eV = 5e+3f; ///< Highest allowed lower cut: photons above 5 keV are never cut.

		/**
		 * @brief Constructor for XRaySpectrumSource.
		 * @param spectrum_probabilities Spectrum probabilities.
		 * @param shape Shape of the radiation source.
		 * @param energy_lower_cut_eV Lower cut-off energy in electron volts; the spectrum below it is not sampled. At most
		 *        max_energy_lower_cut_eV, so no photon above 5 keV is ever cut.
		 * @param max_energy_eV Upper end of the generated-spectrum histogram; must reach the spectrum's maximum.
		 */
		XRaySpectrumSource(std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> spectrum_probabilities, std::unique_ptr<ISourceShape> shape, float energy_lower_cut_eV = 0.f, float max_energy_eV = 0.f);

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
		 * @brief Draws the energy of the radiation from the spectrum above the lower cut.
		 * @param uniform Source of the random numbers.
		 * @return The energy of the radiation in electron volts.
		 */
		virtual float drawEnergy_eV(const UniformRandom& uniform) override;

		/** The loaded spectrum. */
		const Statistics::ProbabilityDensityFunction<float>& getSpectrum() const { return *this->spectrum_probabilities; }

		/** Number of bins of the generated spectrum. */
		size_t getGeneratedSpectrumBins() const { return this->generated_counts.size(); }

		/** Width of the generated spectrum's bins in eV. */
		float getGeneratedSpectrumBinWidth_eV() const { return generated_bin_width_eV; }

		/** Snapshot of the generated energies per bin; safe while workers keep drawing. */
		std::vector<uint64_t> getGeneratedCounts() const;
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