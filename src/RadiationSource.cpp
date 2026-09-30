#include "RadiationSource.hpp"
#include <glm/gtc/constants.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>


using namespace RadiationSimulation;

RadiationSource::RadiationSource(float energy_eV, std::string particle_name, std::unique_ptr<ISourceShape> shape)
	: particle_name(particle_name),
	  energy_eV(energy_eV),
	  shape(std::move(shape))
{
}

void RadiationSource::setTransform(const glm::vec3& location, const glm::vec3& orientation)
{
	this->location = location;
	this->rotation = glm::quat_cast(glm::mat4(1.0f));
	const glm::vec3 normalizedOrientation = glm::normalize(orientation);
	const glm::vec3 up(0.0f, 0.0f, -1.0f);

	if (normalizedOrientation == up) {
		return;
	}
	glm::vec3 axis = glm::cross(up, normalizedOrientation);
	// exactly opposite to the reference axis the shortest rotation has no axis of its own: turn about X
	if (glm::length(axis) < 1e-6f) {
		this->rotation = glm::angleAxis(glm::pi<float>(), glm::vec3(1.0f, 0.0f, 0.0f));
		return;
	}
	float angle = acos(glm::clamp(glm::dot(up, normalizedOrientation), -1.0f, 1.0f));
	this->rotation = glm::angleAxis(angle, glm::normalize(axis));
}

glm::quat RadiationSource::getCArmRotation() const
{
	// base pose: beam along +Y, which the source reaches from its reference axis -Z by a quarter turn about X
	const glm::quat base_pose = glm::angleAxis(glm::half_pi<float>(), glm::vec3(1.0f, 0.0f, 0.0f));
	return this->rotation * glm::inverse(base_pose);
}

glm::vec3 RadiationSimulation::RadiationSource::drawRayDirection(const UniformRandom& uniform)
{
	return glm::vec3(this->rotation * glm::vec4(this->shape->drawRayDirection(uniform), 0.f));
}

XRaySource::XRaySource(float energy_eV, std::unique_ptr<ISourceShape> shape)
	: RadiationSource(energy_eV, "gamma", std::move(shape))
{
}

RadiationSimulation::XRaySpectrumSource::XRaySpectrumSource(std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> spectrum_probabilities, std::unique_ptr<ISourceShape> shape, float energy_lower_cut_eV, float max_energy_eV)
	: XRaySource(0.0f, std::move(shape)),
	  energy_lower_cut_eV(energy_lower_cut_eV),
	  spectrum_probabilities(spectrum_probabilities)
{
	if (!(energy_lower_cut_eV >= 0.f && energy_lower_cut_eV <= max_energy_lower_cut_eV))
		throw std::invalid_argument("The lower energy cut of " + std::to_string(energy_lower_cut_eV) + " eV lies outside [0, " + std::to_string(max_energy_lower_cut_eV) + "] eV; photons above 5 keV must never be cut.");

	if (max_energy_eV <= 0.f)
		max_energy_eV = spectrum_probabilities->max();

	// relative tolerance: bin edges are computed in double from float energies
	if (max_energy_eV > 0.f && max_energy_eV * (1.0 + 1e-6) < spectrum_probabilities->max()) {
		std::cerr << "Max energy value is lower than the maximum value in the spectrum: " << max_energy_eV << " < " << spectrum_probabilities->max() << std::endl;
		throw std::runtime_error("Max energy value is lower than the maximum value in the spectrum");
	}

	this->lower_cut_probability = spectrum_probabilities->cdf(energy_lower_cut_eV);
	if (this->lower_cut_probability >= 1.0)
		throw std::runtime_error("The spectrum has no energies above the lower cut of " + std::to_string(energy_lower_cut_eV) + " eV");

	this->generated_counts = std::vector<std::atomic<uint64_t>>(std::max<size_t>(static_cast<size_t>(max_energy_eV / generated_bin_width_eV), 1));
}

RadiationSimulation::XRaySpectrumSource::~XRaySpectrumSource()
{
}

std::vector<uint64_t> RadiationSimulation::XRaySpectrumSource::getGeneratedCounts() const
{
	std::vector<uint64_t> counts(this->generated_counts.size());
	for (size_t i = 0; i < counts.size(); i++)
		counts[i] = this->generated_counts[i].load(std::memory_order_relaxed);
	return counts;
}

size_t RadiationSimulation::XRaySpectrumSource::getPossibilitiesCount() const
{
	return this->spectrum_probabilities->get_point_count();
}

float RadiationSimulation::XRaySpectrumSource::drawEnergy_eV(const UniformRandom& uniform)
{
	// inverse transform restricted to the part of the spectrum above the lower cut: exact, no rejection
	const double p = this->lower_cut_probability + uniform() * (1.0 - this->lower_cut_probability);
	const float energy_eV = this->spectrum_probabilities->quantile(p);

	const size_t bin = std::min(static_cast<size_t>((energy_eV + 0.5f * generated_bin_width_eV) / generated_bin_width_eV), this->generated_counts.size() - 1);
	this->generated_counts[bin].fetch_add(1, std::memory_order_relaxed);

	return energy_eV;
}

std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> RadiationSimulation::SpectrumLoader::LoadSpectrum(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file.is_open()) {
		throw std::runtime_error("Failed to open file: " + filename);
	}

	std::string line;
	float energy_unit = 0.f;
	bool second_column_is_fluence = false;
	// Read first line as header
	if (!std::getline(file, line)) {
		throw std::runtime_error("File is empty");
	}
	else {
		while (line.starts_with("#") || line.starts_with("//")) {
			if (line.starts_with("# energy / keV")) {
				energy_unit = 1e+3;
				second_column_is_fluence = true;
			}
			std::getline(file, line);
		}
		std::replace(line.begin(), line.end(), ',', ' ');
		// check if it is indeed a spectrum file
		std::istringstream iss(line);
		std::string column_name;
		size_t column_count = 0;
		while (iss >> column_name) {
			column_count++;

			size_t start = column_name.find('[');
			size_t end = column_name.find(']');
			if (start != std::string::npos && end != std::string::npos && end > start) {
				std::string plain_column_name = column_name.substr(0, start);
				if (plain_column_name.compare("Energy") == 0 && column_count == 1) {
					std::string unit = column_name.substr(start + 1, end - start - 1);
					if (unit.compare("MeV") == 0) {
						energy_unit = 1e+6;
					}
					else if (unit.compare("keV") == 0) {
						energy_unit = 1e+3;
					}
					else if (unit.compare("eV") == 0) {
						energy_unit = 1.f;
					}
					else {
						throw std::runtime_error("Unknown energy unit: " + unit);
					}
				}
				if (plain_column_name.compare("Fluence") == 0 && column_count == 2) {
					second_column_is_fluence = true;
				}
			}
		}
	}

	if (!second_column_is_fluence)
		throw std::runtime_error("Second column is not Fluence");

	if (energy_unit <= 0.f)
		throw std::runtime_error("Energy unit is not set");

	std::map<float, float> energy_fluence_map;
	// Read the rest of the lines
	while (std::getline(file, line)) {
		if (line.starts_with("#") || line.starts_with("//")) {
			std::getline(file, line);
		}
		std::replace(line.begin(), line.end(), ',', ' ');
		std::istringstream iss(line);
		double energy, fluence;
		if (!(iss >> energy >> fluence)) {
			std::cerr << "Failed to parse line: " << line << std::endl;
			continue;
		}
		float f_energy = static_cast<float>(energy * energy_unit);
		if (energy_fluence_map.find(f_energy) != energy_fluence_map.end()) {
			energy_fluence_map[f_energy] += static_cast<float>(fluence);
		}
		else {
			energy_fluence_map[f_energy] = static_cast<float>(fluence);
		}
	}

	std::vector<std::pair<float, float>> energy_fluence_points;

	for (auto& itr : energy_fluence_map) {
		energy_fluence_points.push_back(itr);
	}

	return std::make_shared<Statistics::ProbabilityDensityFunction<float>>(energy_fluence_points);

}

RadiationSimulation::ConeSourceShape::ConeSourceShape(float opening_angle_deg)
	: opening_angle_radians(glm::radians(opening_angle_deg))
{
}

glm::vec2 RadiationSimulation::ConeSourceShape::getHalfTangents() const
{
	// the opening angle is the cone's half angle (see drawRayDirection)
	const float t = (this->opening_angle_radians < glm::half_pi<float>()) ? std::tan(this->opening_angle_radians) : INFINITY;
	return glm::vec2(t);
}

glm::vec3 RadiationSimulation::ConeSourceShape::drawRayDirection(const UniformRandom& uniform)
{
	// uniform per solid angle within the cone's half opening angle
	const float azimuth = 2.0f * glm::pi<float>() * static_cast<float>(uniform());
	const float polar = std::acos(1.0f - static_cast<float>(uniform()) * (1.0f - std::cos(this->opening_angle_radians)));

	return glm::vec3(
		std::sin(polar) * std::cos(azimuth),
		std::sin(polar) * std::sin(azimuth),
		-std::cos(polar)
	);
}

RadiationSimulation::RectangleSourceShape::RectangleSourceShape(const glm::vec2& size, float distance)
	: size(size),
	  distance(distance)
{
	if (!(distance > 0.f))
		throw std::invalid_argument("The distance of a rectangle source's field size must be positive.");
}

glm::vec3 RadiationSimulation::RectangleSourceShape::drawRayDirection(const UniformRandom& uniform)
{
	// A point source behind a rectangular collimator: uniform per solid angle inside the pyramid through the rectangle.
	// Such a source puts cos³θ photons per area onto the flat rectangle, so a uniform point on it is kept with that
	// probability.
	const double distance = static_cast<double>(this->distance);
	while (true) {
		const double x = (uniform() - 0.5) * static_cast<double>(this->size.x);
		const double y = (uniform() - 0.5) * static_cast<double>(this->size.y);
		const double length = std::sqrt(x * x + y * y + distance * distance);
		const double cos_theta = distance / length;
		if (uniform() < cos_theta * cos_theta * cos_theta)
			return glm::vec3(static_cast<float>(x / length), static_cast<float>(y / length), static_cast<float>(-cos_theta));
	}
}

RadiationSimulation::EllipsoidSourceShape::EllipsoidSourceShape(const glm::vec2& half_angles)
	: half_angles_degrees(half_angles)
{
	if (!(half_angles.x > 0.f && half_angles.x < 90.f && half_angles.y > 0.f && half_angles.y < 90.f))
		throw std::invalid_argument("The half opening angles of an ellipsoid source must lie in (0, 90) degrees.");
	this->tan_half_angles = glm::dvec2(std::tan(glm::radians(static_cast<double>(half_angles.x))), std::tan(glm::radians(static_cast<double>(half_angles.y))));
	this->cos_enclosing_angle = std::cos(glm::radians(static_cast<double>(std::max(half_angles.x, half_angles.y))));
}

glm::vec3 RadiationSimulation::EllipsoidSourceShape::drawRayDirection(const UniformRandom& uniform)
{
	// Uniform per solid angle in the circular cone around the larger half angle, kept if inside the elliptical cone:
	// rejection keeps the distribution uniform per solid angle.
	while (true) {
		const double azimuth = 2.0 * glm::pi<double>() * uniform();
		const double cos_polar = 1.0 - uniform() * (1.0 - this->cos_enclosing_angle);
		const double sin_polar = std::sqrt(std::max(0.0, 1.0 - cos_polar * cos_polar));
		// slopes of the direction against the beam axis, i.e. its point on the plane at unit distance
		const double slope_x = sin_polar * std::cos(azimuth) / cos_polar;
		const double slope_y = sin_polar * std::sin(azimuth) / cos_polar;
		const double ex = slope_x / this->tan_half_angles.x;
		const double ey = slope_y / this->tan_half_angles.y;
		if (ex * ex + ey * ey <= 1.0)
			return glm::vec3(static_cast<float>(sin_polar * std::cos(azimuth)), static_cast<float>(sin_polar * std::sin(azimuth)), static_cast<float>(-cos_polar));
	}
}
