#include "RadiationSource.hpp"
#include <fstream>
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
	float angle = acos(glm::clamp(glm::dot(up, normalizedOrientation), -1.0f, 1.0f));
	this->rotation = glm::angleAxis(angle, glm::normalize(axis));
}

glm::vec3 RadiationSimulation::RadiationSource::drawRayDirection()
{
	return glm::vec3(this->rotation * glm::vec4(this->shape->drawRayDirection(), 0.f));
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
	if (max_energy_eV <= 0.f)
		max_energy_eV = spectrum_probabilities->max();
	else 
		if (max_energy_eV < spectrum_probabilities->max()) {
			std::cerr << "Max energy value is lower than the maximum value in the spectrum: " << max_energy_eV << " < " << spectrum_probabilities->max() << std::endl;
			throw std::runtime_error("Max energy value is lower than the maximum value in the spectrum");
		}

	this->hist_buffer = std::vector<float>(std::max<size_t>(max_energy_eV / 1e+3, 1), 0.0f);
	this->hist = RadFiled3D::HistogramVoxel(this->hist_buffer.size(), 1e+3, this->hist_buffer.data());
	this->hist.clear();
}

RadiationSimulation::XRaySpectrumSource::~XRaySpectrumSource()
{
	// No need to delete hist_buffer, unique_ptr will handle it
}

size_t RadiationSimulation::XRaySpectrumSource::getPossibilitiesCount() const
{
	return this->spectrum_probabilities->get_point_count();
}

float RadiationSimulation::XRaySpectrumSource::drawEnergy_eV()
{
	float energy_eV = 0.0f;
	size_t draw_count = 0;
	while(energy_eV <= 0.0f || energy_eV < this->energy_lower_cut_eV) {
		energy_eV = this->spectrum_probabilities->draw_sample();
		if(++draw_count > 10000) {
			throw std::runtime_error("Failed to draw a valid energy value");
		}
	}

	if (energy_eV > spectrum_probabilities->max())
		std::cout << "Energy value is higher than the maximum value in the spectrum: " << energy_eV << " > " << spectrum_probabilities->max() << std::endl;

	std::unique_lock lock(this->hist_mutex);
	this->hist.add_value(energy_eV);

	return energy_eV;
}

std::shared_ptr<Statistics::ProbabilityDensityFunction<float>> RadiationSimulation::SpectrumLoader::LoadSpectrum(const std::string& filename)
{
	std::ifstream file(filename);
	if (!file.is_open()) {
		throw std::runtime_error("Failed to open file");
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
	: rot_angle_distribution(
		  std::uniform_real_distribution<float>(
			  0.f,
			  1.f
		  )
	  ),
	  opening_angle_radians(glm::radians(opening_angle_deg)),
	  rand_engine(std::random_device{}())
{
}

glm::vec3 RadiationSimulation::ConeSourceShape::drawRayDirection()
{
	float theta = 2.0f * glm::pi<float>() * this->rot_angle_distribution(this->rand_engine); // Azimuthal angle
	float phi = acos(1.0f - this->rot_angle_distribution(this->rand_engine) * (1.0f - cos(this->opening_angle_radians))); // Polar angle

	return glm::vec3(
		std::sin(phi) * std::cos(theta),
		std::sin(phi) * std::sin(theta),
		-std::cos(phi)
	);
}

RadiationSimulation::RectangleSourceShape::RectangleSourceShape(const glm::vec2& size, float distance)
	: size(size),
	  distance(distance),
	  distribution_x(-0.5f * size.x, 0.5f * size.x),
	  distribution_y(-0.5f * size.y, 0.5f * size.y),
	  rand_engine(std::random_device{}())
{
	
}

glm::vec3 RadiationSimulation::RectangleSourceShape::drawRayDirection()
{
	float x = this->distribution_x(this->rand_engine);
	float y = this->distribution_y(this->rand_engine);

	return glm::normalize(glm::vec3(x, y, -this->distance));
}

RadiationSimulation::EllipsoidSourceShape::EllipsoidSourceShape(const glm::vec2& angles)
	: angles(
		glm::vec2(
			glm::radians(angles.x),
			glm::radians(angles.y)
		)
	),
	distribution(-0.5f, 0.5f),
	rand_engine(std::random_device{}())
{
}

glm::vec3 RadiationSimulation::EllipsoidSourceShape::drawRayDirection()
{
	// draw ellipsoid direction distribution from 2D-angles
	float x = this->distribution(this->rand_engine) * this->angles.x;
	float y = this->distribution(this->rand_engine) * this->angles.y;

	glm::vec3 direction = glm::vec3(0.f, 0.f, -1.0f);
	glm::quat rotation = glm::angleAxis(x, glm::vec3(1.0f, 0.0f, 0.0f)) * glm::angleAxis(y, glm::vec3(0.0f, 1.0f, 0.0f));

	return rotation * direction;
}
