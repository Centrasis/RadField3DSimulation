#include "RadiationFieldDetector.hpp"
#include <cassert>
#include "RadiationFieldDetector.hpp"
#include <World.hpp>
#include <G4ios.hh>
#include "RadiationSource.hpp"


using namespace RadiationSimulation;

RadiationSimulation::RadiationDetector::RadiationDetector(size_t spectrum_bins, double spectrum_bin_width, std::shared_ptr<RadFiled3D::VoxelGridBuffer> voxels, const glm::uvec3& coord)
	: energy(voxels->get_voxel<RadFiled3D::ScalarVoxel<float>>("energy", coord.x, coord.y, coord.z)),
	  hits(voxels->get_voxel<RadFiled3D::ScalarVoxel<float>>("hits", coord.x, coord.y, coord.z)),
	  direction(voxels->get_voxel<RadFiled3D::ScalarVoxel<glm::vec3>>("direction", coord.x, coord.y, coord.z)),
	  spectrum(voxels->get_voxel<RadFiled3D::HistogramVoxel>("spectrum", coord.x, coord.y, coord.z)),
	  energy_squared_sum(0.f)
{

}

void RadiationSimulation::RadiationDetector::normalize(size_t particle_count)
{
	this->energy.normalize(particle_count);
	this->hits.normalize(particle_count);
	this->direction.normalize(particle_count);
	this->spectrum.normalize(particle_count);
}


RadiationSimulation::RadiationFieldDetector::RadiationFieldDetector(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, size_t spectra_bins, double spectra_bin_width)
	: field(std::make_shared<RadFiled3D::CartesianRadiationField>(radiation_field_dimensions, radiation_field_voxel_dimensions)),
	  spectra_bins(spectra_bins),
	  spectra_bin_width(spectra_bin_width)
{
	assert(this->field->get_voxel_counts().x > 0 && this->field->get_voxel_counts().y > 0 && this->field->get_voxel_counts().z > 0);
}

std::shared_ptr<RadFiled3D::IRadiationField> RadiationSimulation::RadiationFieldDetector::evaluate()
{
	this->evaluate_field();
	return this->field;
}

void RadiationSimulation::SpectrumVoxel::normalize(size_t particle_count)
{
	float sum = 0.f;
	for (size_t i = 0; i < this->voxel.get_bins(); i++)
		sum += (&this->voxel.get_data())[i];
	for (size_t i = 0; i < this->voxel.get_bins(); i++)
		(&this->voxel.get_data())[i] /= sum;
}
