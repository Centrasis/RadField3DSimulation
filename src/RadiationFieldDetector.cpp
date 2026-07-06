#include "RadiationFieldDetector.hpp"
#include <cassert>
#include "RadiationFieldDetector.hpp"
#include <World.hpp>
#include <G4ios.hh>
#include "RadiationSource.hpp"


using namespace RadiationSimulation;

RadiationSimulation::RadiationFieldDetector::RadiationFieldDetector(const glm::vec3& radiation_field_dimensions, const glm::vec3& radiation_field_voxel_dimensions, size_t spectra_bins, double spectra_bin_width)
	: field(std::make_shared<RadFiled3D::CartesianRadiationField>(radiation_field_dimensions, radiation_field_voxel_dimensions)),
	  spectra_bins(spectra_bins),
	  spectra_bin_width(spectra_bin_width)
{
	assert(this->field->get_voxel_counts().x > 0 && this->field->get_voxel_counts().y > 0 && this->field->get_voxel_counts().z > 0);
}

std::shared_ptr<RadFiled3D::IRadiationField> RadiationSimulation::RadiationFieldDetector::evaluate()
{
	// The scoring field accumulates in double; every consumer gets the normalized fp32 copy.
	this->evaluate_field();
	return this->get_normalized_field_copy();
}

void RadiationSimulation::SpectrumVoxel::normalize(size_t particle_count)
{
	float sum = 0.f;
	for (size_t i = 0; i < this->voxel.get_bins(); i++)
		sum += (&this->voxel.get_data())[i];
	for (size_t i = 0; i < this->voxel.get_bins(); i++)
		(&this->voxel.get_data())[i] /= sum;
}
