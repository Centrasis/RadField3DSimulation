#include "Statistics.hpp"
#include <math.h>
#include <algorithm>
#include <numeric>


void Statistics::Variance::add(float value)
{
	this->count++;
	float delta = value - this->mean;
	this->mean += delta / this->count;
	float delta2 = value - this->mean;
	this->M2 += delta * delta2;
}

float Statistics::Variance::get_variance() const
{
	return (this->count > 0) ? this->M2 / this->count : 0.f;
}

float Statistics::Variance::get_mean() const
{
	return this->mean;
}

float Statistics::Variance::get_relative_error() const
{
	float mean = this->get_mean();
	return (mean != 0.f) ? this->get_standard_error() / mean : this->get_standard_error();
}

float Statistics::Variance::get_standard_error() const
{
	float variance = this->get_variance();
	if (variance <= 0.f)
		return 0.f;
	return (this->count > 0) ? sqrt(variance) : 1.f;
}

Statistics::HistogramMeanDistributionVariance::HistogramMeanDistributionVariance(float bin_width, size_t bins, size_t score_every_n)
	: bin_width(bin_width),
	  bins(bins),
	  score_every_n(score_every_n)
{
	this->reset();
}

void Statistics::HistogramMeanDistributionVariance::add(const RadFiled3D::HistogramVoxel<float>& vx)
{
	this->add_count++;
	if (this->add_count % this->score_every_n != 0)
		return;
	this->add_count = 0;
	this->count++;

	float sum = 0.f;
	for (size_t i = 0; i < vx.get_bins(); i++)
		sum += vx.get_histogram()[i];

	float intersectionSum = 0.f;
	float unionSum = 0.f;
	for (size_t i = 0; i < vx.get_bins(); ++i) {
		float hist_val = vx.get_histogram()[i] / sum;
		this->cumulative_histogram[i] += hist_val;
		float mean_hist_val = this->cumulative_histogram[i] / static_cast<float>(this->count);
		intersectionSum += std::min(mean_hist_val, hist_val);
		unionSum += std::max(mean_hist_val, hist_val);
	}

	this->cumulative_error += 1.f - ((unionSum > 0.f) ? intersectionSum / unionSum : 1.f);
}

void Statistics::HistogramMeanDistributionVariance::reset()
{
	this->add_count = 0;
	this->count = 0;
	this->cumulative_error = 0.f;
	this->cumulative_histogram = std::vector<float>(this->bins, 0.f);
}

float Statistics::HistogramMeanDistributionVariance::get_relative_error() const
{
	if (this->count < 2) {
		return 1.f; // Keine statistische Aussage, wenn weniger als zwei Histogramme vorhanden sind
	}

	return this->cumulative_error / static_cast<float>(this->count - 1);
}

Statistics::HistogramDistributionVariance::HistogramDistributionVariance(size_t bins, size_t score_every_n)
	: bins(bins),
	score_every_n(score_every_n),
	variances(bins)
{
	this->reset();
}

void Statistics::HistogramDistributionVariance::add(const RadFiled3D::HistogramVoxel<float>& vx)
{
	this->add_count++;
	if (this->add_count % this->score_every_n != 0)
		return;
	this->add_count = 0;
	this->count++;

	float sum = 0.f;
	for (size_t i = 0; i < vx.get_bins(); i++)
		sum += vx.get_histogram()[i];

	for (size_t i = 0; i < vx.get_bins(); i++) {
		float val = (sum > 0.f) ? vx.get_histogram()[i] / sum : 0.f;
		this->variances[i].add(val);
	}
}

void Statistics::HistogramDistributionVariance::reset()
{
	this->add_count = 0;
	this->count = 0;
	for (size_t i = 0; i < this->bins; i++)
		this->variances[i].reset();
}

float Statistics::HistogramDistributionVariance::get_variance() const
{
	if (this->count < 2) {
		return 0.25f; // Keine statistische Aussage, wenn weniger als zwei Histogramme vorhanden sind
	}

	float sum = 0.f;
	for (size_t i = 0; i < this->bins; i++)
		sum += this->variances[i].get_variance();

	return sum / static_cast<float>(this->bins);
}

float Statistics::HistogramDistributionVariance::get_relative_error() const
{
	if (this->count < 2) {
		return 1.f; // Keine statistische Aussage, wenn weniger als zwei Histogramme vorhanden sind
	}

	float sum = 0.f;
	for (size_t i = 0; i < this->bins; i++)
		sum += this->variances[i].get_variance();

	return (sum / static_cast<float>(this->bins)) * 4.f;
}


Statistics::VoxelSpectraVariance::VoxelSpectraVariance(size_t voxel_count, size_t bins, size_t score_every_n)
	: bins(bins),
	  voxel_count(voxel_count),
	  score_every_n(score_every_n),
	  add_counts(voxel_count, 0),
	  counts(voxel_count, 0),
	  means(voxel_count * bins, 0.f),
	  m2s(voxel_count * bins, 0.f)
{
}

void Statistics::VoxelSpectraVariance::add(size_t voxel_idx, const RadFiled3D::HistogramVoxel<double>& vx)
{
	this->add_counts[voxel_idx]++;
	if (this->add_counts[voxel_idx] % this->score_every_n != 0)
		return;
	this->add_counts[voxel_idx] = 0;
	const uint32_t count = ++this->counts[voxel_idx];

	double sum = 0.0;
	for (size_t i = 0; i < this->bins; i++)
		sum += vx.get_histogram()[i];

	float* mean = this->means.data() + voxel_idx * this->bins;
	float* m2 = this->m2s.data() + voxel_idx * this->bins;
	for (size_t i = 0; i < this->bins; i++) {
		const float val = (sum > 0.0) ? static_cast<float>(vx.get_histogram()[i] / sum) : 0.f;
		const float delta = val - mean[i];
		mean[i] += delta / static_cast<float>(count);
		m2[i] += delta * (val - mean[i]);
	}
}

void Statistics::VoxelSpectraVariance::reset()
{
	std::fill(this->add_counts.begin(), this->add_counts.end(), 0);
	std::fill(this->counts.begin(), this->counts.end(), 0);
	std::fill(this->means.begin(), this->means.end(), 0.f);
	std::fill(this->m2s.begin(), this->m2s.end(), 0.f);
}

float Statistics::VoxelSpectraVariance::get_relative_error(size_t voxel_idx) const
{
	const uint32_t count = this->counts[voxel_idx];
	if (count < 2)
		return 1.f; // Keine statistische Aussage, wenn weniger als zwei Histogramme vorhanden sind

	const float* m2 = this->m2s.data() + voxel_idx * this->bins;
	float sum = 0.f;
	for (size_t i = 0; i < this->bins; i++)
		sum += m2[i] / static_cast<float>(count);

	return (sum / static_cast<float>(this->bins)) * 4.f;
}
