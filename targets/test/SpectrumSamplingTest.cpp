#include "RadiationSource.hpp"
#include "ProbabilityFunctions.hpp"
#include <gtest/gtest.h>
#include <Randomize.hh>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <thread>
#include <vector>


using namespace RadiationSimulation;

namespace {
	// Evenly spaced "random" numbers (k + 0.5) / n: sampling with them reproduces the distribution exactly.
	UniformRandom stratified(size_t n)
	{
		auto k = std::make_shared<size_t>(0);
		return [k, n] { return (static_cast<double>((*k)++ % n) + 0.5) / static_cast<double>(n); };
	}

	Statistics::ProbabilityDensityFunction<float> three_bins()
	{
		return Statistics::ProbabilityDensityFunction<float>({ { 10.f, 1.f }, { 20.f, 1.f }, { 30.f, 2.f } });
	}
}

TEST(SpectrumSampling, BinsAreCentredOnTheListedEnergies) {
	const auto pdf = three_bins();
	EXPECT_FLOAT_EQ(pdf.min(), 5.f);
	EXPECT_FLOAT_EQ(pdf.max(), 35.f);
	EXPECT_DOUBLE_EQ(pdf.cdf(5.0), 0.0);
	EXPECT_DOUBLE_EQ(pdf.cdf(15.0), 0.25);
	EXPECT_DOUBLE_EQ(pdf.cdf(25.0), 0.5);
	EXPECT_DOUBLE_EQ(pdf.cdf(35.0), 1.0);
	EXPECT_FLOAT_EQ(pdf.quantile(0.125), 10.f);
	EXPECT_FLOAT_EQ(pdf.quantile(0.375), 20.f);
	EXPECT_FLOAT_EQ(pdf.quantile(0.75), 30.f);
	EXPECT_FLOAT_EQ(pdf.quantile(1.0), 35.f);
	EXPECT_DOUBLE_EQ(pdf.probability(25.0, 35.0), 0.5);
}

TEST(SpectrumSampling, LastBinIsSpreadNotASpike) {
	const auto pdf = three_bins();
	const size_t n = 10000;
	const UniformRandom uniform = stratified(n);
	size_t upper_half_of_last_bin = 0, at_listed_energy = 0;
	for (size_t i = 0; i < n; i++) {
		const float e = pdf.sample(uniform());
		EXPECT_GE(e, 5.f);
		EXPECT_LT(e, 35.f);
		upper_half_of_last_bin += (e > 30.f);
		at_listed_energy += (e == 30.f);
	}
	EXPECT_EQ(upper_half_of_last_bin, n / 4);
	EXPECT_LE(at_listed_energy, 1u);
}

TEST(SpectrumSampling, ZeroWeightRowsAtTheEndsAreNeverSampled) {
	const Statistics::ProbabilityDensityFunction<float> pdf({ { 0.f, 0.f }, { 1.f, 0.f }, { 2.f, 5.f }, { 3.f, 5.f }, { 4.f, 0.f }, { 5.f, 0.f } });
	EXPECT_FLOAT_EQ(pdf.min(), 1.5f);
	EXPECT_FLOAT_EQ(pdf.max(), 3.5f);
	EXPECT_FLOAT_EQ(pdf.quantile(0.0), 1.5f);
	EXPECT_FLOAT_EQ(pdf.quantile(1.0), 3.5f);
}

TEST(SpectrumSampling, MidBinSpectrumKeepsItsMeanEnergy) {
	// SpekPy-like: 0.5 keV bins listed by their centres 10.25, 10.75, ... keV
	std::vector<std::pair<float, float>> points;
	double weighted_sum = 0.0, weight_total = 0.0;
	for (int i = 0; i < 100; i++) {
		const float centre = 10250.f + 500.f * static_cast<float>(i);
		const float weight = static_cast<float>(1 + (i * 37) % 11);
		points.push_back({ centre, weight });
		weighted_sum += centre * weight;
		weight_total += weight;
	}
	const Statistics::ProbabilityDensityFunction<float> pdf(points);
	EXPECT_FLOAT_EQ(pdf.max(), 10250.f + 500.f * 99.f + 250.f);

	const size_t n = 200000;
	const UniformRandom uniform = stratified(n);
	double sum = 0.0;
	for (size_t i = 0; i < n; i++)
		sum += pdf.sample(uniform());
	// the former lower-edge convention was biased by +250 eV
	EXPECT_NEAR(sum / n, weighted_sum / weight_total, 1.0);
}

TEST(SpectrumSampling, LowerCutIsAppliedExactly) {
	XRaySpectrumSource source(std::make_shared<Statistics::ProbabilityDensityFunction<float>>(std::vector<std::pair<float, float>>{ { 10.f, 1.f }, { 20.f, 1.f }, { 30.f, 2.f } }), std::make_unique<ConeSourceShape>(5.f), 20.f, 40000.f);
	const size_t n = 10000;
	const UniformRandom uniform = stratified(n);
	size_t below_25 = 0;
	for (size_t i = 0; i < n; i++) {
		const float e = source.drawEnergy_eV(uniform);
		EXPECT_GE(e, 20.f);
		below_25 += (e < 25.f);
	}
	// [20, 25) holds 0.125 of the 0.625 above the cut
	EXPECT_EQ(below_25, n / 5);

	EXPECT_THROW(XRaySpectrumSource(std::make_shared<Statistics::ProbabilityDensityFunction<float>>(std::vector<std::pair<float, float>>{ { 10.f, 1.f } }), std::make_unique<ConeSourceShape>(5.f), 20.f, 40000.f), std::runtime_error);
}

TEST(SpectrumSampling, PhotonsAboveFiveKeVAreNeverCut) {
	auto pdf = std::make_shared<Statistics::ProbabilityDensityFunction<float>>(std::vector<std::pair<float, float>>{ { 4000.f, 1.f }, { 30000.f, 1.f } });
	EXPECT_NO_THROW(XRaySpectrumSource(pdf, std::make_unique<ConeSourceShape>(5.f), XRaySpectrumSource::max_energy_lower_cut_eV, 50000.f));
	EXPECT_THROW(XRaySpectrumSource(pdf, std::make_unique<ConeSourceShape>(5.f), 5001.f, 50000.f), std::invalid_argument);
	EXPECT_THROW(XRaySpectrumSource(pdf, std::make_unique<ConeSourceShape>(5.f), -1.f, 50000.f), std::invalid_argument);

	// the default keeps the whole spectrum
	XRaySpectrumSource source(pdf, std::make_unique<ConeSourceShape>(5.f));
	const UniformRandom uniform = stratified(1000);
	float lowest = 1e9f;
	for (int i = 0; i < 1000; i++)
		lowest = std::min(lowest, source.drawEnergy_eV(uniform));
	EXPECT_LT(lowest, 5000.f);
}

TEST(SpectrumSampling, LoaderReadsBinCentres) {
	const std::string path = (std::filesystem::temp_directory_path() / "rf3_spectrum_test.csv").string();
	{
		std::ofstream csv(path);
		csv << "Energy[keV]    Fluence[]\n10.25, 1\n10.75, 3\n11.25, 0\n";
	}
	const auto pdf = SpectrumLoader::LoadSpectrum(path);
	EXPECT_FLOAT_EQ(pdf->min(), 10000.f);
	EXPECT_FLOAT_EQ(pdf->max(), 11000.f);
	EXPECT_DOUBLE_EQ(pdf->probability(10500.0, 11000.0), 0.75);
	std::remove(path.c_str());
}

TEST(SpectrumSampling, SameGeant4SeedDrawsTheSameEnergies) {
	XRaySpectrumSource source(std::make_shared<Statistics::ProbabilityDensityFunction<float>>(std::vector<std::pair<float, float>>{ { 10.f, 1.f }, { 20.f, 1.f }, { 30.f, 2.f } }), std::make_unique<ConeSourceShape>(5.f), 0.f, 40000.f);
	const UniformRandom uniform = [] { return G4UniformRand(); };
	auto draw = [&](long seed) {
		G4Random::setTheSeed(seed);
		std::vector<float> energies;
		for (int i = 0; i < 100; i++)
			energies.push_back(source.drawEnergy_eV(uniform));
		return energies;
	};
	EXPECT_EQ(draw(42), draw(42));
	EXPECT_NE(draw(42), draw(43));
}

namespace {
	XRaySpectrumSource delta_source(float energy_eV)
	{
		return XRaySpectrumSource(std::make_shared<Statistics::ProbabilityDensityFunction<float>>(std::vector<std::pair<float, float>>{ { energy_eV, 1.f } }), std::make_unique<ConeSourceShape>(5.f), 0.f, 40000.f);
	}
}

TEST(TubeSpectrum, CountsPastTheFloatLimitAndBinsByCentre) {
	// a float bin stops growing at 2^24 counts; the counters must not
	XRaySpectrumSource source = delta_source(10400.f);
	const UniformRandom uniform = [] { return 0.5; };
	const uint64_t n = (uint64_t(1) << 24) + 1000;
	for (uint64_t i = 0; i < n; i++)
		source.drawEnergy_eV(uniform);
	const std::vector<uint64_t> counts = source.getGeneratedCounts();
	ASSERT_EQ(counts.size(), 40u);
	EXPECT_FLOAT_EQ(source.getGeneratedSpectrumBinWidth_eV(), 1000.f);
	// bin i is centred at i keV: 10.4 keV lies in bin 10
	EXPECT_EQ(counts[10], n);

	XRaySpectrumSource upper = delta_source(10600.f);
	upper.drawEnergy_eV(uniform);
	EXPECT_EQ(upper.getGeneratedCounts()[11], 1u);
}

TEST(TubeSpectrum, ConcurrentDrawsAreAllCounted) {
	XRaySpectrumSource source(std::make_shared<Statistics::ProbabilityDensityFunction<float>>(std::vector<std::pair<float, float>>{ { 10000.f, 1.f }, { 20000.f, 1.f }, { 30000.f, 2.f } }), std::make_unique<ConeSourceShape>(5.f), 0.f, 40000.f);
	const size_t threads = 8, per_thread = 200000;
	std::vector<std::thread> workers;
	for (size_t t = 0; t < threads; t++)
		workers.emplace_back([&source] {
			const UniformRandom uniform = stratified(per_thread);
			for (size_t i = 0; i < per_thread; i++)
				source.drawEnergy_eV(uniform);
		});
	for (auto& worker : workers)
		worker.join();
	uint64_t total = 0;
	for (uint64_t count : source.getGeneratedCounts())
		total += count;
	EXPECT_EQ(total, threads * per_thread);
}
