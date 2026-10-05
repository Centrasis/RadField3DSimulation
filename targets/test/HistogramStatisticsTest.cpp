#include "Statistics.hpp"
#include <iostream>
#include "gtest/gtest.h"
#include <math.h>
#include <vector>
#include <Randomize.hh>

namespace {
	TEST(Variance, TestBounds) {
		Statistics::HistogramDistributionVariance variance(1, 1);
		float db[1] = { 1.f };

		radfiled3d::HistogramVoxel<float> histogram(static_cast<size_t>(1), 1.f, (float*)&db);
		variance.add(histogram);
		variance.add(histogram);
		float var = variance.get_variance();
		EXPECT_FLOAT_EQ(var, 0.f);
		variance.add(histogram);
		var = variance.get_variance();
		EXPECT_FLOAT_EQ(var, 0.f);

		variance.reset();
		variance.add(histogram);
		db[0] = 0.f;
		variance.add(histogram);
		var = variance.get_variance();
		EXPECT_FLOAT_EQ(var, 0.25f);

		db[0] = 0.5f;
		variance.add(histogram);
		var = variance.get_variance();
		EXPECT_TRUE(std::abs(var - 0.222f) < 0.001f);

		variance.reset();
		for (size_t i = 0; i < 100; i++) {
			db[0] = (i % 2 == 0) ? 0.f : 1.f;
			variance.add(histogram);
		}
		var = variance.get_variance();
		EXPECT_FLOAT_EQ(var, 0.25f);
	}

	TEST(RelError, TestBounds) {
		Statistics::HistogramDistributionVariance variance(1, 1);
		float db[1] = { 1.f };

		radfiled3d::HistogramVoxel<float> histogram(static_cast<size_t>(1), 1.f, (float*)&db);
		variance.add(histogram);
		variance.add(histogram);
		float err = variance.get_relative_error();
		EXPECT_FLOAT_EQ(err, 0.f);
		variance.add(histogram);
		err = variance.get_relative_error();
		EXPECT_FLOAT_EQ(err, 0.f);

		variance.reset();
		variance.add(histogram);
		db[0] = 0.f;
		variance.add(histogram);
		err = variance.get_relative_error();
		EXPECT_FLOAT_EQ(err, 1.f);

		db[0] = 0.5f;
		variance.add(histogram);
		err = variance.get_relative_error();
		EXPECT_FLOAT_EQ(err, 0.88888884f);

		variance.reset();
		for (size_t i = 0; i < 100; i++) {
			db[0] = (i % 2 == 0) ? 0.f : 1.f;
			variance.add(histogram);
		}
		err = variance.get_relative_error();
		EXPECT_FLOAT_EQ(err, 1.f);
	}

	TEST(RelError, Convergence) {
		Statistics::HistogramDistributionVariance variance(20, 1);
		float db[20] = { 0.f };

		radfiled3d::HistogramVoxel<float> histogram(static_cast<size_t>(20), 1.f, (float*)&db);
		G4Random::setTheSeed(1234);

		variance.add(histogram);
		float last_err = variance.get_relative_error();

		for (size_t epoch = 9; epoch <= 0; epoch--) {
			size_t repeats = (epoch <= 1) ? 5 : 1;
			for (size_t repeat = 0; repeat < repeats; repeat++) {
				for (size_t i = 0; i < 100; i++) {
					histogram.get_histogram()[G4RandFlat::shootInt(static_cast<long>(10 - epoch), static_cast<long>(10 + epoch + 1))]++;
				}

				variance.add(histogram);

				float err = variance.get_relative_error();
				EXPECT_GE(err, 0.f);
				EXPECT_LE(err, 1.f);

				if (err == 1.f) {
					EXPECT_FLOAT_EQ(err, last_err);
				}
				else {
					EXPECT_LE(err, last_err);
					last_err = err;
				}
			}
		}
	}

	TEST(RelError, Divergence) {
		Statistics::HistogramDistributionVariance variance(20, 1);
		float db[20] = { 0.f };

		radfiled3d::HistogramVoxel<float> histogram(static_cast<size_t>(20), 1.f, (float*)&db);
		G4Random::setTheSeed(1234);

		variance.add(histogram);
		variance.add(histogram);
		float last_err = 0.f;

		for (size_t epoch = 0; epoch < 10; epoch++) {
			size_t repeats = 2;
			for (size_t repeat = 0; repeat < repeats; repeat++) {
				for (size_t i = 0; i < 100; i++) {
					histogram.get_histogram()[G4RandFlat::shootInt(static_cast<long>(10 - epoch), static_cast<long>(10 + epoch + 1))]++;
				}
				variance.add(histogram);

				float err = variance.get_relative_error();
				EXPECT_GE(err, 0.f);
				EXPECT_LE(err, 1.f);
				EXPECT_GT(err, last_err * 0.7f);
				last_err = err;
			}
		}
	}
};

TEST(VoxelHistoryVariance, MatchesTheSampleStandardErrorOfTheMean) {
	// histories scoring 1, 0, 0, 1: mean 0.5, sample variance 1/3, standard error sqrt(1/12), relative 1/sqrt(3)
	Statistics::VoxelHistoryVariance variance(1);
	Statistics::VoxelHistoryVariance::History history;
	double sum = 0.0;
	for (double score : { 1.0, 0.0, 0.0, 1.0 }) {
		if (score > 0.0)
			history[0] += score;
		sum += score;
		variance.end_history(history);
	}
	EXPECT_TRUE(history.empty());
	EXPECT_DOUBLE_EQ(variance.get_sum_of_squares(0), 2.0);
	EXPECT_NEAR(Statistics::VoxelHistoryVariance::relative_error(sum, variance.get_sum_of_squares(0), 4), 1.0 / std::sqrt(3.0), 1e-12);
}

TEST(VoxelHistoryVariance, ScoresOfOneHistoryAreNotIndependent) {
	// every history crosses the voxel twice with 0.5 each: a total of 1 per history, so the flux is exact (R = 0);
	// squaring the single crossings would wrongly report an error
	Statistics::VoxelHistoryVariance variance(1);
	Statistics::VoxelHistoryVariance::History history;
	const size_t n = 1000;
	for (size_t i = 0; i < n; i++) {
		history[0] += 0.5;
		history[0] += 0.5;
		variance.end_history(history);
	}
	EXPECT_NEAR(Statistics::VoxelHistoryVariance::relative_error(static_cast<double>(n), variance.get_sum_of_squares(0), n), 0.0, 1e-12);
}

TEST(VoxelHistoryVariance, NothingScoredIsFullyUncertain) {
	EXPECT_EQ(Statistics::VoxelHistoryVariance::relative_error(0.0, 0.0, 1000), 1.0);
	EXPECT_EQ(Statistics::VoxelHistoryVariance::relative_error(1.0, 1.0, 1), 1.0);
}

TEST(VoxelHistoryVariance, PredictsTheSpreadOfIndependentRuns) {
	// a voxel hit with probability p per history, scoring an exponential path length: the predicted R must match the
	// relative spread of the mean over independent runs, and fall as 1 / sqrt(N)
	G4Random::setTheSeed(17);
	const double p = 0.02;
	const size_t histories = 20000, runs = 400;
	std::vector<double> means;
	double predicted = 0.0;
	for (size_t r = 0; r < runs; r++) {
		Statistics::VoxelHistoryVariance variance(1);
		Statistics::VoxelHistoryVariance::History history;
		double sum = 0.0;
		for (size_t i = 0; i < histories; i++) {
			if (G4UniformRand() < p) {
				const double score = -std::log(1.0 - G4UniformRand());
				history[0] += score;
				sum += score;
			}
			variance.end_history(history);
		}
		means.push_back(sum / histories);
		predicted += Statistics::VoxelHistoryVariance::relative_error(sum, variance.get_sum_of_squares(0), histories) / runs;
	}
	double mean = 0.0, spread = 0.0;
	for (double m : means)
		mean += m / runs;
	for (double m : means)
		spread += (m - mean) * (m - mean) / (runs - 1);
	const double observed = std::sqrt(spread) / mean;
	// analytic: R = sqrt((E[x^2] / (p E[x]^2) - 1) / N) = sqrt((2 / p - 1) / N) for unit exponential scores
	EXPECT_NEAR(predicted, std::sqrt((2.0 / p - 1.0) / histories), 0.02 * predicted);
	EXPECT_NEAR(observed, predicted, 0.1 * predicted);
}
