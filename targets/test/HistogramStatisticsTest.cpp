#include "Statistics.hpp"
#include <iostream>
#include "gtest/gtest.h"
#include <math.h>
#include <vector>
#include <random>

namespace {
	TEST(Variance, TestBounds) {
		Statistics::HistogramDistributionVariance variance(1, 1);
		float db[1] = { 1.f };

		RadFiled3D::HistogramVoxel histogram(static_cast<size_t>(1), 1.f, (float*)&db);
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

		RadFiled3D::HistogramVoxel histogram(static_cast<size_t>(1), 1.f, (float*)&db);
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

		RadFiled3D::HistogramVoxel histogram(static_cast<size_t>(20), 1.f, (float*)&db);
		std::random_device rd;
		std::mt19937 gen(rd());

		variance.add(histogram);
		float last_err = variance.get_relative_error();

		for (size_t epoch = 9; epoch <= 0; epoch--) {
			std::uniform_int_distribution<size_t> dis(10 - epoch, 10 + epoch);
			size_t repeats = (epoch <= 1) ? 5 : 1;
			for (size_t repeat = 0; repeat < repeats; repeat++) {
				for (size_t i = 0; i < 100; i++) {
					histogram.get_histogram()[dis(gen)]++;
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

		RadFiled3D::HistogramVoxel histogram(static_cast<size_t>(20), 1.f, (float*)&db);
		std::random_device rd;
		std::mt19937 gen(rd());

		variance.add(histogram);
		variance.add(histogram);
		float last_err = 0.f;

		for (size_t epoch = 0; epoch < 10; epoch++) {
			std::uniform_int_distribution<size_t> dis(10 - epoch, 10 + epoch);
			size_t repeats = 2;
			for (size_t repeat = 0; repeat < repeats; repeat++) {
				for (size_t i = 0; i < 100; i++) {
					histogram.get_histogram()[dis(gen)]++;
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
