#include "VMFTrainer.hpp"
#include <gtest/gtest.h>
#include <Randomize.hh>
#include <glm/geometric.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>


using namespace RadiationSimulation;

namespace {
	// vMF sample around mean with concentration kappa (Wood 1994, 3D case)
	glm::vec3 sample_vmf(const glm::vec3& mean, double kappa)
	{
		const double u = G4UniformRand();
		const double w = 1.0 + std::log(u + (1.0 - u) * std::exp(-2.0 * kappa)) / kappa;
		const double phi = 2.0 * std::numbers::pi * G4UniformRand();
		const glm::vec3 helper = (std::abs(mean.x) < 0.9f) ? glm::vec3(1.f, 0.f, 0.f) : glm::vec3(0.f, 1.f, 0.f);
		const glm::vec3 b = glm::normalize(glm::cross(mean, helper));
		const glm::vec3 c = glm::cross(mean, b);
		const double s = std::sqrt(std::max(0.0, 1.0 - w * w));
		return glm::normalize(static_cast<float>(w) * mean + static_cast<float>(s * std::cos(phi)) * b + static_cast<float>(s * std::sin(phi)) * c);
	}

	struct Lobe { float weight; glm::vec3 mean; float kappa; };

	std::vector<Lobe> lobes_of(const VMFTrainer& trainer, size_t voxel)
	{
		std::vector<float> raw(trainer.get_lobes() * VMFTrainer::VALUES_PER_LOBE);
		trainer.write_lobes(voxel, raw.data());
		std::vector<Lobe> out;
		for (size_t k = 0; k < trainer.get_lobes(); k++)
			out.push_back({ raw[k * 5], glm::vec3(raw[k * 5 + 1], raw[k * 5 + 2], raw[k * 5 + 3]), raw[k * 5 + 4] });
		std::sort(out.begin(), out.end(), [](const Lobe& a, const Lobe& b) { return a.weight > b.weight; });
		return out;
	}
}

TEST(VMFTrainer, RecoversAMixtureInDoublingPasses) {
	G4Random::setTheSeed(7);
	const glm::vec3 beam = glm::normalize(glm::vec3(1.f, 0.2f, 0.f));
	const glm::vec3 room = glm::normalize(glm::vec3(-0.3f, 1.f, 0.5f));
	// the initial guess points elsewhere on purpose
	VMFTrainer trainer(1, 3, [](size_t) { return glm::vec3(0.f, 0.f, 1.f); });
	size_t pass_length = 1000;
	for (int pass = 0; pass < 8; pass++, pass_length *= 2) {
		for (size_t i = 0; i < pass_length; i++)
			trainer.add(0, (G4UniformRand() < 0.7) ? sample_vmf(beam, 200.0) : sample_vmf(room, 10.0));
		trainer.m_step();
	}
	const std::vector<Lobe> lobes = lobes_of(trainer, 0);
	EXPECT_NEAR(lobes[0].weight, 0.7, 0.03);
	EXPECT_GT(glm::dot(lobes[0].mean, beam), 0.999f);
	EXPECT_NEAR(lobes[0].kappa, 200.0, 20.0);
	// the wide component may be split over the remaining two lobes: together they carry 30 % and point to the room
	const glm::vec3 rest = lobes[1].weight * lobes[1].mean + lobes[2].weight * lobes[2].mean;
	EXPECT_NEAR(lobes[1].weight + lobes[2].weight, 0.3, 0.03);
	EXPECT_GT(glm::dot(glm::normalize(rest), room), 0.95f);
	float total = 0.f;
	for (const Lobe& l : lobes)
		total += l.weight;
	EXPECT_NEAR(total, 1.f, 1e-5f);
}

TEST(VMFTrainer, EightLobesRecoverTheSameMixture) {
	G4Random::setTheSeed(9);
	const glm::vec3 beam = glm::normalize(glm::vec3(1.f, 0.2f, 0.f));
	const glm::vec3 room = glm::normalize(glm::vec3(-0.3f, 1.f, 0.5f));
	VMFTrainer trainer(1, 8, [](size_t) { return glm::vec3(0.f, 0.f, 1.f); });
	size_t pass_length = 1000;
	for (int pass = 0; pass < 8; pass++, pass_length *= 2) {
		for (size_t i = 0; i < pass_length; i++)
			trainer.add(0, (G4UniformRand() < 0.7) ? sample_vmf(beam, 200.0) : sample_vmf(room, 10.0));
		trainer.m_step();
	}
	const std::vector<Lobe> lobes = lobes_of(trainer, 0);
	// the density must match the mixture, however the lobes split it
	double err = 0.0, norm = 0.0;
	for (const glm::vec3& d : { beam, room, glm::normalize(beam + room), glm::vec3(0.f, 0.f, -1.f), glm::normalize(glm::vec3(1.f, 0.25f, 0.05f)) }) {
		auto vmf = [](const glm::vec3& d, const glm::vec3& mu, double k) { return k / (2 * std::numbers::pi * (1 - std::exp(-2 * k))) * std::exp(k * (glm::dot(d, mu) - 1)); };
		const double truth = 0.7 * vmf(d, beam, 200.0) + 0.3 * vmf(d, room, 10.0);
		double fitted = 0.0;
		for (const Lobe& l : lobes)
			if (l.weight > 0.f)
				fitted += l.weight * vmf(d, l.mean, l.kappa);
		err += std::abs(fitted - truth);
		norm += truth;
	}
	EXPECT_LT(err / norm, 0.1);
	float total = 0.f;
	for (const Lobe& l : lobes)
		total += l.weight;
	EXPECT_NEAR(total, 1.f, 1e-5f);
}

TEST(VMFTrainer, OneSourceEndsAsOneLobe) {
	G4Random::setTheSeed(8);
	const glm::vec3 source = glm::normalize(glm::vec3(0.3f, -1.f, 0.2f));
	VMFTrainer trainer(1, 3, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); });
	size_t pass_length = 1000;
	for (int pass = 0; pass < 8; pass++, pass_length *= 2) {
		for (size_t i = 0; i < pass_length; i++)
			trainer.add(0, sample_vmf(source, 50.0));
		trainer.m_step();
	}
	const std::vector<Lobe> lobes = lobes_of(trainer, 0);
	EXPECT_GT(lobes[0].weight, 0.99f);
	EXPECT_GT(glm::dot(lobes[0].mean, source), 0.999f);
	EXPECT_NEAR(lobes[0].kappa, 50.0, 5.0);
	EXPECT_LT(lobes[1].weight + lobes[2].weight, 0.01f);

	// fixed layout: all slots present, strongest first, unused slots entirely zero
	std::vector<float> raw(15);
	trainer.write_lobes(0, raw.data());
	EXPECT_EQ(raw[0], lobes[0].weight);
	EXPECT_GE(raw[0], raw[5]);
	EXPECT_GE(raw[5], raw[10]);
	for (size_t k = 1; k < 3; k++)
		if (raw[k * 5] == 0.f)
			for (size_t c = 0; c < 5; c++)
				EXPECT_EQ(raw[k * 5 + c], 0.f);
}

TEST(VMFTrainer, VoxelsWithoutEnoughSamplesStoreTheUniformDistribution) {
	const std::vector<float> uniform = { 1.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };
	VMFTrainer trainer(3, 3, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); });
	std::vector<float> raw(15);
	// before any sample every voxel is uniform
	trainer.write_lobes(0, raw.data());
	EXPECT_EQ(raw, uniform);

	G4Random::setTheSeed(12);
	const glm::vec3 source = glm::normalize(glm::vec3(0.f, 1.f, 1.f));
	for (int i = 0; i < 19; i++)
		trainer.add(1, sample_vmf(source, 50.0));
	for (int i = 0; i < 20; i++)
		trainer.add(2, sample_vmf(source, 50.0));
	trainer.m_step();
	trainer.write_lobes(0, raw.data());
	EXPECT_EQ(raw, uniform);
	trainer.write_lobes(1, raw.data());                   // 19 of the 20 effective samples needed
	EXPECT_EQ(raw, uniform);
	trainer.write_lobes(2, raw.data());                   // enough: fitted towards the source
	EXPECT_NEAR(raw[0] + raw[5] + raw[10], 1.f, 1e-6f);
	EXPECT_GT(glm::dot(glm::vec3(raw[1], raw[2], raw[3]), source), 0.9f);

	// unequal weights count less: 30 samples, one of them carrying almost all the weight
	VMFTrainer weighted(1, 2, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); });
	weighted.add(0, source, 100.0);
	for (int i = 0; i < 29; i++)
		weighted.add(0, sample_vmf(source, 50.0), 1.0);
	weighted.m_step();
	weighted.write_lobes(0, raw.data());
	EXPECT_EQ(std::vector<float>(raw.begin(), raw.begin() + 10), std::vector<float>(uniform.begin(), uniform.begin() + 10));

	// the threshold is a parameter; 0 stores the fit of a single sample
	VMFTrainer eager(1, 2, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); }, 0.0);
	eager.add(0, source);
	eager.m_step();
	eager.write_lobes(0, raw.data());
	EXPECT_GT(glm::dot(glm::vec3(raw[1], raw[2], raw[3]), source), 0.99f);

	EXPECT_THROW(VMFTrainer(1, 0, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); }), std::invalid_argument);
}

TEST(VMFTrainer, WeightedSamplesCountLikeRepeatedOnes) {
	const glm::vec3 a = glm::normalize(glm::vec3(1.f, 0.2f, 0.f)), b = glm::normalize(glm::vec3(-0.3f, 1.f, 0.5f));
	VMFTrainer weighted(1, 2, [](size_t) { return glm::vec3(0.f, 0.f, 1.f); });
	VMFTrainer repeated(1, 2, [](size_t) { return glm::vec3(0.f, 0.f, 1.f); });
	G4Random::setTheSeed(11);
	for (int pass = 0; pass < 4; pass++) {
		for (int i = 0; i < 2000; i++) {
			const glm::vec3 d = sample_vmf((i % 3 == 0) ? b : a, 80.0);
			weighted.add(0, d, 2.0);
			repeated.add(0, d);
			repeated.add(0, d);
		}
		weighted.m_step();
		repeated.m_step();
	}
	std::vector<float> w(10), r(10);
	weighted.write_lobes(0, w.data());
	repeated.write_lobes(0, r.data());
	for (size_t i = 0; i < w.size(); i++)
		EXPECT_NEAR(w[i], r[i], 1e-4f * std::max(1.f, std::abs(r[i])));
	// a zero weight adds nothing: the voxel stays uniform
	VMFTrainer empty(1, 2, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); });
	empty.add(0, a, 0.0);
	empty.m_step();
	empty.write_lobes(0, w.data());
	EXPECT_EQ(w, (std::vector<float>{ 1.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f }));
}

TEST(VMFTrainer, CompactSourceAndWideBackgroundStaySeparate) {
	// the situation behind a shield: most radiation from the patient in a compact lobe, a few percent diffuse room
	// scatter around it; merging them would make one medium lobe and lose the background
	G4Random::setTheSeed(12);
	const glm::vec3 source = glm::normalize(glm::vec3(1.f, 0.1f, 0.f));
	VMFTrainer trainer(1, 4, [](size_t) { return glm::vec3(0.f, 0.f, 1.f); });
	size_t pass_length = 2000;
	for (int pass = 0; pass < 7; pass++, pass_length *= 2) {
		for (size_t i = 0; i < pass_length; i++)
			trainer.add(0, (G4UniformRand() < 0.92) ? sample_vmf(source, 40.0) : sample_vmf(glm::normalize(glm::vec3(0.4f, 1.f, 0.f)), 1.0));
		trainer.m_step();
	}
	const std::vector<Lobe> lobes = lobes_of(trainer, 0);
	EXPECT_GT(glm::dot(lobes[0].mean, source), 0.99f);
	EXPECT_GT(lobes[0].kappa, 25.f);   // stays compact
	float wide = 0.f;
	for (const Lobe& l : lobes)
		if (l.weight > 0.f && l.kappa < 5.f)
			wide += l.weight;
	EXPECT_GT(wide, 0.03f);            // the background keeps its own wide lobe(s)
}

TEST(VMFTrainer, AcceptsMoreThanEightLobes) {
	G4Random::setTheSeed(11);
	const glm::vec3 source = glm::normalize(glm::vec3(-0.2f, 0.4f, 1.f));
	VMFTrainer trainer(1, 20, [](size_t) { return glm::vec3(1.f, 0.f, 0.f); });
	ASSERT_EQ(trainer.get_lobes(), 20u);
	size_t pass_length = 1000;
	for (int pass = 0; pass < 8; pass++, pass_length *= 2) {
		for (size_t i = 0; i < pass_length; i++)
			trainer.add(0, sample_vmf(source, 50.0));
		trainer.m_step();
	}
	const std::vector<Lobe> lobes = lobes_of(trainer, 0);
	ASSERT_EQ(lobes.size(), 20u);
	float total = 0.f;
	for (const Lobe& l : lobes)
		total += l.weight;
	EXPECT_NEAR(total, 1.f, 1e-5f);
	EXPECT_GT(glm::dot(lobes[0].mean, source), 0.99f);
	for (size_t k = 1; k < lobes.size(); k++)
		EXPECT_GE(lobes[k - 1].weight, lobes[k].weight);
}
