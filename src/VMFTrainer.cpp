#include "VMFTrainer.hpp"
#include <glm/geometric.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>


using namespace RadiationSimulation;

namespace {
	constexpr float INITIAL_KAPPA = 2.f;
	constexpr double MAX_KAPPA = 1e4;
	constexpr double MIN_LOBE_SHARE = 1e-4;   // below this share of the voxel's counts a lobe keeps its mean and kappa

	// log of the vMF normalization kappa / (4 pi sinh kappa), written for exp(kappa * (mu . d - 1)); 1 / (4 pi) at 0
	double log_vmf_norm(double kappa)
	{
		if (kappa < 1e-4)
			return -std::log(4.0 * std::numbers::pi);
		return std::log(kappa) - std::log(2.0 * std::numbers::pi) - std::log1p(-std::exp(-2.0 * kappa));
	}

	// mean resultant length of a vMF lobe, coth(kappa) - 1 / kappa
	double resultant_from_kappa(double kappa)
	{
		if (kappa < 1e-3)
			return kappa / 3.0 - kappa * kappa * kappa / 45.0;
		return 1.0 / std::tanh(kappa) - 1.0 / kappa;
	}

	// kappa of a vMF lobe from its mean resultant length: the approximation of Banerjee et al. (2005), which
	// overestimates kappa by a few percent in the mid range, refined by Newton steps (as RadFiled3D's merge does)
	double kappa_from_resultant(double r)
	{
		r = std::clamp(r, 1e-6, 0.99999);
		double kappa = r * (3.0 - r * r) / (1.0 - r * r);
		for (int i = 0; i < 4; i++) {
			const double a = resultant_from_kappa(kappa);
			const double da = 1.0 - a * a - 2.0 * a / kappa;
			if (!(da > 0.0))
				break;
			const double next = kappa - (a - r) / da;
			kappa = (next > 0.0) ? next : kappa * 0.5;
		}
		return std::min(kappa, MAX_KAPPA);
	}

	// Lobes of clearly different width describe different things (e.g. a compact source and the diffuse room
	// background) and are never merged, however close their means are.
	constexpr double MAX_MERGE_KAPPA_RATIO = 4.0;

	// Merges lobes that describe the same source: similar widths (kappa ratio at most MAX_MERGE_KAPPA_RATIO) and an
	// angle between their means below the narrower lobe's width (1 / sqrt(kappa)). Moment-preserving; the freed lobe
	// gets weight 0 and is not used any more.
	void merge_overlapping(float* lobes, uint32_t count, size_t stride)
	{
		while (true) {
			double best = 0.0;
			int bi = -1, bj = -1;
			for (uint32_t i = 0; i < count; i++) {
				const float* a = lobes + i * stride;
				if (a[0] <= 0.f)
					continue;
				for (uint32_t j = i + 1; j < count; j++) {
					const float* b = lobes + j * stride;
					if (b[0] <= 0.f)
						continue;
					const double k_min = std::max(1e-6, static_cast<double>(std::min(a[4], b[4])));
					const double k_max = std::max(1e-6, static_cast<double>(std::max(a[4], b[4])));
					if (k_max > MAX_MERGE_KAPPA_RATIO * k_min)
						continue;
					const double cosine = std::clamp(static_cast<double>(a[1] * b[1] + a[2] * b[2] + a[3] * b[3]), -1.0, 1.0);
					const double width = 1.0 / std::sqrt(k_max);
					const double overlap = width - std::acos(cosine);
					if (overlap > best) {
						best = overlap;
						bi = static_cast<int>(i);
						bj = static_cast<int>(j);
					}
				}
			}
			if (bi < 0)
				return;
			float* a = lobes + bi * stride;
			float* b = lobes + bj * stride;
			const double wa = a[0], wb = b[0], ra = resultant_from_kappa(a[4]), rb = resultant_from_kappa(b[4]);
			double r[3];
			for (int c = 0; c < 3; c++)
				r[c] = (wa * ra * a[1 + c] + wb * rb * b[1 + c]) / (wa + wb);
			const double length = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
			a[0] = static_cast<float>(wa + wb);
			if (length > 0.0)
				for (int c = 0; c < 3; c++)
					a[1 + c] = static_cast<float>(r[c] / length);
			a[4] = static_cast<float>(kappa_from_resultant(length));
			b[0] = 0.f;
		}
	}
}

VMFTrainer::VMFTrainer(size_t voxel_count, uint32_t lobes, const std::function<glm::vec3(size_t)>& initial_direction)
	: voxel_count(voxel_count),
	  lobes(lobes),
	  model(voxel_count * lobes * VALUES_PER_LOBE, 0.f),
	  log_norm(voxel_count * lobes, 0.f),
	  current(voxel_count * lobes * STATS_PER_LOBE, 0.0),
	  previous(voxel_count * lobes * STATS_PER_LOBE, 0.0)
{
	if (lobes < 1 || lobes > 8)
		throw std::invalid_argument("VMFTrainer: the number of lobes must lie in [1, 8]");
	for (size_t v = 0; v < voxel_count; v++) {
		glm::vec3 a = initial_direction(v);
		a = (glm::length(a) > 0.f) ? glm::normalize(a) : glm::vec3(1.f, 0.f, 0.f);
		const glm::vec3 helper = (std::abs(a.x) < 0.9f) ? glm::vec3(1.f, 0.f, 0.f) : glm::vec3(0.f, 1.f, 0.f);
		const glm::vec3 b = glm::normalize(glm::cross(a, helper));
		const glm::vec3 c = glm::cross(a, b);
		// distinct start directions (identical lobes would never separate in EM): the 6 axes of the voxel's frame,
		// then the 8 diagonals
		const float d = 1.f / std::sqrt(3.f);
		const glm::vec3 seeds[14] = { a, -a, b, -b, c, -c,
			d * (a + b + c), d * (a + b - c), d * (a - b + c), d * (a - b - c),
			d * (-a + b + c), d * (-a + b - c), d * (-a - b + c), d * (-a - b - c) };
		for (uint32_t k = 0; k < lobes; k++) {
			float* lobe = &this->model[(v * lobes + k) * VALUES_PER_LOBE];
			const glm::vec3 mean = seeds[k];
			lobe[0] = 1.f / static_cast<float>(lobes);
			lobe[1] = mean.x;
			lobe[2] = mean.y;
			lobe[3] = mean.z;
			lobe[4] = INITIAL_KAPPA;
		}
		this->update_log_norm(v);
	}
}

void VMFTrainer::update_log_norm(size_t voxel)
{
	for (uint32_t k = 0; k < this->lobes; k++) {
		const float* lobe = &this->model[(voxel * this->lobes + k) * VALUES_PER_LOBE];
		this->log_norm[voxel * this->lobes + k] = (lobe[0] > 0.f)
			? static_cast<float>(std::log(static_cast<double>(lobe[0])) + log_vmf_norm(lobe[4]))
			: -INFINITY;
	}
}

void VMFTrainer::add(size_t voxel, const glm::vec3& direction, double weight)
{
	if (!(weight > 0.0))
		return;
	double logp[8];
	double max_logp = -INFINITY;
	const float* model = &this->model[voxel * this->lobes * VALUES_PER_LOBE];
	for (uint32_t k = 0; k < this->lobes; k++) {
		const float* lobe = model + k * VALUES_PER_LOBE;
		const double cosine = lobe[1] * direction.x + lobe[2] * direction.y + lobe[3] * direction.z;
		logp[k] = this->log_norm[voxel * this->lobes + k] + lobe[4] * (cosine - 1.0);
		max_logp = std::max(max_logp, logp[k]);
	}
	if (!std::isfinite(max_logp))
		return;
	double sum = 0.0;
	for (uint32_t k = 0; k < this->lobes; k++) {
		logp[k] = std::exp(logp[k] - max_logp);
		sum += logp[k];
	}
	double* stats = &this->current[voxel * this->lobes * STATS_PER_LOBE];
	for (uint32_t k = 0; k < this->lobes; k++) {
		const double r = weight * logp[k] / sum;
		stats[k * STATS_PER_LOBE + 0] += r;
		stats[k * STATS_PER_LOBE + 1] += r * direction.x;
		stats[k * STATS_PER_LOBE + 2] += r * direction.y;
		stats[k * STATS_PER_LOBE + 3] += r * direction.z;
	}
}

void VMFTrainer::fit(const double* stats, const float* fallback, float* out, bool merge) const
{
	double total = 0.0;
	for (uint32_t k = 0; k < this->lobes; k++)
		total += stats[k * STATS_PER_LOBE];
	if (total <= 0.0) {
		std::copy(fallback, fallback + this->lobes * VALUES_PER_LOBE, out);
		return;
	}
	for (uint32_t k = 0; k < this->lobes; k++) {
		const double* s = stats + k * STATS_PER_LOBE;
		const float* old = fallback + k * VALUES_PER_LOBE;
		float* lobe = out + k * VALUES_PER_LOBE;
		const double n = s[0];
		const double length = std::sqrt(s[1] * s[1] + s[2] * s[2] + s[3] * s[3]);
		lobe[0] = static_cast<float>(n / total);
		if (n / total < MIN_LOBE_SHARE || length <= 0.0) {
			std::copy(old + 1, old + VALUES_PER_LOBE, lobe + 1);
			continue;
		}
		lobe[1] = static_cast<float>(s[1] / length);
		lobe[2] = static_cast<float>(s[2] / length);
		lobe[3] = static_cast<float>(s[3] / length);
		lobe[4] = static_cast<float>(kappa_from_resultant(length / n));
	}
	if (merge)
		merge_overlapping(out, this->lobes, VALUES_PER_LOBE);
}

void VMFTrainer::m_step()
{
	std::vector<float> fitted(this->lobes * VALUES_PER_LOBE);
	for (size_t v = 0; v < this->voxel_count; v++) {
		float* model = &this->model[v * this->lobes * VALUES_PER_LOBE];
		// training keeps all lobes: early lobes are wide and would merge before they separate
		this->fit(&this->current[v * this->lobes * STATS_PER_LOBE], model, fitted.data(), false);
		std::copy(fitted.begin(), fitted.end(), model);
		this->update_log_norm(v);
	}
	this->previous.swap(this->current);
	std::fill(this->current.begin(), this->current.end(), 0.0);
	this->passes++;
}

void VMFTrainer::write_lobes(size_t voxel, float* out) const
{
	double stats[8 * STATS_PER_LOBE];
	double total = 0.0;
	for (size_t i = 0; i < this->lobes * STATS_PER_LOBE; i++) {
		stats[i] = this->current[voxel * this->lobes * STATS_PER_LOBE + i] + this->previous[voxel * this->lobes * STATS_PER_LOBE + i];
		if (i % STATS_PER_LOBE == 0)
			total += stats[i];
	}
	if (total <= 0.0) {
		std::fill(out, out + this->lobes * VALUES_PER_LOBE, 0.f);
		return;
	}
	this->fit(stats, &this->model[voxel * this->lobes * VALUES_PER_LOBE], out, true);

	// canonical slots for datasets: strongest lobe first, unused (weight 0) slots all zero
	float sorted[8 * VALUES_PER_LOBE];
	uint32_t order[8];
	for (uint32_t k = 0; k < this->lobes; k++)
		order[k] = k;
	std::stable_sort(order, order + this->lobes, [out](uint32_t a, uint32_t b) { return out[a * VALUES_PER_LOBE] > out[b * VALUES_PER_LOBE]; });
	for (uint32_t k = 0; k < this->lobes; k++) {
		const float* lobe = out + order[k] * VALUES_PER_LOBE;
		if (lobe[0] > 0.f)
			std::copy(lobe, lobe + VALUES_PER_LOBE, sorted + k * VALUES_PER_LOBE);
		else
			std::fill(sorted + k * VALUES_PER_LOBE, sorted + (k + 1) * VALUES_PER_LOBE, 0.f);
	}
	std::copy(sorted, sorted + this->lobes * VALUES_PER_LOBE, out);
}
