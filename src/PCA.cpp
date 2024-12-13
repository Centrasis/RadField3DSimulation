#include <utils/PCA.hpp>

OnlinePCA::OnlinePCA()
    : covariance(Eigen::Matrix3f::Zero())
{
    this->reset();
}

void OnlinePCA::addVector(const glm::vec3& vector)
{
    glm::vec3 normalized = glm::normalize(vector);
    
    // Ensure the vector is normalizable
    if (glm::length(normalized) == 0.0f) {
        return;
    }

    Eigen::Vector3f vec(normalized.x, normalized.y, normalized.z);

    if (count == 0) {
        this->mean = normalized;
        this->count = 1;
    }
    else {
        glm::vec3 delta = normalized - this->mean;
        this->mean += delta / static_cast<float>(++this->count);
        glm::vec3 deltaMean = normalized - this->mean;
        covariance += Eigen::Vector3f(delta.x, delta.y, delta.z) * Eigen::Vector3f(deltaMean.x, deltaMean.y, deltaMean.z).transpose();
    }
}

glm::vec3 OnlinePCA::getPrincipalDirection() const
{
    if (count < 2) {
        return this->mean;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver(covariance / static_cast<float>(this->count - 1));
    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to compute eigenvalues and eigenvectors");
    }

    Eigen::Vector3f::Index maxIndex;
    eigensolver.eigenvalues().maxCoeff(&maxIndex);
    Eigen::Vector3f principalDirection = eigensolver.eigenvectors().col(maxIndex).normalized();

    return glm::vec3(principalDirection(0), principalDirection(1), principalDirection(2));
}

void OnlinePCA::reset()
{
	this->count = 0;
	this->mean = glm::vec3(0.f);
	this->covariance = Eigen::Matrix3f::Zero();
}
