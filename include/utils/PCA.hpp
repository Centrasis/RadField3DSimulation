#pragma once
#include <glm/glm.hpp>
#include <Eigen/Dense>


class OnlinePCA {
public:
    OnlinePCA();

    void addVector(const glm::vec3& vector);

    glm::vec3 getPrincipalDirection() const;

    void reset();

private:
    size_t count;
    glm::vec3 mean;
    Eigen::Matrix3f covariance;
};
