#include "Geometry.hpp"
#include <Eigen/Dense>
#include <glm/gtc/quaternion.hpp>
#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>


using namespace RadiationSimulation::Geometry;


Mesh::Mesh(const std::vector<glm::vec3>& vertices, const std::vector<Face*>& faces, const std::string& name)
	:	vertices(vertices),
		faces(faces),
		name(name)
{
}

Mesh::~Mesh()
{
	for (Face* face : this->faces)
		delete face;
}

void Mesh::setType(const std::string& type)
{
	if (type.empty())
		throw std::invalid_argument("A mesh type must not be empty.");
	if (type.length() > Mesh::MAX_TYPE_LENGTH)
		throw std::invalid_argument("The mesh type \"" + type + "\" is longer than " + std::to_string(Mesh::MAX_TYPE_LENGTH) + " characters.");
	std::string lower_case(type);
	std::transform(lower_case.begin(), lower_case.end(), lower_case.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	this->type = lower_case;
}

void Mesh::attachMaterialName(const std::string name)
{
	this->material_name = name;
}

const std::string& Mesh::getMaterialName() const
{
	return this->material_name;
}

const std::vector<glm::vec3>& Mesh::getVertices() const
{
	return this->vertices;
}

const std::vector<Face*>& Mesh::getFaces() const
{
	return this->faces;
}

const std::string& Mesh::getName() const
{
	return this->name;
}

size_t Mesh::vertexCount() const
{
	return this->vertices.size();
}

size_t Mesh::faceCount() const
{
	return this->faces.size();
}

std::shared_ptr<OrientedBoundingBox> Mesh::getOrientedBoundingBox() const
{
	std::vector<glm::vec3> transformed_vertices(this->vertices.size());
	for (size_t i = 0; i < this->vertices.size(); i++) {
		glm::vec3 vertex = this->vertices[i];
		glm::quat rotation = glm::quat(this->rotation);
		vertex = glm::vec3(glm::mat4_cast(rotation) * glm::vec4(vertex, 1.0f));
		vertex += this->position;
		transformed_vertices[i] = vertex;
	}

	Eigen::MatrixXf points(3, transformed_vertices.size());
	size_t i = 0;
	for (auto& v : transformed_vertices)
		points.col(i++) = Eigen::Vector3f(v.x, v.y, v.z);

	Eigen::Vector3f centroid = points.colwise().mean();
	points.colwise() -= centroid;
	Eigen::Matrix3f cov = points.transpose() * points;
	Eigen::EigenSolver<Eigen::Matrix3f> solver(cov);
	Eigen::Matrix3f eigenvectors = solver.eigenvectors().real();
	Eigen::Vector3f eigenvalues = solver.eigenvalues().real();
	std::vector<int> indices = { 0, 1, 2 };
	std::sort(indices.begin(), indices.end(), [&eigenvalues](int i1, int i2) { return eigenvalues[i1] > eigenvalues[i2]; });

	Eigen::Matrix3f axes(3, 3);
	axes.col(0) = eigenvectors.col(indices[0]);
	axes.col(1) = eigenvectors.col(indices[1]);
	axes.col(2) = eigenvectors.col(indices[2]);

	Eigen::MatrixXf points_in_obb = points * axes;
	Eigen::Vector3f min_coeffs = points_in_obb.colwise().minCoeff();
	Eigen::Vector3f max_coeffs = points_in_obb.colwise().maxCoeff();

	std::array<glm::vec3, 3> axes_list = {
		glm::vec3(axes(0, 0), axes(1, 0), axes(2, 0)),
		glm::vec3(axes(0, 1), axes(1, 1), axes(2, 1)),
		glm::vec3(axes(0, 2), axes(1, 2), axes(2, 2))
	};

	return std::make_shared<OrientedBoundingBox>(glm::vec3(centroid.x(), centroid.y(), centroid.z()), axes_list, glm::vec3(min_coeffs.x(), min_coeffs.y(), min_coeffs.z()), glm::vec3(max_coeffs.x(), max_coeffs.y(), max_coeffs.z()));
}

std::shared_ptr<Collisions::Capsule> Mesh::getBoundingCapsule() const
{
	auto obb = this->getOrientedBoundingBox();

	glm::vec3 p1 = obb->centroid + obb->axes[0] * obb->min_coeffs.x;
	glm::vec3 p2 = obb->centroid + obb->axes[0] * obb->max_coeffs.x;

	// The radius is half the length of the second longest axis of the OBB
	float radius = std::max(
		glm::length(obb->axes[1] * obb->max_coeffs.y - obb->axes[1] * obb->min_coeffs.y),
		glm::length(obb->axes[2] * obb->max_coeffs.z - obb->axes[2] * obb->min_coeffs.z)
	) / 2.0f;

	return std::make_shared<Collisions::Capsule>(p1, p2, radius);
}

OrientedBoundingBox::OrientedBoundingBox(const glm::vec3& centroid, const std::array<glm::vec3, 3>& axes, const glm::vec3& min_coeffs, const glm::vec3& max_coeffs)
	: centroid(centroid),
	  axes(axes),
	  min_coeffs(min_coeffs),
	  max_coeffs(max_coeffs)
{
}
