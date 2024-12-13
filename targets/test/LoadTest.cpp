#include "RadiationSimulation.hpp"
#include "GeometryLoader.hpp"
#include <stdio.h>
#include <iostream>
#include <G4ios.hh>
#include <chrono>
#include "RadiationFieldStore.hpp"
#include <stdexcept>
#include "World.hpp"


using namespace RadiationSimulation;

int main() {
	RadiationSimulator::initialize(RadiationHandlerType::Geant4);

	RadiationSimulator::set_world_info(std::make_shared<WorldInfo>("Air", glm::vec3(1.f)));

	auto meshes = GeometryLoader::Load("../data/monkey.obj");
	G4cout << "Meshes loaded: " << meshes.size() << G4endl;

	G4cout << "Add meshes...." << meshes.size() << G4endl;
	RadiationSimulator::add_geometry(meshes);
	G4cout << "Add radiation...." << meshes.size() << G4endl;
	auto source = std::make_shared<XRaySource>(10.0e+6);
	source->setTransform(glm::vec3(-1.f, 0.f, 0.f), glm::vec3(1.f, 0.f, 0.f));
	RadiationSimulator::add_radiation_source(source);

	RadiationSimulator::set_radiation_field_resolution(glm::vec3(1.f), glm::vec3(0.1f));

	auto storage = std::make_shared<BinaryRadiationFieldStore>();
	RadiationSimulator::set_exporter(storage);
	RadiationSimulator::set_importer(storage);

	size_t count = 0;

	while (count < 100) {
		std::shared_ptr<DirectedSourceRadiationField> test_field = RadiationSimulator::simulate_radiation_field(0.1f)->get_future().get();
		double overall_hits = test_field->get_overall_hits();
		std::cout << "Overall hits: " << overall_hits << std::endl;
		RadiationSimulator::export_field(test_field, "./test.rf");
		auto test_field_loaded = RadiationSimulator::import_field("./test.rf");

		if (overall_hits <= 0) {
			std::cout << "Error: Hits = " << overall_hits << " <= 0" << std::endl;
			return 1;
		}

		if (test_field_loaded->get_overall_hits() != overall_hits) {
			std::cout << "Error: " << test_field_loaded->get_overall_hits() << " (loaded) != " << overall_hits << " (stored)" << std::endl;
			return 1;
		}

		auto steps = test_field_loaded->getSteps();

		// test pca
		bool pca_valid = true;
		std::cout << "Test all PCA results..." << std::endl;
		for (size_t x = 0; x < steps.x; x++)
			for (size_t y = 0; y < steps.y; y++)
				for (size_t z = 0; z < steps.z; z++)
				{
					auto& vx1 = test_field->get_voxel(x, y, z);
					auto vec1 = vx1.getPrincipalDirection();
					float len1 = vec1.x * vec1.x + vec1.y * vec1.y + vec1.z * vec1.z;

					auto& vx2 = test_field_loaded->get_voxel(x, y, z);
					auto vec2 = vx1.getPrincipalDirection();
					float len2 = vec2.x * vec2.x + vec2.y * vec2.y + vec2.z * vec2.z;

					if (len1 != len2) {
						std::cout << "Error: Simulated and Loaded principal direction of voxel: (" << x << ", " << y << ", " << z << ") was mismatched (" << len1 << " != " << len2 << ")!" << std::endl;
						std::cout << "Voxel properties:" << std::endl;
						std::cout << "Doserate:\t" << vx1.getDoserate() << std::endl;
						std::cout << "Energy:\t\t" << vx1.getEnergy() << std::endl;
						std::cout << "Hits:\t\t" << vx1.getHits() << std::endl;
						std::cout << "Bin Width:\t" << vx1.getSpectrumBinWidth() << std::endl << std::endl;
						pca_valid = false;
					}

					if (len1 == 0.0 || len1 != len1)
						continue;

					len1 = std::sqrt(len1);
					if (std::abs(len1 - 1.0) > 0.001) {
						std::cout << "Error: Simulated principal direction of voxel: (" << x << ", " << y << ", " << z << ") was invalid with length (" << len1 << ") != 1.0!" << std::endl << std::endl;
						pca_valid = false;
					}
				}
		if (!pca_valid)
		{
			std::cout << "PCA failed!" << std::endl;
			return -1;
		}
		std::cout << "PCA PASSED!" << std::endl << std::endl << std::endl;

		std::cout << "Test all PCA results of the loaded field..." << std::endl;
		for (size_t x = 0; x < steps.x; x++)
			for (size_t y = 0; y < steps.y; y++)
				for (size_t z = 0; z < steps.z; z++)
				{
					auto& vx = test_field_loaded->get_voxel(x, y, z);
					auto vec = vx.getPrincipalDirection();
					float len = vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
					if (len == 0.0)
						continue;
					len = std::sqrt(len);
					if (std::abs(len - 1.0) > 0.01) {
						std::cout << "Error: Loaded principal direction of voxel: (" << x << ", " << y << ", " << z << ") was invalid with length (" << len << ") != 1.0!" << std::endl;
						return 1;
					}
				}
		std::cout << "PCA PASSED!" << std::endl << std::endl << std::endl;

		count++;
	}

	std::cout << "Test passed" << std::endl;

	RadiationSimulator::deinitialize();
	return 0;
}