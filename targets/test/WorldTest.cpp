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
	std::shared_ptr<RadiationSource> source = std::make_shared<XRaySource>(10.0e+6);
	source->setTransform(glm::vec3(-1.f, 0.f, 0.f), glm::vec3(1.f, 0.f, 0.f));
	RadiationSimulator::add_radiation_source(source);

	RadiationSimulator::set_radiation_field_resolution(glm::vec3(1.f), glm::vec3(0.1f), 150.f, 4.f);

	auto storage = std::make_shared<RadFiled3D::FieldStoreV1>();
	RadiationSimulator::set_exporter(storage);
	RadiationSimulator::set_importer(storage);

	size_t count = 0;

	RadFiled3D::RadiationFieldMetadata metadata;
	while (count < 100) {
		std::shared_ptr<RadFiled3D::IRadiationField> test_field = RadiationSimulator::simulate_radiation_field(0.1f)->get_future().get();
		RadiationSimulator::export_field(test_field, metadata, "./test.rf");
		auto test_field_loaded = RadiationSimulator::import_field("./test.rf");

		count++;
	}

	std::cout << "Test passed" << std::endl;

	RadiationSimulator::deinitialize();
	return 0;
}