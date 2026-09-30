#include "GeometryLoader.hpp"
#include "Geometry.hpp"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>


using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;

namespace {
	// Writes an OBJ with one tetrahedron per object name and the given description next to it.
	std::string write_scene(const std::string& basename, const std::vector<std::string>& objects, const std::string& description)
	{
		const std::string obj_path = (std::filesystem::temp_directory_path() / (basename + ".obj")).string();
		std::ofstream obj(obj_path);
		for (size_t i = 0; i < objects.size(); i++) {
			const float offset = static_cast<float>(i);
			obj << "o " << objects[i] << "\n";
			obj << "v " << offset << " 0 0\nv " << offset + 0.1f << " 0 0\nv " << offset << " 0.1 0\nv " << offset << " 0 0.1\n";
			const size_t base = i * 4 + 1;
			obj << "f " << base << " " << base + 2 << " " << base + 1 << "\n";
			obj << "f " << base << " " << base + 1 << " " << base + 3 << "\n";
			obj << "f " << base << " " << base + 3 << " " << base + 2 << "\n";
			obj << "f " << base + 1 << " " << base + 2 << " " << base + 3 << "\n";
		}
		std::ofstream(std::filesystem::temp_directory_path() / (basename + ".desc")) << description;
		return obj_path;
	}

	void collect_types(const std::shared_ptr<Mesh>& mesh, std::map<std::string, std::string>& types)
	{
		types[mesh->getName()] = mesh->getType();
		for (const auto& child : mesh->getChildren())
			collect_types(child, types);
	}

	std::map<std::string, std::string> load_types(const std::string& obj_path)
	{
		std::map<std::string, std::string> types;
		for (const auto& mesh : GeometryLoader::Load(obj_path))
			collect_types(mesh, types);
		return types;
	}
}

TEST(GeometryLoader, TypesAreCaseInsensitiveAndInherited) {
	const std::string obj = write_scene("rf3_types", { "body", "lung", "bronchus", "implant", "table", "plate" }, R"({
		"body": { "Type": "Patient", "Children": {
			"lung": { "Children": { "bronchus": {} } },
			"implant": { "Type": "METAL" }
		} },
		"table": { "Type": "Shield", "Children": { "plate": { "Type": "sHiElD" } } }
	})");
	const auto types = load_types(obj);
	EXPECT_EQ(types.at("body"), "patient");
	EXPECT_EQ(types.at("lung"), "patient");
	EXPECT_EQ(types.at("bronchus"), "patient");
	EXPECT_EQ(types.at("implant"), "metal");
	EXPECT_EQ(types.at("table"), "shield");
	EXPECT_EQ(types.at("plate"), "shield");
}

TEST(GeometryLoader, UntypedRootsAreUnknownAndLegacyPatientFlagIsRead) {
	const std::string obj = write_scene("rf3_legacy", { "patient", "lung", "stand" }, R"({
		"patient": { "Patient": true, "Children": { "lung": { "Patient": false } } },
		"stand": {}
	})");
	const auto types = load_types(obj);
	EXPECT_EQ(types.at("patient"), Mesh::PATIENT_TYPE);
	EXPECT_EQ(types.at("lung"), Mesh::PATIENT_TYPE);
	EXPECT_EQ(types.at("stand"), Mesh::DEFAULT_TYPE);
}

TEST(GeometryLoader, OnlyOnePatientMayBeDeclared) {
	const std::string obj = write_scene("rf3_two_patients", { "a", "b" }, R"({
		"a": { "Type": "patient" },
		"b": { "Type": "PATIENT" }
	})");
	EXPECT_THROW(GeometryLoader::Load(obj), std::runtime_error);
}

TEST(GeometryLoader, TypeMustFitALayerName) {
	const std::string obj = write_scene("rf3_long_type", { "a" }, "{ \"a\": { \"Type\": \"" + std::string(Mesh::MAX_TYPE_LENGTH + 1, 'x') + "\" } }");
	EXPECT_THROW(GeometryLoader::Load(obj), std::runtime_error);

	const std::string empty = write_scene("rf3_empty_type", { "a" }, R"({ "a": { "Type": "" } })");
	EXPECT_THROW(GeometryLoader::Load(empty), std::runtime_error);
}

TEST(GeometryLoader, IsocenterDistanceOnlyOnRootImageDetectors) {
	const std::string obj = write_scene("rf3_detector_range", { "plate", "housing" }, R"({
		"plate": { "Type": "ImageDetector", "IsocenterDistance": { "Min": 0.45, "Max": 0.62 } },
		"housing": { "Type": "XRayTube" }
	})");
	const auto meshes = GeometryLoader::Load(obj);
	for (const auto& mesh : meshes) {
		if (mesh->getName() == "plate") {
			ASSERT_TRUE(mesh->getIsocenterDistanceRange().has_value());
			EXPECT_FLOAT_EQ(mesh->getIsocenterDistanceRange()->x, 0.45f);
			EXPECT_FLOAT_EQ(mesh->getIsocenterDistanceRange()->y, 0.62f);
		}
		else {
			EXPECT_FALSE(mesh->getIsocenterDistanceRange().has_value());
		}
	}

	const std::string on_tube = write_scene("rf3_tube_range", { "a" }, R"({ "a": { "Type": "XRayTube", "IsocenterDistance": { "Min": 0.4, "Max": 0.6 } } })");
	EXPECT_THROW(GeometryLoader::Load(on_tube), std::runtime_error);
	const std::string on_child = write_scene("rf3_child_range", { "a", "b" }, R"({ "a": { "Type": "ImageDetector", "Children": { "b": { "IsocenterDistance": { "Min": 0.4, "Max": 0.6 } } } } })");
	EXPECT_THROW(GeometryLoader::Load(on_child), std::runtime_error);
	const std::string reversed = write_scene("rf3_reversed_range", { "a" }, R"({ "a": { "Type": "ImageDetector", "IsocenterDistance": { "Min": 0.6, "Max": 0.4 } } })");
	EXPECT_THROW(GeometryLoader::Load(reversed), std::runtime_error);
}
