#include "gtest/gtest.h"
#include "RadiationSimulation.hpp"
#include "World.hpp"
#include "Geometry.hpp"
#include <Geant4/G4RadiationFieldDetector.hpp>
#include <radfiled3d/radiation_field.hpp>
#include <G4Material.hh>
#include <G4SystemOfUnits.hh>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Short simulation of a patient (water with a bone child) and a concrete block side by side in the beam: nothing may be
// scored inside the patient's materials, while the concrete is scored and unscattered primaries crossing it stay in
// direct_beam. Needs the Geant4 data environment (G4LEDATA etc.).

using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;

namespace {
	std::shared_ptr<Mesh> make_box(const std::string& name, const glm::vec3& half)
	{
		const std::vector<glm::vec3> vertices = {
			{ -half.x, -half.y, -half.z }, { half.x, -half.y, -half.z }, { half.x, half.y, -half.z }, { -half.x, half.y, -half.z },
			{ -half.x, -half.y, half.z }, { half.x, -half.y, half.z }, { half.x, half.y, half.z }, { -half.x, half.y, half.z }
		};
		const std::vector<Face*> faces = {
			new TriFace(glm::uvec3(0, 3, 2)), new TriFace(glm::uvec3(0, 2, 1)),
			new TriFace(glm::uvec3(4, 5, 6)), new TriFace(glm::uvec3(4, 6, 7)),
			new TriFace(glm::uvec3(0, 1, 5)), new TriFace(glm::uvec3(0, 5, 4)),
			new TriFace(glm::uvec3(2, 3, 7)), new TriFace(glm::uvec3(2, 7, 6)),
			new TriFace(glm::uvec3(1, 2, 6)), new TriFace(glm::uvec3(1, 6, 5)),
			new TriFace(glm::uvec3(0, 4, 7)), new TriFace(glm::uvec3(0, 7, 3))
		};
		return std::make_shared<Mesh>(vertices, faces, name);
	}

	struct Box {
		glm::vec3 centre, half;
		bool contains(const glm::vec3& low, float size) const {
			return glm::all(glm::greaterThanEqual(low, this->centre - this->half)) && glm::all(glm::lessThanEqual(low + size, this->centre + this->half));
		}
	};
}

TEST(PatientScoring, NothingIsScoredInsideThePatientsMaterials) {
	RadiationSimulator::initialize(2);
	const glm::vec3 world_dim(1.f);
	const float voxel_dim = 0.05f;
	RadiationSimulator::set_world_info(std::make_unique<WorldInfo>("Air", world_dim));

	// no box face lies on a voxel face: steps starting on a surface touch the voxels on both sides of it
	const Box patient_box{ glm::vec3(-0.1f, 0.f, 0.f), glm::vec3(0.11f, 0.16f, 0.11f) };
	const Box bone_box{ patient_box.centre, glm::vec3(0.06f) };
	const Box concrete_box{ glm::vec3(0.13f, 0.f, 0.f), glm::vec3(0.09f, 0.16f, 0.11f) };
	auto patient = make_box("patient", patient_box.half);
	patient->setType(Mesh::PATIENT_TYPE);
	patient->attachMaterialName("G4_WATER");
	patient->setPosition(patient_box.centre);
	auto bone = make_box("bone", bone_box.half);   // no Type: inherits the patient's
	bone->attachMaterialName("G4_BONE_COMPACT_ICRU");
	patient->addChild(bone);
	auto concrete = make_box("concrete", concrete_box.half);
	concrete->setType("Room");
	concrete->attachMaterialName("G4_CONCRETE");
	concrete->setPosition(concrete_box.centre);
	RadiationSimulator::add_geometry(std::vector<std::shared_ptr<Mesh>>{ patient, concrete });

	auto source = std::make_shared<XRaySource>(50e3f, std::make_unique<RectangleSourceShape>(glm::vec2(0.5f), 1.f));
	source->setTransform(glm::vec3(0.f, 0.f, 1.f), glm::vec3(0.f, 0.f, -1.f));
	RadiationSimulator::add_radiation_source(source);
	RadiationSimulator::set_radiation_field_resolution(world_dim, glm::vec3(voxel_dim), 60e3f * eV, 1e3f * eV, 0.f, 1.f);

	auto field = std::static_pointer_cast<radfiled3d::CartesianRadiationField>(RadiationSimulator::simulate_radiation_field(20000));

	std::set<std::string> patient_materials;
	for (const G4Material* material : World::Get()->get_radiation_field_detector()->get_patient_materials())
		patient_materials.insert(material->GetName());
	EXPECT_EQ(patient_materials, (std::set<std::string>{ "G4_WATER", "G4_BONE_COMPACT_ICRU" }));

	auto scatter = field->get_channel("scatter_field");
	auto direct = field->get_channel("direct_beam");
	const glm::uvec3 counts = scatter->get_voxel_counts();
	size_t in_patient = 0, in_bone = 0, in_concrete = 0;
	double concrete_scatter = 0.0, concrete_direct = 0.0;
	for (unsigned z = 0; z < counts.z; z++) {
		for (unsigned y = 0; y < counts.y; y++) {
			for (unsigned x = 0; x < counts.x; x++) {
				const glm::vec3 low = glm::vec3(x, y, z) * voxel_dim - glm::vec3(counts) * voxel_dim / 2.f;
				const float s = scatter->get_voxel<radfiled3d::ScalarVoxel<float>>("flux", x, y, z).get_data();
				const float d = direct->get_voxel<radfiled3d::ScalarVoxel<float>>("flux", x, y, z).get_data();
				if (patient_box.contains(low, voxel_dim)) {
					in_patient++;
					in_bone += bone_box.contains(low, voxel_dim);
					EXPECT_EQ(s + d, 0.f) << "voxel (" << x << ", " << y << ", " << z << ") lies inside the patient but has flux";
				}
				if (concrete_box.contains(low, voxel_dim)) {
					in_concrete++;
					concrete_scatter += s;
					concrete_direct += d;
				}
			}
		}
	}
	ASSERT_GT(in_patient, 0u);
	ASSERT_GT(in_bone, 0u);
	ASSERT_GT(in_concrete, 0u);
	EXPECT_GT(concrete_scatter, 0.0);
	EXPECT_GT(concrete_direct, 0.0);

	RadiationSimulator::deinitialize();
}
