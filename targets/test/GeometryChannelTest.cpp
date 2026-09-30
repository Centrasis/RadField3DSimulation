#include "gtest/gtest.h"
#include "RadiationSimulation.hpp"
#include "World.hpp"
#include "Geometry.hpp"
#include <Geant4/G4Geometry.hpp>
#include <Geant4/G4RadiationFieldDetector.hpp>
#include <radfiled3d/radiation_field.hpp>
#include <G4Navigator.hh>
#include <G4TransportationManager.hh>
#include <G4TouchableHistory.hh>
#include <G4SystemOfUnits.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <Geant4/G4World.hpp>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

// Runs a short simulation of a nested, rotated and translated box scene and checks the voxelized geometry channel
// against Geant4's own navigation of the scene and against the exact box/voxel overlap. Needs the Geant4 data environment (G4LEDATA etc.).

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

	struct PlacedBox {
		std::string type;
		glm::dmat3 rotation;
		glm::dvec3 center;
		glm::dvec3 half;
	};

	// Collects the boxes as Geant4 placed them below `volume` (world frame, metres).
	void collect_boxes(const G4LogicalVolume& volume, const G4RotationMatrix& rotation, const G4ThreeVector& translation, const std::map<std::string, glm::vec3>& half_extents, std::vector<PlacedBox>& boxes)
	{
		for (size_t i = 0; i < volume.GetNoDaughters(); i++) {
			const G4VPhysicalVolume* daughter = volume.GetDaughter(i);
			const G4RotationMatrix r = rotation * daughter->GetObjectRotationValue();
			const G4ThreeVector t = rotation * daughter->GetObjectTranslation() + translation;
			if (const auto* mesh = dynamic_cast<const Geant4::Mesh*>(daughter->GetLogicalVolume()->GetSolid())) {
				glm::dmat3 m3;
				for (int c = 0; c < 3; c++)
					for (int row = 0; row < 3; row++)
						m3[c][row] = r(row, c);
				boxes.push_back({ mesh->getMesh()->getType(), m3, glm::dvec3(t.x(), t.y(), t.z()) / m, glm::dvec3(half_extents.at(mesh->getMesh()->getName())) });
			}
			collect_boxes(*daughter->GetLogicalVolume(), r, t, half_extents, boxes);
		}
	}

	// Exact separating axis test: does the box overlap the voxel with non-zero volume (touching does not count)?
	bool box_overlaps_voxel(const PlacedBox& box, const glm::dvec3& voxel_center, double voxel_half)
	{
		const glm::dvec3 a[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
		const glm::dvec3 b[3] = { box.rotation[0], box.rotation[1], box.rotation[2] };
		const glm::dvec3 offset = box.center - voxel_center;
		auto separated = [&](const glm::dvec3& axis) {
			if (glm::dot(axis, axis) < 1e-12)
				return false;
			double r_voxel = 0.0, r_box = 0.0;
			for (int i = 0; i < 3; i++) {
				r_voxel += voxel_half * std::abs(glm::dot(a[i], axis));
				r_box += box.half[i] * std::abs(glm::dot(b[i], axis));
			}
			// the placements come from float geometry: faces meant to lie on a voxel plane may be off by ~1e-8 m
			return std::abs(glm::dot(offset, axis)) >= r_voxel + r_box - 1e-6 * glm::length(axis);
		};
		for (int i = 0; i < 3; i++) {
			if (separated(a[i]) || separated(b[i]))
				return false;
			for (int j = 0; j < 3; j++)
				if (separated(glm::cross(a[i], b[j])))
					return false;
		}
		return true;
	}
}

TEST(GeometryChannel, MatchesGeant4Navigation) {
	RadiationSimulator::initialize(2);
	const glm::vec3 world_dim(1.f);
	const float voxel_dim = 0.05f;
	RadiationSimulator::set_world_info(std::make_unique<WorldInfo>("Air", world_dim));

	const glm::vec3 patient_half(0.25f, 0.15f, 0.2f), organ_half(0.08f, 0.06f, 0.1f), shield_half(0.1f, 0.2f, 0.01f), second_shield_half(0.02f, 0.1f, 0.1f);
	auto patient = make_box("patient", patient_half);
	patient->setType(Mesh::PATIENT_TYPE);
	patient->attachMaterialName("G4_WATER");
	patient->setRotation(glm::vec3(0.3f, 0.5f, 0.2f));
	patient->setPosition(glm::vec3(0.1f, -0.05f, 0.08f));
	auto organ = make_box("organ", organ_half);
	organ->setType("Organ");
	organ->attachMaterialName("G4_LUNG_ICRP");
	organ->setRotation(glm::vec3(0.f, 0.f, 0.6f));
	organ->setPosition(glm::vec3(0.1f, 0.02f, 0.f));
	patient->addChild(organ);
	auto shield = make_box("shield", shield_half);
	shield->setType("Shield");
	shield->attachMaterialName("G4_Pb");
	shield->setPosition(glm::vec3(-0.35f, 0.f, 0.f));
	auto second_shield = make_box("second_shield", second_shield_half);
	second_shield->setType("SHIELD");
	second_shield->attachMaterialName("G4_Pb");
	second_shield->setRotation(glm::vec3(0.f, 0.4f, 0.f));
	second_shield->setPosition(glm::vec3(0.1f, 0.35f, -0.3f));
	RadiationSimulator::add_geometry(std::vector<std::shared_ptr<Mesh>>{ patient, shield, second_shield });

	auto source = std::make_shared<XRaySource>(50e3f, std::make_unique<RectangleSourceShape>(glm::vec2(0.2f), 1.f));
	source->setTransform(glm::vec3(0.f, 0.f, 1.f), glm::vec3(0.f, 0.f, -1.f));
	RadiationSimulator::add_radiation_source(source);
	RadiationSimulator::set_radiation_field_resolution(world_dim, glm::vec3(voxel_dim), 60e3f * eV, 1e3f * eV, 0.f, 1.f);

	auto field = std::static_pointer_cast<radfiled3d::CartesianRadiationField>(RadiationSimulator::simulate_radiation_field(1000));
	World::Get()->get_radiation_field_detector()->add_geometry_channel_to(*field);
	ASSERT_TRUE(field->has_channel("geometry"));
	auto geometry = field->get_channel("geometry");
	ASSERT_TRUE(geometry->has_layer("patient"));
	ASSERT_TRUE(geometry->has_layer("organ"));
	EXPECT_EQ(geometry->get_layers().size(), 3u);   // patient, organ and one merged shield layer

	G4Navigator navigator;
	navigator.SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());

	const glm::uvec3 counts = geometry->get_voxel_counts();
	const glm::vec3 half_field = glm::vec3(counts) * voxel_dim / 2.f;
	// types of all mesh volumes (the located one and its ancestors) containing the point
	auto types_at = [&](const glm::vec3& point) {
		std::set<std::string> types;
		navigator.LocateGlobalPointAndSetup(G4ThreeVector(point.x * m, point.y * m, point.z * m), nullptr, false, true);
		std::unique_ptr<G4TouchableHistory> touchable(navigator.CreateTouchableHistory());
		for (G4int depth = 0; depth <= touchable->GetHistoryDepth(); depth++) {
			const auto* mesh = dynamic_cast<const Geant4::Mesh*>(touchable->GetVolume(depth)->GetLogicalVolume()->GetSolid());
			if (mesh != nullptr)
				types.insert(mesh->getMesh()->getType());
		}
		return types;
	};

	const std::map<std::string, glm::vec3> half_extents = {
		{ "patient", patient_half }, { "organ", organ_half }, { "shield", shield_half }, { "second_shield", second_shield_half }
	};
	std::vector<PlacedBox> boxes;
	collect_boxes(*Geant4::World::Get()->get_volume(), G4RotationMatrix(), G4ThreeVector(), half_extents, boxes);
	ASSERT_EQ(boxes.size(), 4u);

	std::map<std::string, size_t> centre_inside, overlapping, marked, missed;
	for (unsigned z = 0; z < counts.z; z++) {
		for (unsigned y = 0; y < counts.y; y++) {
			for (unsigned x = 0; x < counts.x; x++) {
				const size_t idx = geometry->get_voxel_idx(x, y, z);
				const glm::vec3 centre = (glm::vec3(x, y, z) + 0.5f) * voxel_dim - half_field;
				// independent check of the placements: Geant4's navigator
				for (const std::string& type : types_at(centre)) {
					centre_inside[type]++;
					EXPECT_EQ(geometry->get_layer<uint8_t>(type)[idx], 255) << type << " voxel (" << x << ", " << y << ", " << z << ") has its centre inside the geometry but is not marked";
				}
				for (const std::string& type : geometry->get_layers()) {
					bool overlaps = false;
					for (const PlacedBox& box : boxes)
						overlaps |= box.type == type && box_overlaps_voxel(box, glm::dvec3(centre), 0.5 * voxel_dim);
					const bool is_marked = geometry->get_layer<uint8_t>(type)[idx] == 255;
					overlapping[type] += overlaps;
					marked[type] += is_marked;
					missed[type] += overlaps && !is_marked;
					EXPECT_FALSE(is_marked && !overlaps) << type << " voxel (" << x << ", " << y << ", " << z << ") is marked but does not overlap the geometry";
				}
			}
		}
	}
	ASSERT_GT(centre_inside["patient"], 0u);
	ASSERT_GT(centre_inside["organ"], 0u);
	for (const auto& [type, count] : overlapping) {
		ASSERT_GT(count, 0u) << type;
		// voxels only cut by a polygon that is not the closest one to their centre may be missed
		EXPECT_LE(missed[type], count / 100) << type;
		std::cout << type << ": " << centre_inside[type] << " voxel centres inside, " << count << " voxels overlapping, " << marked[type] << " marked, " << missed[type] << " missed" << std::endl;
	}

	// The patient is its material and the materials of its children (the organ), whatever Type they declare.
	std::set<std::string> patient_materials;
	for (const G4Material* material : World::Get()->get_radiation_field_detector()->get_patient_materials())
		patient_materials.insert(material->GetName());
	EXPECT_EQ(patient_materials, (std::set<std::string>{ "G4_WATER", "G4_LUNG_ICRP" }));

	// Nothing is scored inside the patient's materials: voxels lying entirely in the patient or the organ have no flux in
	// any channel (PatientScoringTest checks that other geometry is scored).
	auto innermost_type = [&](const glm::vec3& point) -> std::string {
		navigator.LocateGlobalPointAndSetup(G4ThreeVector(point.x * m, point.y * m, point.z * m), nullptr, false, true);
		std::unique_ptr<G4TouchableHistory> touchable(navigator.CreateTouchableHistory());
		const auto* mesh = dynamic_cast<const Geant4::Mesh*>(touchable->GetVolume(0)->GetLogicalVolume()->GetSolid());
		return mesh != nullptr ? mesh->getMesh()->getType() : "";
	};
	size_t inside_patient = 0;
	for (unsigned z = 0; z < counts.z; z++) {
		for (unsigned y = 0; y < counts.y; y++) {
			for (unsigned x = 0; x < counts.x; x++) {
				const glm::vec3 low = glm::vec3(x, y, z) * voxel_dim - half_field;
				bool all_patient = true;
				for (int s = 0; s < 27 && all_patient; s++) {
					const std::string type = innermost_type(low + glm::vec3(s % 3, (s / 3) % 3, s / 9) * (voxel_dim / 2.f));
					all_patient = type == Mesh::PATIENT_TYPE || type == "organ";
				}
				if (!all_patient)
					continue;
				inside_patient++;
				for (const char* channel : { "scatter_field", "direct_beam" })
					EXPECT_EQ(field->get_channel(channel)->get_voxel<radfiled3d::ScalarVoxel<float>>("flux", x, y, z).get_data(), 0.f)
						<< channel << " voxel (" << x << ", " << y << ", " << z << ") lies inside the patient but has flux";
			}
		}
	}
	ASSERT_GT(inside_patient, 0u);

	RadiationSimulator::deinitialize();
}
