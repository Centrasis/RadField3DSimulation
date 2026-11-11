#include "gtest/gtest.h"
#include "RadiationSimulation.hpp"
#include "World.hpp"
#include "Geometry.hpp"
#include "CArm.hpp"
#include "RadiationSource.hpp"
#include "Geant4/G4Geometry.hpp"
#include <G4PhysicalVolumeStore.hh>
#include <G4VPhysicalVolume.hh>
#include <G4Navigator.hh>
#include <G4TransportationManager.hh>
#include <G4SystemOfUnits.hh>
#include <glm/gtc/quaternion.hpp>
#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

// A short simulation with an image detector and a tube housing, tube turned away from its base pose: both must be
// placed where the C-arm puts them. Needs the Geant4 data environment (G4LEDATA etc.).

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

	glm::vec3 to_m(const G4ThreeVector& v) { return glm::vec3(v.x(), v.y(), v.z()) / static_cast<float>(m); }
}

TEST(CArmPlacement, DetectorAndTubeTurnWithTheBeam) {
	RadiationSimulator::initialize(2);
	RadiationSimulator::set_world_info(std::make_unique<WorldInfo>("Air", glm::vec3(2.f)));

	auto patient = make_box("patient", glm::vec3(0.1f, 0.1f, 0.3f));
	patient->setType(Mesh::PATIENT_TYPE);
	patient->attachMaterialName("G4_WATER");
	// detector: a plate 0.5 m above the isocentre in the base pose, deliberately 5 cm off the beam axis
	auto detector = make_box("detector", glm::vec3(0.1f, 0.001f, 0.1f));
	detector->setType("ImageDetector");
	detector->attachMaterialName("CsI");
	const glm::vec3 detector_base(0.05f, 0.5f, 0.f);
	detector->setPosition(detector_base);
	// tube housing: modelled around the focal spot at the origin, behind it as seen along the beam
	auto tube = make_box("tube", glm::vec3(0.05f, 0.03f, 0.05f));
	tube->setType("XRayTube");
	tube->attachMaterialName("G4_Pb");
	const glm::vec3 tube_base(0.f, -0.05f, 0.f);
	tube->setPosition(tube_base);
	RadiationSimulator::add_geometry(std::vector<std::shared_ptr<Mesh>>{ patient, detector, tube });

	// tube turned like RadField3D does for --source-phi 30 --source-theta 60
	const glm::quat angles = glm::angleAxis(glm::radians(60.f), glm::vec3(1.f, 0.f, 0.f)) * glm::angleAxis(glm::radians(30.f), glm::vec3(0.f, 1.f, 0.f));
	const glm::vec3 dir = angles * glm::vec3(0.f, 0.f, -1.f);
	auto source = std::make_shared<XRaySource>(50e3f, std::make_unique<RectangleSourceShape>(glm::vec2(0.1f), 1.f));
	source->setTransform(-dir * 0.6f, dir);
	RadiationSimulator::add_radiation_source(source);
	RadiationSimulator::set_radiation_field_resolution(glm::vec3(1.6f), glm::vec3(0.05f), 60e3f * eV, 1e3f * eV, 0.f, 1.f);
	RadiationSimulator::simulate_radiation_field(1000);

	const glm::quat c_arm = source->getCArmRotation();
	auto placed = [](const std::string& name) {
		G4VPhysicalVolume* pv = G4PhysicalVolumeStore::GetInstance()->GetVolume(name, false);
		EXPECT_NE(pv, nullptr) << name;
		return pv;
	};

	// detector: turned about the isocentre, so it keeps its distance, lies on the beam's far side and faces the beam
	const G4VPhysicalVolume* det = placed("detector");
	const glm::vec3 det_centre = to_m(det->GetObjectTranslation());
	const glm::vec3 expected_centre = c_arm * detector_base;
	EXPECT_LT(glm::length(det_centre - expected_centre), 1e-5f);
	EXPECT_NEAR(glm::length(det_centre), glm::length(detector_base), 1e-5f);
	EXPECT_GT(glm::dot(det_centre, dir), 0.49f);                                   // beyond the isocentre, along the beam
	const G4ThreeVector normal = det->GetObjectRotationValue() * G4ThreeVector(0., 1., 0.);
	EXPECT_NEAR(glm::dot(glm::vec3(normal.x(), normal.y(), normal.z()), dir), 1.f, 1e-5f);   // plane orthogonal to the beam
	// the 5 cm offset turned with it: it is perpendicular to the beam, 5 cm long
	const glm::vec3 offset = det_centre - dir * glm::dot(det_centre, dir);
	EXPECT_NEAR(glm::length(offset), 0.05f, 1e-5f);

	// tube: turned like the beam and moved to the focal spot
	const G4VPhysicalVolume* tb = placed("tube");
	EXPECT_LT(glm::length(to_m(tb->GetObjectTranslation()) - (source->getLocation() + c_arm * tube_base)), 1e-5f);
	const G4ThreeVector tube_axis = tb->GetObjectRotationValue() * G4ThreeVector(0., 1., 0.);
	EXPECT_NEAR(glm::dot(glm::vec3(tube_axis.x(), tube_axis.y(), tube_axis.z()), dir), 1.f, 1e-5f);

	// Geant4's own navigation agrees: the detector's centre point is inside the detector
	G4Navigator navigator;
	navigator.SetWorldVolume(G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
	const G4VPhysicalVolume* at_centre = navigator.LocateGlobalPointAndSetup(G4ThreeVector(expected_centre.x, expected_centre.y, expected_centre.z) * m, nullptr, false, true);
	ASSERT_NE(at_centre, nullptr);
	EXPECT_EQ(at_centre->GetName(), "detector");

	RadiationSimulator::deinitialize();
}

TEST(CArmPlacement, DetectorIsFittedToTheBeam) {
	// a 34 x 34 cm plate at the Artis Q focus-isocentre distance of 78.5 cm; field sizes are given at the isocentre
	const double sod = 0.785 * m;
	const double min_distance = 0.45 * m, max_distance = 0.62 * m;
	auto entrance_after_fit = [&](const glm::vec3& centre, const glm::vec2& field, double lo, double hi, double& distance) {
		auto plate = make_box("plate", glm::vec3(0.17f, 0.0003f, 0.17f));
		plate->setType("ImageDetector");
		plate->setPosition(centre);
		Geant4::Mesh g4plate(plate, m);
		const RectangleSourceShape shape(field, static_cast<float>(sod / m));
		distance = g4plate.fitToBeam(sod, shape.getHalfTangents(), lo, hi);
		double entrance = INFINITY;
		for (const G4ThreeVector& c : g4plate.placedBoundingBoxCorners())
			entrance = std::min(entrance, c.y());
		return entrance;
	};
	double distance = 0.0;

	// a small field fits at any distance: the detector goes to its farthest position
	EXPECT_NEAR(entrance_after_fit(glm::vec3(0.f, 0.5f, 0.f), glm::vec2(0.10f), min_distance, max_distance, distance), max_distance, 1e-6);
	EXPECT_NEAR(distance, max_distance, 1e-6);
	// in between: the beam exactly covers the plate, 20 cm at the isocentre -> 34 cm at 78.5 + 54.95 cm
	EXPECT_NEAR(entrance_after_fit(glm::vec3(0.f, 0.5f, 0.f), glm::vec2(0.20f), min_distance, max_distance, distance), 0.5495 * m, 1e-3);
	// the longer field side limits it; the source frame's y runs along the plate's Z
	EXPECT_NEAR(entrance_after_fit(glm::vec3(0.f, 0.5f, 0.f), glm::vec2(0.10f, 0.20f), min_distance, max_distance, distance), 0.5495 * m, 1e-3);
	// a large field pulls the detector to its nearest position
	EXPECT_NEAR(entrance_after_fit(glm::vec3(0.f, 0.5f, 0.f), glm::vec2(0.30f), min_distance, max_distance, distance), min_distance, 1e-6);
	// off the axis, the nearest edge limits the beam: 15 cm reach -> 15 / 10 * 78.5 - 78.5 = 39.25 cm
	EXPECT_NEAR(entrance_after_fit(glm::vec3(0.02f, 0.5f, 0.f), glm::vec2(0.20f), 0.3 * m, max_distance, distance), 0.3925 * m, 1e-3);
	// a beam axis that misses the plate cannot be fitted
	EXPECT_THROW(entrance_after_fit(glm::vec3(0.2f, 0.5f, 0.f), glm::vec2(0.20f), min_distance, max_distance, distance), std::runtime_error);
}


TEST(CArmPlacement, GeometryFreePlacementMatchesTheGeant4Placement) {
	// the analytical simulator places meshes with CArm::place, the Geant4 binary with Geant4::Mesh::fitToBeam and
	// turnWithCArm: both must put a detector (fitted or not), a tube housing and the arm at the same corners
	struct Case { float phi, theta, psi; glm::vec3 detector_at; glm::vec2 field; bool fitted; };
	const Case cases[] = {
		{ 0.f, 90.f, 0.f, glm::vec3(0.f, 0.5f, 0.f), glm::vec2(0.10f), false },
		{ 30.f, 60.f, 25.f, glm::vec3(0.05f, 0.5f, -0.02f), glm::vec2(0.10f, 0.20f), false },
		{ 200.f, 120.f, -90.f, glm::vec3(0.02f, 0.5f, 0.03f), glm::vec2(0.20f, 0.12f), true },
		{ 90.f, 45.f, 180.f, glm::vec3(0.f, 0.48f, 0.f), glm::vec2(0.05f), true },
		{ 359.f, 179.f, 10.f, glm::vec3(-0.03f, 0.45f, 0.f), glm::vec2(0.21f), true },
	};
	const float source_distance = 0.785f;
	for (const Case& c : cases) {
		const glm::quat angles = glm::angleAxis(glm::radians(c.theta), glm::vec3(1.f, 0.f, 0.f)) * glm::angleAxis(glm::radians(c.phi), glm::vec3(0.f, 1.f, 0.f));
		const glm::vec3 dir = angles * glm::vec3(0.f, 0.f, -1.f);
		XRaySource source(50e3f, std::make_unique<RectangleSourceShape>(c.field, source_distance));
		source.setTransform(-dir * source_distance, dir);
		source.setCArmPsi(glm::radians(c.psi));
		const glm::vec2 half_tangents = source.getShape()->getHalfTangents();

		auto plate = [&](const std::string& name) {
			auto mesh = make_box(name, glm::vec3(0.17f, 0.0003f, 0.17f));
			mesh->setType("ImageDetector");
			mesh->setPosition(c.detector_at);
			if (c.fitted)
				mesh->setIsocenterDistanceRange(glm::vec2(0.45f, 0.62f));
			return mesh;
		};
		auto housing = [&](const std::string& name) {
			auto mesh = make_box(name, glm::vec3(0.05f, 0.03f, 0.05f));
			mesh->setType("XRayTube");
			mesh->setPosition(glm::vec3(0.01f, -0.05f, 0.f));
			return mesh;
		};
		auto arm = [&](const std::string& name) {
			auto mesh = make_box(name, glm::vec3(0.06f, 0.10f, 0.45f));
			mesh->setType("CArm");
			mesh->setPosition(glm::vec3(0.f, 0.f, 0.9f));
			return mesh;
		};
		const auto base_corners = [](const Mesh& mesh) { return CArm::placed_bounding_box_corners(mesh); };

		// Geant4 path, in millimetres
		auto g4_plate_mesh = plate("g4_plate");
		auto g4_tube_mesh = housing("g4_tube");
		auto g4_arm_mesh = arm("g4_arm");
		Geant4::Mesh g4_plate(g4_plate_mesh, m), g4_tube(g4_tube_mesh, m), g4_arm(g4_arm_mesh, m);
		const glm::vec3 focal_spot = source.getLocation();
		if (c.fitted)
			g4_plate.fitToBeam(source_distance * m, half_tangents, 0.45 * m, 0.62 * m);
		g4_plate.turnWithCArm(source.getCArmRotation(), G4ThreeVector());
		g4_tube.turnWithCArm(source.getCArmRotation(), G4ThreeVector(focal_spot.x, focal_spot.y, focal_spot.z) * m);
		g4_arm.turnWithCArm(source.getCArmStructureRotation(), G4ThreeVector());

		// Geant4-free path, in metres
		auto plate_mesh = plate("plate");
		auto tube_mesh = housing("tube");
		auto arm_mesh = arm("arm");
		const CArm::Pose pose{ source.getCArmRotation(), focal_spot, half_tangents, source.getCArmStructureRotation() };
		const std::optional<float> distance = CArm::place(*plate_mesh, pose);
		EXPECT_EQ(distance.has_value(), c.fitted);
		EXPECT_FALSE(CArm::place(*tube_mesh, pose).has_value());
		EXPECT_FALSE(CArm::place(*arm_mesh, pose).has_value());

		const auto expect_same_corners = [&](const std::vector<G4ThreeVector>& g4, const Mesh& mesh) {
			const std::vector<glm::vec3> corners = base_corners(mesh);
			ASSERT_EQ(g4.size(), corners.size());
			for (size_t i = 0; i < corners.size(); i++)
				EXPECT_LT(glm::length(to_m(g4[i]) - corners[i]), 1e-5f) << "phi " << c.phi << " theta " << c.theta << " corner " << i << " of " << mesh.getName();
		};
		expect_same_corners(g4_plate.placedBoundingBoxCorners(), *plate_mesh);
		expect_same_corners(g4_tube.placedBoundingBoxCorners(), *tube_mesh);
		expect_same_corners(g4_arm.placedBoundingBoxCorners(), *arm_mesh);
		// psi turns the arm about the beam axis: its centre keeps its distance from the beam and from the isocentre
		const glm::vec3 arm_centre = arm_mesh->getPosition();
		EXPECT_NEAR(glm::length(arm_centre - dir * glm::dot(arm_centre, dir)), 0.9f, 1e-4f);
		const glm::vec3 unturned = source.getCArmRotation() * glm::vec3(0.f, 0.f, 0.9f);
		EXPECT_NEAR(glm::degrees(std::atan2(glm::length(glm::cross(arm_centre, unturned)), glm::dot(arm_centre, unturned))), std::abs(std::remainder(c.psi, 360.f)), 1e-2f);
		if (c.fitted) {
			// the entrance face lies at the fitted distance along the beam from the isocentre
			float entrance = 1e9f;
			for (const glm::vec3& corner : base_corners(*plate_mesh))
				entrance = std::min(entrance, glm::dot(corner, dir));
			EXPECT_NEAR(entrance, *distance, 1e-4f);
		}
	}
}
