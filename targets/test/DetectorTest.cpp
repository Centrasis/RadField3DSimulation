#include "gtest/gtest.h"
#include <Geant4/G4RadiationFieldDetector.hpp>
#include <RadFiled3D/RadiationField.hpp>
#include <G4RunManager.hh>
#include <G4VUserPhysicsList.hh>
#include <G4Geantino.hh>
#include <glm/glm.hpp>
#include <malloc.h>
#include <memory>

using namespace RadiationSimulation;

// G4UserSteppingAction (the detector's base) refuses construction unless a physics list is
// registered with a run manager. Needs the Geant4 data environment (G4ENSDFSTATEDATA etc.).
namespace {
	class MinimalPhysicsList : public G4VUserPhysicsList {
	protected:
		void ConstructParticle() override { G4Geantino::GeantinoDefinition(); }
		void ConstructProcess() override { AddTransportation(); }
	};

	class G4Bootstrap : public ::testing::Environment {
	public:
		void SetUp() override {
			this->run_manager = new G4RunManager();
			this->run_manager->SetUserInitialization(new MinimalPhysicsList());
		}
		G4RunManager* run_manager = nullptr;
	};

	const auto* g4_bootstrap = ::testing::AddGlobalTestEnvironment(new G4Bootstrap());
}


TEST(NormalizedFieldCopy, WorksWithoutAngularFlux) {
	// The periodic auto-save path: must not throw when no angular resolution was requested
	// (the angular_flux layer does not exist then).
	G4RadiationFieldDetector det(glm::vec3(0.4f), glm::vec3(0.02f), 32, 150000.0 / 32.0);
	std::shared_ptr<RadFiled3D::IRadiationField> copy;
	ASSERT_NO_THROW(copy = det.get_normalized_field_copy());
	auto cart = std::dynamic_pointer_cast<RadFiled3D::CartesianRadiationField>(copy);
	ASSERT_NE(cart, nullptr);
	for (const char* channel : {"scatter_field", "direct_beam"}) {
		auto ch = cart->get_channel(channel);
		EXPECT_TRUE(ch->has_layer("flux"));
		EXPECT_TRUE(ch->has_layer("spectrum"));
		EXPECT_TRUE(ch->has_layer("error"));
		EXPECT_FALSE(ch->has_layer("angular_flux"));
	}
}

TEST(NormalizedFieldCopy, RepeatedCopiesDoNotCorruptTheSource) {
	// Normalization must only touch the copy — a second call has to work identically.
	G4RadiationFieldDetector det(glm::vec3(0.4f), glm::vec3(0.02f), 32, 150000.0 / 32.0);
	ASSERT_NO_THROW(det.get_normalized_field_copy());
	ASSERT_NO_THROW(det.get_normalized_field_copy());
}

TEST(FieldDimensions, VoxelCountsSurviveFloatTruncation) {
	// 2.3f / 0.1f = 22.9999971 in float: without the epsilon guard in CartesianRadiationField
	// this truncates to 22 voxels per axis. The user's production case 2.3 m / 0.02 m (= exactly
	// 115.0f) is safe either way; this pins the guard for the ratios that are not.
	G4RadiationFieldDetector det(glm::vec3(2.3f), glm::vec3(0.1f), 8, 150000.0 / 8.0);
	auto copy = std::dynamic_pointer_cast<RadFiled3D::CartesianRadiationField>(det.get_normalized_field_copy());
	ASSERT_NE(copy, nullptr);
	EXPECT_EQ(copy->get_voxel_counts().x, 23u);
	EXPECT_EQ(copy->get_voxel_counts().y, 23u);
	EXPECT_EQ(copy->get_voxel_counts().z, 23u);
}

TEST(DetectorMemory, FootprintWithinLayerBudget) {
	// The field layers of a 48^3 two-channel detector are ~48 MB. The detector bookkeeping
	// (per-voxel spectrum-variance trackers, per-voxel mutexes, construction-time duplication)
	// must not dwarf the data it annotates; 3x the layer bytes is a generous budget. This test
	// documents the current per-voxel-statistics memory bloat (see docs/bugscan report): it FAILS
	// until ChannelBuffers is slimmed down and passes afterwards.
	const size_t n = 48ull * 48ull * 48ull;
	// 2 channels x (double flux + float error + 32-bin double spectrum): scoring accumulates in
	// double (float += 1 saturates at 2^24 counts); the stored field is fp32.
	const size_t layer_bytes_per_voxel = 2 * (8 + 4 + 32 * 8);
	const size_t budget = 3 * n * layer_bytes_per_voxel;

	struct mallinfo2 before = mallinfo2();
	auto det = std::make_unique<G4RadiationFieldDetector>(glm::vec3(0.96f), glm::vec3(0.02f), 32, 150000.0 / 32.0);
	struct mallinfo2 after = mallinfo2();

	const size_t growth = after.uordblks - before.uordblks;
	EXPECT_LT(growth, budget) << "detector allocated " << growth / (1024.0 * 1024.0)
		<< " MiB for a field whose layers need " << (n * layer_bytes_per_voxel) / (1024.0 * 1024.0) << " MiB";
}
