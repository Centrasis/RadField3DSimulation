#pragma once
#include "World.hpp"

class G4Box;
class G4Material;
class G4LogicalVolume;
class RadiationFieldDetector;

namespace RadiationSimulation {

	/**
	 * @class G4World
	 * @brief A class to abstract Geant4 specific components and integrate them into the simulation framework.
	 *
	 * The G4World class is responsible for managing Geant4 specific objects such as G4Box, G4Material, and G4LogicalVolume.
	 * It extends the World class and provides an interface to set and get radiation sources, geometries, and radiation field detectors.
	 * This class ensures that the underlying Geant4 components are properly initialized and accessible within the simulation framework.
	 */
	class G4World : public World {
	protected:
		const std::shared_ptr<World> raw_world;
		std::shared_ptr<G4Box> box;
		std::shared_ptr<G4Material> material;
		std::shared_ptr<G4LogicalVolume> volume;

	public:
		/**
		 * @brief Constructor to initialize the G4World with a raw World object.
		 * @param raw_world A shared pointer to the raw World object.
		 */
		G4World(std::shared_ptr<World> raw_world);

		/**
		 * @brief Static method to initialize the Geant4 components.
		 * @param box A shared pointer to a G4Box object.
		 * @param material A shared pointer to a G4Material object.
		 * @param volume A shared pointer to a G4LogicalVolume object.
		 */
		static void initialize(std::shared_ptr<G4Box> box, std::shared_ptr<G4Material> material, std::shared_ptr<G4LogicalVolume> volume);

		/**
		 * @brief Static method to get the instance of G4World.
		 * @return A shared pointer to the G4World instance.
		 */
		static std::shared_ptr<G4World> Get();

		virtual void set_radiation_source(std::shared_ptr<RadiationSource> source) override { this->raw_world->set_radiation_source(source); }
		virtual void set_geometries(const std::vector<std::shared_ptr<Mesh>>& geoms) override { this->raw_world->set_geometries(geoms); }
		virtual std::shared_ptr<RadiationSource> get_radiation_source() const override { return this->raw_world->get_radiation_source(); };
		virtual std::shared_ptr<Mesh> get_patient() const override { return this->raw_world->get_patient(); };
		virtual const std::vector<std::shared_ptr<Mesh>>& get_geometries() const override { return this->raw_world->get_geometries(); };
		virtual void set_radiation_field_detector(std::shared_ptr<RadiationFieldDetector> field_detector) override { this->raw_world->set_radiation_field_detector(field_detector); }
		virtual std::shared_ptr<RadiationFieldDetector> get_radiation_field_detector() const override { return this->raw_world->get_radiation_field_detector(); }

		/**
		 * @brief Get the G4Box object.
		 * @return A shared pointer to the G4Box object.
		 */
		std::shared_ptr<G4Box> get_box() const { return this->box; }

		/**
		 * @brief Get the G4Material object.
		 * @return A shared pointer to the G4Material object.
		 */
		std::shared_ptr<G4Material> get_material() const { return this->material; }

		/**
		 * @brief Get the G4LogicalVolume object.
		 * @return A shared pointer to the G4LogicalVolume object.
		 */
		std::shared_ptr<G4LogicalVolume> get_volume() const { return this->volume; }
	};
}