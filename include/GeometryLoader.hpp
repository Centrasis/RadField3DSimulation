#pragma once
#include "Geometry.hpp"
#include <string>
#include <memory>


namespace RadiationSimulation {
	class GeometryLoader {
	public:
		/** 
		 * Loads a geometry from a file and returns a vector of Mesh objects.
		 * @param path The path to the geometry file.
		 * @param description_file Optional JSON file containing mesh descriptions. Default is an empty string which means an optional description file with the same basename like path, but ending with '.desc' will be used.
		 * @return A vector of shared pointers to Mesh objects.
		 * @throws std::runtime_error if the file cannot be loaded or if there are issues with the mesh descriptions.
		 */
		static std::vector<std::shared_ptr<Mesh>> Load(const std::string& path, std::string description_file = "");
	};
}