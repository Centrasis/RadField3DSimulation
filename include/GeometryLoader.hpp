#pragma once
#include "Geometry.hpp"
#include <string>
#include <memory>


namespace RadiationSimulation {
	class GeometryLoader {
	public:
		static std::vector<std::shared_ptr<Mesh>> Load(const std::string& path);
	};
}