// RadField3D-Analytical — analytical ray-tracing X-ray field simulator.
//
// Same CLI and .rf3 output contract as the Geant4 RadField3D binary, so it is a drop-in
// --binary for tools/create_dataset.py. It produces LOW-ACCURACY fields (analytical direct beam
// + single-scatter estimate) in a fraction of a second on the GPU, for pre-flighting network
// trainability before a full Monte-Carlo run. Analytic fields are marked in the metadata software
// section as "RadField3D-Analytical" so they can never be confused with Monte-Carlo ones.
//
// The tube/source setup (spectrum loading + random energy sampling, source shapes, placement) reuses
// the exact classes the Geant4 binary uses (RadField3D.cpp); only the transport is new (GPU tracer).
#include "analytical/AnalyticalTracer.hpp"
#include "analytical/AnalyticalVoxelizer.hpp"
#include "GeometryLoader.hpp"
#include "RadiationSource.hpp"
#include "Voxelization.hpp"

#include <radfiled3d/radiation_field.hpp>
#include <radfiled3d/voxel_grid.hpp>
#include <radfiled3d/voxel.hpp>
#include <radfiled3d/storage/radiation_field_store.hpp>
#include <radfiled3d/storage/types.hpp>
#include <radfiled3d/helpers/typing.hpp>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#if defined _WIN32 || defined _WIN64
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;
using namespace RadiationSimulation::Analytical;

namespace {

struct Cli {
    std::string out, spectrum, geom, geom_desc, source_shape = "rectangle";
    double max_energy = 0.0, energy_resolution = 1000.0;
    float source_phi = 0.f, source_theta = 0.f, source_distance = 1.f, voxel_dim = 0.1f;
    glm::vec3 world_dim = glm::vec3(1.f);
    glm::vec2 opening = glm::vec2(0.f);  // cone: (deg, 0); rectangle: (w, h) in metres
    float particles = 1e6f;
    float focal_spot = 0.01f;  // effective focal-spot size (m) — beam-edge penumbra blur
};

// Channel "geometry" on the field's grid: one 8-bit layer per mesh Type (255 = voxel overlaps a mesh of that type).
// Meshes are placed as for the density grid (AnalyticalVoxelizer): scale and rotation about the mesh origin,
// translations composed additively through the parents.
void add_geometry_channel(radfiled3d::CartesianRadiationField& field, const std::vector<std::shared_ptr<Mesh>>& meshes)
{
    if (meshes.empty())
        return;
    auto channel = std::static_pointer_cast<radfiled3d::VoxelGridBuffer>(field.add_channel("geometry"));
    const glm::uvec3 counts = channel->get_voxel_counts();
    const double voxel_size = channel->get_voxel_dimensions().x;
    const Voxelization::VoxelGrid grid{ counts, voxel_size, -0.5 * glm::dvec3(counts) * voxel_size };

    std::function<void(const std::shared_ptr<Mesh>&, const glm::vec3&)> place = [&](const std::shared_ptr<Mesh>& mesh, const glm::vec3& parent_position) {
        const glm::vec3 position = parent_position + mesh->getPosition();
        const glm::dmat4 local_to_world = glm::translate(glm::dmat4(1.0), glm::dvec3(position))
            * glm::dmat4(glm::mat4_cast(mesh->getRotation()))
            * glm::scale(glm::dmat4(1.0), glm::dvec3(mesh->getScale()));
        if (!channel->has_layer(mesh->getType()))
            channel->add_layer<uint8_t>(mesh->getType(), 0, "occupancy");
        Voxelization::mark_overlapping_voxels(Voxelization::TransformedMeshSurface(*mesh, local_to_world), grid, channel->get_layer<uint8_t>(mesh->getType()));
        for (const auto& child : mesh->getChildren())
            place(child, position);
    };
    for (const auto& mesh : meshes)
        place(mesh, glm::vec3(0.f));
}

}  // namespace

int main(int argc, char* argv[])
{
try {
    const long long start_time =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    Cli cli;
    std::vector<std::string> ignored;
    for (int i = 1; i < argc; ) {
        std::string arg = argv[i];
        auto next = [&]() -> std::string { return (i + 1 < argc) ? std::string(argv[i + 1]) : std::string(); };
        if (arg == "--out") { cli.out = next(); i += 2; }
        else if (arg == "--max-energy") { cli.max_energy = std::stod(next()); i += 2; }
        else if (arg == "--source-phi") { cli.source_phi = std::stof(next()); i += 2; }
        else if (arg == "--source-theta") { cli.source_theta = std::stof(next()); i += 2; }
        else if (arg == "--source-distance") { cli.source_distance = std::stof(next()); i += 2; }
        else if (arg == "--voxel-dim") { cli.voxel_dim = std::stof(next()); i += 2; }
        else if (arg == "--energy-resolution") { cli.energy_resolution = std::stod(next()); i += 2; }
        else if (arg == "--source-shape") { cli.source_shape = next(); i += 2; }
        else if (arg == "--spectrum") { cli.spectrum = next(); i += 2; }
        else if (arg == "--geom") { cli.geom = next(); i += 2; }
        else if (arg == "--geom-desc") { cli.geom_desc = next(); i += 2; }
        else if (arg == "--particles") { cli.particles = std::stof(next()); i += 2; }
        else if (arg == "--focal-spot") { cli.focal_spot = std::stof(next()); i += 2; }
        else if (arg == "--world-dim") {
            std::stringstream ss(next()); float x, y, z; ss >> x >> y >> z; cli.world_dim = glm::vec3(x, y, z); i += 2;
        }
        else if (arg == "--source-opening-angle") {
            std::stringstream ss(next()); float x = 0.f, y = 0.f; ss >> x; ss >> y; cli.opening = glm::vec2(x, y); i += 2;
        }
        // accepted for CLI parity with the Geant4 binary, but analytically irrelevant:
        else if (arg == "--angular-resolution" || arg == "--tracing-algorithm" || arg == "--world-material") { i += 2; }
        else { ignored.push_back(arg); i += 1; }
    }
    if (!ignored.empty()) {
        std::cout << "Ignoring unsupported arguments:";
        for (const auto& a : ignored) std::cout << " " << a;
        std::cout << std::endl;
    }
    if (cli.out.empty() || cli.spectrum.empty() || cli.max_energy <= 0.0)
        throw std::runtime_error("--out, --spectrum and --max-energy are required");

    fs::path out_path = cli.out;
    if (out_path.extension() != ".rf3") out_path.replace_extension(".rf3");

    const int nx = (int)std::llround(cli.world_dim.x / cli.voxel_dim);
    const int ny = (int)std::llround(cli.world_dim.y / cli.voxel_dim);
    const int nz = (int)std::llround(cli.world_dim.z / cli.voxel_dim);
    const int bins = std::max(1, (int)std::llround(cli.max_energy / cli.energy_resolution));
    const double bin_width_ev = cli.max_energy / bins;

    // --- source shape (reuses the Geant4 binary's shapes) ---
    std::unique_ptr<ISourceShape> shape;
    radfiled3d::typing::FieldShape field_shape;
    if (cli.source_shape == "rectangle") {
        shape = std::make_unique<RectangleSourceShape>(glm::vec2(cli.opening.x, cli.opening.y), cli.source_distance);
        field_shape = radfiled3d::typing::FieldShape::Rectangle;
    } else if (cli.source_shape == "cone") {
        shape = std::make_unique<ConeSourceShape>(cli.opening.x);
        field_shape = radfiled3d::typing::FieldShape::Cone;
    } else {
        throw std::runtime_error("analytical simulator supports source-shape rectangle or cone, got " + cli.source_shape);
    }

    // --- spectrum: the loaded distribution, integrated exactly over the histogram bins (no sampling) ---
    auto pdf = SpectrumLoader::LoadSpectrum(cli.spectrum);
    auto source = std::make_shared<XRaySpectrumSource>(pdf, std::move(shape), 0.f, (float)cli.max_energy);
    const int tube_bins = (int)source->getGeneratedSpectrumBins();
    const float tube_bin_width_ev = source->getGeneratedSpectrumBinWidth_eV();

    // tube spectrum metadata: same binning as the Monte-Carlo binary records (bin i centred at i * width, the last
    // bin also holds everything above)
    std::vector<float> tube_weights(tube_bins, 0.f);
    for (int tb = 0; tb < tube_bins; ++tb) {
        const double upper = (tb + 1 < tube_bins) ? (tb + 0.5) * tube_bin_width_ev : std::numeric_limits<double>::max();
        tube_weights[tb] = (float)pdf->probability((tb - 0.5) * tube_bin_width_ev, upper);
    }

    // field-bin weights (physics): bin fb covers [fb, fb + 1) * bin width, the last bin also holds everything above
    std::vector<float> field_weights(bins, 0.f);
    for (int fb = 0; fb < bins; ++fb) {
        const double upper = (fb + 1 < bins) ? (fb + 1) * bin_width_ev : std::numeric_limits<double>::max();
        field_weights[fb] = (float)pdf->probability(fb * bin_width_ev, upper);
    }
    const float spectrum_max_energy_eV = (float)pdf->max();

    // --- source placement, identical to RadField3D.cpp (phi about Y, theta about X) ---
    const glm::quat alpha = glm::angleAxis(glm::radians(cli.source_phi), glm::vec3(0.f, 1.f, 0.f));
    const glm::quat beta = glm::angleAxis(glm::radians(cli.source_theta), glm::vec3(1.f, 0.f, 0.f));
    const glm::vec3 dir = (beta * alpha) * glm::vec3(0.f, 0.f, -1.f);
    source->setTransform(-dir * cli.source_distance, dir);  // reuse: origin + collimation frame
    const glm::vec3 origin = source->getLocation();
    const glm::mat3 rot = glm::mat3_cast(source->getRotation());
    const glm::vec3 e1 = rot * glm::vec3(1.f, 0.f, 0.f);
    const glm::vec3 e2 = rot * glm::vec3(0.f, 1.f, 0.f);

    // --- geometry -> material density grid ---
    std::vector<float> density((size_t)nx * ny * nz, 0.f);
    std::vector<std::shared_ptr<Mesh>> meshes;
    if (!cli.geom.empty()) {
        meshes = GeometryLoader::Load(cli.geom, cli.geom_desc);
        density = AnalyticalVoxelizer::voxelize_density(meshes, glm::ivec3(nx, ny, nz), cli.voxel_dim);
    }

    AnalyticalParams p;
    p.nx = nx; p.ny = ny; p.nz = nz; p.voxel_m = cli.voxel_dim; p.bins = bins; p.bin_width_ev = (float)bin_width_ev;
    p.origin[0] = origin.x; p.origin[1] = origin.y; p.origin[2] = origin.z;
    p.direction[0] = dir.x; p.direction[1] = dir.y; p.direction[2] = dir.z;
    p.e1[0] = e1.x; p.e1[1] = e1.y; p.e1[2] = e1.z;
    p.e2[0] = e2.x; p.e2[1] = e2.y; p.e2[2] = e2.z;
    p.distance = cli.source_distance;
    p.particles = cli.particles;
    p.focal_spot_m = cli.focal_spot;
    if (field_shape == radfiled3d::typing::FieldShape::Rectangle) {
        const glm::vec2 rect = static_cast<const RectangleSourceShape*>(source->getShape())->getFieldSizeMeters();
        p.shape = 0; p.rect_w = rect.x; p.rect_h = rect.y;
    } else {
        const float deg = static_cast<const ConeSourceShape*>(source->getShape())->getOpeningAngleDegrees();
        p.shape = 1; p.cone_half_rad = glm::radians(deg) * 0.5f;
    }

    // --- allocate the .rf3 field first (same channel/layer layout the Geant4 detector writes) and
    //     let the tracer write its results straight into the field-owned VoxelLayer buffers ---
    auto field = std::make_shared<radfiled3d::CartesianRadiationField>(
        glm::vec3(cli.world_dim.x, cli.world_dim.y, cli.world_dim.z),
        glm::vec3(cli.voxel_dim, cli.voxel_dim, cli.voxel_dim));
    const float bin_width_mev = (float)(bin_width_ev * 1e-6);

    auto make_channel = [&](const char* name) -> AnalyticalChannelBuffers {
        auto* buf = static_cast<radfiled3d::VoxelGridBuffer*>(field->add_channel(name).get());
        buf->add_layer<float>("error", 1.f, "Variance");
        buf->add_layer<float>("flux", 0.f, "counts");
        buf->add_custom_layer<radfiled3d::HistogramVoxel<float>>(
            "spectrum", radfiled3d::HistogramVoxel<float>(bins, bin_width_mev, nullptr), 0.f, "eV");
        // VoxelLayer data buffers are contiguous (x-fastest); voxel 0's data is the buffer base.
        AnalyticalChannelBuffers b;
        b.flux = &buf->get_voxel_flat<radfiled3d::ScalarVoxel<float>>("flux", 0).get_data();
        b.error = &buf->get_voxel_flat<radfiled3d::ScalarVoxel<float>>("error", 0).get_data();
        b.spectrum = &buf->get_voxel_flat<radfiled3d::HistogramVoxel<float>>("spectrum", 0).get_data();
        return b;
    };
    AnalyticalOutput out;
    out.scatter = make_channel("scatter_field");
    out.direct = make_channel("direct_beam");

    run_analytical(p, density, field_weights, out);

    // --- metadata (marks these as analytical fields) ---
    auto metadata = std::make_shared<radfiled3d::storage::v1::RadiationFieldMetadata>(
        radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation(
            (size_t)cli.particles, cli.geom, "analytic-raytracing",
            radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Simulation::XRayTube(
                dir, origin, spectrum_max_energy_eV, cli.spectrum)),
        radfiled3d::storage::filed_types::v1::RadiationFieldMetadataHeader::Software(
            "RadField3D-Analytical", "1.0.0-analytic", "https://github.com/Centrasis/RadField3DSimulation", "analytic"));

    metadata->set_dynamic_custom_metadata<radfiled3d::HistogramVoxel<float>>(
        "tube_spectrum", radfiled3d::HistogramVoxel<float>(tube_bins, tube_bin_width_ev, nullptr));
    std::memcpy(metadata->get_dynamic_metadata<radfiled3d::HistogramVoxel<float>>("tube_spectrum").get_histogram().data(),
                tube_weights.data(), sizeof(float) * tube_bins);
    const long long end_time =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    metadata->add_dynamic_metadata<uint64_t>("simulation_duration_s", (uint64_t)(end_time - start_time));
    metadata->add_dynamic_metadata<uint8_t>("xray_field_shape", static_cast<uint8_t>(field_shape));
    if (field_shape == radfiled3d::typing::FieldShape::Rectangle)
        metadata->add_dynamic_metadata<glm::vec2>("xray_tube_field_rect_dimensions_m", glm::vec2(p.rect_w, p.rect_h));
    else
        metadata->add_dynamic_metadata<float>("xray_tube_opening_angle_deg", cli.opening.x);

    add_geometry_channel(*field, meshes);

    if (fs::exists(out_path)) fs::remove(out_path);
    radfiled3d::storage::FieldStore::store(field, metadata, out_path.string(), radfiled3d::storage::StoreVersion::V1);
    std::cout << "Wrote field to: " << out_path.string() << std::endl;
    return 0;
}
catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
}
