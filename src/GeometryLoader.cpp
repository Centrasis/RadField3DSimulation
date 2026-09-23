#include "GeometryLoader.hpp"
#include <assimp/cimport.h>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <glm/vec2.hpp>
#if defined _WIN32 || defined _WIN64
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif
#include <stdexcept>

using json = nlohmann::json;
using namespace RadiationSimulation;
using namespace RadiationSimulation::Geometry;


// Applies a mesh description to `mesh` and its children. A mesh without its own Type takes the type of its parent
// (`parent_type`, empty for root meshes, which then keep Mesh::DEFAULT_TYPE). `declared_patients` counts the meshes that
// declare the patient type themselves.
void SetupMesh(json& mesh_desc, std::shared_ptr<Mesh> mesh, const std::map<std::string, std::shared_ptr<Mesh>>& all_meshes, const std::string& parent_type, size_t& declared_patients) {
    if (mesh_desc.find("Transform") != mesh_desc.end()) {
        auto& transform_info = mesh_desc["Transform"];
        if (transform_info.find("Rotation") != transform_info.end()) {
            auto& info = transform_info["Rotation"];
            float x = info["X"].get<float>();
            float y = info["Y"].get<float>();
            float z = info["Z"].get<float>();
            mesh->setRotation(glm::vec3(x, y, z));
        }
        if (transform_info.find("Translation") != transform_info.end()) {
            auto& info = transform_info["Translation"];
            float x = info["X"].get<float>();
            float y = info["Y"].get<float>();
            float z = info["Z"].get<float>();
            mesh->setPosition(glm::vec3(x, y, z));
        }
        if (transform_info.find("Scale") != transform_info.end()) {
            auto& info = transform_info["Scale"];
            float x = info["X"].get<float>();
            float y = info["Y"].get<float>();
            float z = info["Z"].get<float>();
            mesh->setScale(glm::vec3(x, y, z));
        }
    }

    if (mesh_desc.find("Type") != mesh_desc.end()) {
        try {
            mesh->setType(mesh_desc["Type"].get<std::string>());
        }
        catch (const std::invalid_argument& e) {
            throw std::runtime_error("Invalid Type of mesh \"" + mesh->getName() + "\": " + e.what());
        }
        if (mesh->isPatient())
            declared_patients++;
    }
    // legacy descriptions flag the patient with a boolean instead of a Type
    else if (mesh_desc.find("Patient") != mesh_desc.end() && mesh_desc["Patient"].get<bool>()) {
        mesh->setType(Mesh::PATIENT_TYPE);
        declared_patients++;
    }
    else if (!parent_type.empty()) {
        mesh->setType(parent_type);
    }

    if (mesh_desc.find("MaterialName") != mesh_desc.end()) {
		std::string material_name = mesh_desc["MaterialName"].get<std::string>();
		mesh->attachMaterialName(material_name);
	}

    if (mesh_desc.find("Children") != mesh_desc.end()) {
        auto& children = mesh_desc["Children"];
        for (auto& [child_name, child] : children.items()) {
            std::shared_ptr<Mesh> child_mesh = all_meshes.find(child_name)->second;
            SetupMesh(child, child_mesh, all_meshes, mesh->getType(), declared_patients);
            mesh->addChild(child_mesh);
        }
    }

    if (mesh_desc.find("Source") != mesh_desc.end() && mesh_desc["Source"] == true) {
		if (mesh_desc.find("SourceOffsets") == mesh_desc.end()) {
			throw std::runtime_error("Source mesh with name: \"" + mesh->getName() + "\" must have SourceOffsets defined!");
		}

		auto& source_offsets_info = mesh_desc["SourceOffsets"];
		auto& translation_info = source_offsets_info["Translation"];
		auto& rotation_info = source_offsets_info["Rotation"];

        mesh->markAsSource(
            translation_info["ConcentricDistance"].get<float>(),
            glm::vec2(rotation_info["Phi"].get<float>(), rotation_info["Theta"].get<float>())
        );
    }
}

std::vector<std::shared_ptr<Mesh>> GeometryLoader::Load(const std::string& path, std::string description_file)
{
    const aiScene* scene = aiImportFile(path.c_str(), aiProcessPreset_TargetRealtime_MaxQuality);
    if (!scene) {
        std::string error_msg = "Could not load file: " + path;
        throw std::runtime_error(error_msg.c_str());
    }
    std::map<std::string, std::shared_ptr<Mesh>> meshes;
    std::vector<std::shared_ptr<Mesh>> root_meshes;

    for (size_t mid = 0; mid < scene->mNumMeshes; mid++) {
        aiMesh* raw_mesh = scene->mMeshes[mid];
        std::vector<glm::vec3> vertices(raw_mesh->mNumVertices);
        for (size_t vid = 0; vid < raw_mesh->mNumVertices; vid++) {
            aiVector3D& v = raw_mesh->mVertices[vid];
            vertices[vid] = glm::vec3(v[0], v[1], v[2]);
        }

        std::vector<Face*> faces(raw_mesh->mNumFaces);
        for (size_t fid = 0; fid < raw_mesh->mNumFaces; fid++) {
            aiFace& face = raw_mesh->mFaces[fid];
            switch (face.mNumIndices) {
            case 3:
                faces[fid] = new TriFace(glm::uvec3(face.mIndices[0], face.mIndices[1], face.mIndices[2]));
                break;
            case 4:
                faces[fid] = new QuadFace(glm::uvec4(face.mIndices[0], face.mIndices[1], face.mIndices[2], face.mIndices[3]));
                break;
            default:
                throw std::runtime_error("Invalid face indices count!");
            }
        }

        const std::string m_name(raw_mesh->mName.C_Str());
        meshes.insert({ m_name, std::make_shared<Mesh>(vertices, faces, m_name) });
    }

    if (description_file.size() == 0)
        description_file = path.substr(0, path.find_last_of(".")) + ".desc";

    if (fs::exists(description_file)) {
        std::ifstream desc_file(description_file);

        json data;
        desc_file >> data;

        size_t declared_patients = 0;
        for (auto& [name, m_data] : data.items()) {
            if (meshes.find(name) == meshes.end()) {
                std::cout << "Meshes were:" << std::endl;
                for (auto& [m_name, m] : meshes) {
					std::cout << m_name << std::endl;
				}
                throw std::runtime_error("Mesh not found: " + name);
            }
            std::shared_ptr<Mesh> mesh = meshes.find(name)->second;
            root_meshes.push_back(mesh);
            SetupMesh(m_data, mesh, meshes, "", declared_patients);
        }

        if (declared_patients > 1)
            throw std::runtime_error("More than one mesh declares the Type \"" + std::string(Mesh::PATIENT_TYPE) + "\"! There cannot be more than one patient!");
    } else {
        for (auto& [name, m] : meshes) {
			root_meshes.push_back(m);
		}
	}

    return root_meshes;
}
