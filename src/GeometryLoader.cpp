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


void SetupMesh(json& mesh_desc, std::shared_ptr<Mesh> mesh, const std::map<std::string, std::shared_ptr<Mesh>>& all_meshes) {
    bool has_a_patient = false;
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

    if (mesh_desc.find("Patient") != mesh_desc.end()) {
        bool is_patient = mesh_desc["Patient"].get<bool>();
        if (is_patient && has_a_patient)
            throw std::runtime_error("A patient mesh was already defined! There cannot be more than one patient!");
        if (is_patient) {
            mesh->markAsPatient();
            has_a_patient = true;
        }
    }

    if (mesh_desc.find("MaterialName") != mesh_desc.end()) {
		std::string material_name = mesh_desc["MaterialName"].get<std::string>();
		mesh->attachMaterialName(material_name);
	}

    if (mesh_desc.find("Children") != mesh_desc.end()) {
        auto& children = mesh_desc["Children"];
        for (auto& [child_name, child] : children.items()) {
            std::shared_ptr<Mesh> child_mesh = all_meshes.find(child_name)->second;
            SetupMesh(child, child_mesh, all_meshes);
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
            glm::vec2(rotation_info["Alpha"].get<float>(), rotation_info["Beta"].get<float>())
        );
    }
}

std::vector<std::shared_ptr<Mesh>> GeometryLoader::Load(const std::string& path)
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

    const std::string desc_file_path = path.substr(0, path.find_last_of(".")) + ".desc";

    if (fs::exists(desc_file_path)) {
        std::ifstream desc_file(desc_file_path);

        json data;
        desc_file >> data;

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
            SetupMesh(m_data, mesh, meshes);
        }
    } else {
        for (auto& [name, m] : meshes) {
			root_meshes.push_back(m);
		}
	}

    return root_meshes;
}
