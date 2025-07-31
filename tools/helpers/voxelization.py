import trimesh
import rtree # needed to be included for voxelization
import numpy as np
import json


class VoxelizationHelper:
    @staticmethod
    def load_mesh(geom_file: str) -> trimesh.Trimesh:
        """
        Load a mesh from a geometry file.

        :param geom_file: Path to the geometry file.
        :return: Loaded mesh as a trimesh object.
        """
        loaded_meshes = trimesh.load(geom_file)
        if isinstance(loaded_meshes, trimesh.Scene):
            return trimesh.util.concatenate(
            [
                trimesh.Trimesh(mesh.vertices, mesh.faces, face_normals=mesh.face_normals, vertex_normals=mesh.vertex_normals)
                for mesh in loaded_meshes.geometry.values()
            ]
        )
        else:
            return loaded_meshes

    @staticmethod
    def generate_voxelgrid_with_geometry(geom_file: str, voxel_size: float, grid_size: tuple[int, int, int], description_file: str = None) -> np.ndarray:
        """
        Generate a voxel grid from a geometry file.

        :param geom_file: Path to the geometry file.
        :param voxel_size: Size of each voxel.
        :return: Voxel grid as a numpy array.
        """
        mesh = VoxelizationHelper.load_mesh(geom_file)
        mesh.apply_translation(-mesh.bounds[0])
        grid_size = np.array(grid_size)
        voxel_size = np.array([voxel_size] * 3)
        grid_dimensions = grid_size * voxel_size

        # Center the mesh inside the grid volume
        mesh_offset = (grid_dimensions - mesh.extents) / 2.0
        mesh.apply_translation(mesh_offset)
        
        if description_file is not None:
            geom_description = json.load(open(description_file, "r"))
            mesh_name = list(geom_description.keys())[0]
            mesh_rotation = np.array([
                geom_description[mesh_name]["Transform"]["Rotation"]["X"],
                geom_description[mesh_name]["Transform"]["Rotation"]["Y"],
                geom_description[mesh_name]["Transform"]["Rotation"]["Z"]
            ])
            if mesh_rotation[0] != 0.0 or mesh_rotation[1] != 0.0 or mesh_rotation[2] != 0.0:
                rotation_matrix = trimesh.transformations.rotation_matrix(
                    np.radians(mesh_rotation[0]),
                    [1, 0, 0]
                ) @ trimesh.transformations.rotation_matrix(
                    np.radians(mesh_rotation[1]),
                    [0, 1, 0]
                ) @ trimesh.transformations.rotation_matrix(
                    np.radians(mesh_rotation[2]),
                    [0, 0, 1]
                )
                mesh.apply_transform(rotation_matrix)

            mesh_translation = np.array([
                geom_description[mesh_name]["Transform"]["Translation"]["X"],
                geom_description[mesh_name]["Transform"]["Translation"]["Y"],
                geom_description[mesh_name]["Transform"]["Translation"]["Z"]
            ])
            mesh.apply_translation(mesh_translation)
        x = np.arange(grid_size[0]) * grid_dimensions[0] + grid_dimensions[0] / 2.0
        y = np.arange(grid_size[1]) * grid_dimensions[1] + grid_dimensions[1] / 2.0
        z = np.arange(grid_size[2]) * grid_dimensions[2] + grid_dimensions[2] / 2.0
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        points = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

        # Check which points are inside the mesh
        contained = mesh.contains(points)
        final_grid = contained.reshape(grid_size)

        return final_grid
