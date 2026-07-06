import trimesh
import rtree # needed to be included for voxelization
import numpy as np
import json
from logging import Logger
from typing import Union, cast


class VoxelizationHelper:
    @staticmethod
    def load_meshes(geom_file: str) -> dict[str, trimesh.Trimesh]:
        """
        Load all meshes from a geometry file, keyed by their scene object name.

        :param geom_file: Path to the geometry file.
        :return: Dict of object name -> trimesh object.
        """
        loaded = trimesh.load(geom_file)
        if isinstance(loaded, trimesh.Scene):
            return {
                name: trimesh.Trimesh(mesh.vertices, mesh.faces, face_normals=mesh.face_normals, vertex_normals=mesh.vertex_normals)
                for name, mesh in loaded.geometry.items()
            }
        return {"mesh": loaded}

    @staticmethod
    def load_mesh(geom_file: str) -> trimesh.Trimesh:
        """
        Load a geometry file as a single concatenated mesh.

        :param geom_file: Path to the geometry file.
        :return: Loaded mesh as a trimesh object.
        """
        meshes = list(VoxelizationHelper.load_meshes(geom_file).values())
        return meshes[0] if len(meshes) == 1 else trimesh.util.concatenate(meshes)

    @staticmethod
    def _apply_transform(mesh: trimesh.Trimesh, transform: dict) -> trimesh.Trimesh:
        """Apply a desc Transform to a mesh copy in the AUTHORED frame, exactly like the
        simulation places it (GeometryLoader/G4Mesh): scale about the authored origin, then the
        X/Y/Z rotation about the authored origin, then the translation. The world origin is the
        world center, matching the detector's coordinate mapping."""
        mesh = mesh.copy()
        scale = transform.get("Scale")
        if scale is not None:
            mesh.vertices = mesh.vertices * np.array([scale["X"], scale["Y"], scale["Z"]])
        rotation = transform.get("Rotation")
        if rotation is not None and any(rotation[a] != 0.0 for a in ("X", "Y", "Z")):
            matrix = trimesh.transformations.rotation_matrix(np.radians(rotation["X"]), [1, 0, 0]) \
                   @ trimesh.transformations.rotation_matrix(np.radians(rotation["Y"]), [0, 1, 0]) \
                   @ trimesh.transformations.rotation_matrix(np.radians(rotation["Z"]), [0, 0, 1])
            mesh.apply_transform(matrix)
        translation = transform.get("Translation")
        if translation is not None:
            mesh.apply_translation([translation["X"], translation["Y"], translation["Z"]])
        return mesh

    @staticmethod
    def _rasterize_into(grid: np.ndarray, mesh: trimesh.Trimesh, voxel_size: np.ndarray, world_center_m: np.ndarray):
        """Mark every grid voxel whose CENTER lies inside the mesh. The mesh is in world
        coordinates (origin = world center); the grid spans [0, grid_size*voxel_size) with the
        world center at ``world_center_m``. Candidate voxels are limited to the mesh bounding box,
        so memory stays proportional to the object, not the field."""
        lo = np.maximum(np.floor((mesh.bounds[0] + world_center_m) / voxel_size).astype(int), 0)
        hi = np.minimum(np.ceil((mesh.bounds[1] + world_center_m) / voxel_size).astype(int), np.array(grid.shape))
        if np.any(hi <= lo):
            return
        axes = [np.arange(lo[i], hi[i]) for i in range(3)]
        xx, yy, zz = np.meshgrid(*axes, indexing="ij")
        centers = (np.stack([xx, yy, zz], axis=-1).reshape(-1, 3) + 0.5) * voxel_size - world_center_m
        inside = mesh.contains(centers).reshape(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2])
        grid[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] |= inside

    @staticmethod
    def _walk_desc(desc: dict, parent_transform: dict = None):
        """Yield (object_name, effective_transform) for every desc node, composing a parent's
        transform onto its Children the way the simulation nests placements."""
        for name, node in desc.items():
            if not isinstance(node, dict):
                continue
            transform = node.get("Transform", {})
            if parent_transform:
                transform = {**parent_transform, **transform}
            yield name, transform
            children = node.get("Children")
            if children:
                yield from VoxelizationHelper._walk_desc(children, transform)

    @staticmethod
    def generate_voxelgrid_with_geometry(geom_file: str, voxel_size: Union[float, np.ndarray], grid_size: Union[np.ndarray, tuple[int, int, int]], description_file: Union[None, str] = None, logger: Union[None, Logger] = None, world_center_m: Union[None, np.ndarray] = None) -> np.ndarray:
        """
        Generate a voxel grid from a geometry file, placing every object exactly like the
        simulation does: authored coordinates + desc Scale/Rotation/Translation, world origin at
        the world center. A voxel is set when its center lies inside a mesh.

        :param geom_file: Path to the geometry file.
        :param voxel_size: Size of each voxel.
        :param grid_size: Size of the grid as a tuple of integers.
        :param description_file: Optional path to a description file
        :param logger: Optional logger to log information during processing
        :param world_center_m: Grid-frame position of the world center in meters. Defaults to the
            grid center; pass the shifted center when the field was cropped asymmetrically.
        :return: Voxel grid as a numpy array.
        """
        meshes = VoxelizationHelper.load_meshes(geom_file)
        if logger is not None:
            logger.debug(f"generate_voxelgrid_with_geometry(): Loaded {len(meshes)} mesh(es) from {geom_file}")

        grid_size = np.array(grid_size) if not isinstance(grid_size, np.ndarray) else grid_size
        voxel_size: np.ndarray = cast(np.ndarray, np.array([voxel_size] * 3) if isinstance(voxel_size, float) else voxel_size)
        world_center_m = (grid_size * voxel_size) / 2.0 if world_center_m is None else np.asarray(world_center_m, dtype=np.float64)

        transforms = {}
        if description_file is not None:
            desc = json.load(open(description_file, "r"))
            transforms = dict(VoxelizationHelper._walk_desc(desc))

        final_grid = np.zeros(tuple(grid_size), dtype=bool)
        matched = [name for name in meshes if name in transforms]
        if transforms and not matched and len(meshes) == 1:
            # single unnamed mesh + desc: apply the first desc object's transform
            only_mesh = next(iter(meshes.values()))
            first_transform = next(iter(transforms.values()))
            VoxelizationHelper._rasterize_into(final_grid, VoxelizationHelper._apply_transform(only_mesh, first_transform), voxel_size, world_center_m)
        else:
            for name, mesh in meshes.items():
                placed = VoxelizationHelper._apply_transform(mesh, transforms.get(name, {}))
                VoxelizationHelper._rasterize_into(final_grid, placed, voxel_size, world_center_m)

        if logger is not None:
            logger.debug("generate_voxelgrid_with_geometry(): Voxelization completed!")

        return final_grid
