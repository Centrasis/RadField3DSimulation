from typing import Tuple

from RadFiled3D.RadFiled3D import (
    FieldStore,
    CartesianRadiationField,
    DType,
    vec3,
)


_SCALAR_DTYPE_BY_NAME = {
    "float16": DType.FLOAT16,
    "float": DType.FLOAT32,
    "double": DType.FLOAT64,
    "int": DType.INT32,
    "char": DType.SCHAR,
    "unsigned char": DType.BYTE,
    "uint8_t": DType.BYTE,
    "uint64_t": DType.UINT64,
    "unsigned long long": DType.UINT64,
    "unsigned long": DType.UINT64,
    "uint32_t": DType.UINT32,
    "unsigned int": DType.UINT32,
    "glm::vec2": DType.VEC2,
    "glm::vec3": DType.VEC3,
    "glm::vec4": DType.VEC4,
}


def crop_rf3_file(path: str, crop_voxel_counts: Tuple[int, int, int]) -> None:
    """Crop the Cartesian radiation field stored at ``path`` to a smaller voxel count.

    Loads the field, checks the requested voxel counts against the current ones and, if
    they are smaller, builds a new field that keeps the voxel edge length but holds fewer
    voxels, cropping each layer's volume by centering the kept region around the middle.
    All channels and layers are preserved. The cropped field is then stored back to the
    same file, overriding it. If the requested counts equal the current ones, nothing is
    done. Requesting a larger count in any dimension raises a ``ValueError``.

    :param path: File path to the stored radiation field.
    :param crop_voxel_counts: Target voxel counts (x, y, z).
    """
    field = FieldStore.load(path)
    if not isinstance(field, CartesianRadiationField):
        raise TypeError("crop_rf3_file only supports Cartesian radiation fields")

    current = field.get_voxel_counts()
    cur = (current.x, current.y, current.z)
    crop = tuple(int(c) for c in crop_voxel_counts)

    if any(c <= 0 for c in crop):
        raise ValueError(f"crop voxel counts must be positive, got {crop}")
    if any(c > o for c, o in zip(crop, cur)):
        raise ValueError(f"crop voxel counts {crop} exceed the field's voxel counts {cur}")
    if crop == cur:
        return

    cx, cy, cz = crop
    start = tuple((o - c) // 2 for c, o in zip(crop, cur))
    sx, sy, sz = start

    vd = field.get_voxel_dimensions()
    cropped = CartesianRadiationField(
        vec3(cx * vd.x, cy * vd.y, cz * vd.z),
        vec3(vd.x, vd.y, vd.z),
    )

    for channel_name in field.get_channel_names():
        src_channel = field.get_channel(channel_name)
        dst_channel = cropped.add_channel(channel_name)

        for layer_name in src_channel.get_layers():
            voxel_type = src_channel.get_layer_voxel_type(layer_name)
            unit = src_channel.get_layer_unit(layer_name)

            if voxel_type == "histogram":
                sample = src_channel.get_voxel_flat(layer_name, 0)
                dst_channel.add_histogram_layer(layer_name, sample.get_bins(), sample.get_histogram_bin_width(), unit)
            elif voxel_type == "spherical":
                sample = src_channel.get_voxel_flat(layer_name, 0)
                dst_channel.add_spherical_layer(layer_name, sample.get_phi_segments(), sample.get_theta_segments(), unit)
            else:
                try:
                    dtype = _SCALAR_DTYPE_BY_NAME[voxel_type]
                except KeyError:
                    raise ValueError(f"unsupported voxel type '{voxel_type}' in layer '{layer_name}'")
                dst_channel.add_layer(layer_name, unit, dtype)

            src = src_channel.get_layer_as_ndarray(layer_name)
            dst = dst_channel.get_layer_as_ndarray(layer_name)
            if voxel_type == "spherical":
                dst[...] = src[..., sx:sx + cx, sy:sy + cy, sz:sz + cz]
            else:
                dst[...] = src[sx:sx + cx, sy:sy + cy, sz:sz + cz]
            dst_channel.set_statistical_error(layer_name, src_channel.get_statistical_error(layer_name))

    metadata = FieldStore.load_metadata(path)
    FieldStore.store(cropped, metadata, path)
