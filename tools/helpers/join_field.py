from typing import Iterable

import numpy as np

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


def _add_layer_like(src_channel, dst_channel, layer_name: str) -> None:
    """Create ``layer_name`` on ``dst_channel`` with the same type/unit as on ``src_channel``."""
    voxel_type = src_channel.get_layer_voxel_type(layer_name)
    unit = src_channel.get_layer_unit(layer_name)
    if voxel_type == "histogram":
        s = src_channel.get_voxel_flat(layer_name, 0)
        dst_channel.add_histogram_layer(layer_name, s.get_bins(), s.get_histogram_bin_width(), unit)
    elif voxel_type == "spherical":
        s = src_channel.get_voxel_flat(layer_name, 0)
        dst_channel.add_spherical_layer(layer_name, s.get_phi_segments(), s.get_theta_segments(), unit)
    else:
        try:
            dtype = _SCALAR_DTYPE_BY_NAME[voxel_type]
        except KeyError:
            raise ValueError(f"unsupported voxel type '{voxel_type}' in layer '{layer_name}'")
        dst_channel.add_layer(layer_name, unit, dtype)


def _copy_channel(src_channel, dst_channel, layer_names: Iterable[str]) -> None:
    """Copy the listed layers from ``src_channel`` to ``dst_channel`` unchanged."""
    for layer_name in layer_names:
        _add_layer_like(src_channel, dst_channel, layer_name)
        dst_channel.get_layer_as_ndarray(layer_name)[...] = src_channel.get_layer_as_ndarray(layer_name)
        dst_channel.set_statistical_error(layer_name, src_channel.get_statistical_error(layer_name))


def join_rf3_file(
    path: str,
    direct_channel: str = "direct_beam",
    scatter_channel: str = "scatter_field",
    joined_channel: str = "joined_beam",
) -> None:
    """Join the direct-beam and scatter channels of the field at ``path`` into one channel.

    Halves the per-field storage by replacing the two beam channels with a single combined one.
    The Monte-Carlo per-primary normalization of the flux is preserved (fluxes add); the spectrum
    is combined as a flux-weighted mix of the two channels' per-voxel distributions, renormalized
    per voxel along the histogram bins and zeroed where both channels are empty; the statistical
    error is the mean of the two. Any further per-voxel layer common to both channels (e.g. an
    angular layer) is summed per primary. Channels other than the two beam channels (e.g. a geometry
    channel) are copied unchanged. The result is stored back to ``path``, overriding it. If either
    beam channel is missing, nothing is done.

    :param path: File path to the stored radiation field.
    :param direct_channel: Name of the direct-beam channel to consume.
    :param scatter_channel: Name of the scatter channel to consume.
    :param joined_channel: Name of the combined channel to create.
    """
    field = FieldStore.load(path)
    if not isinstance(field, CartesianRadiationField):
        raise TypeError("join_rf3_file only supports Cartesian radiation fields")

    channel_names = list(field.get_channel_names())
    if direct_channel not in channel_names or scatter_channel not in channel_names:
        return

    beam = field.get_channel(direct_channel)
    scatter = field.get_channel(scatter_channel)

    vd = field.get_voxel_dimensions()
    vc = field.get_voxel_counts()
    out = CartesianRadiationField(
        vec3(vc.x * vd.x, vc.y * vd.y, vc.z * vd.z),
        vec3(vd.x, vd.y, vd.z),
    )
    dst = out.add_channel(joined_channel)

    beam_layers = set(beam.get_layers())
    scatter_layers = set(scatter.get_layers())

    # --- flux: per-primary, additive (keeps the MC per-primary normalization) ---
    beam_flux = np.asarray(beam.get_layer_as_ndarray("flux"))
    scatter_flux = np.asarray(scatter.get_layer_as_ndarray("flux"))
    total_flux = beam_flux + scatter_flux
    _add_layer_like(beam, dst, "flux")
    dst.get_layer_as_ndarray("flux")[...] = total_flux.astype(dst.get_layer_as_ndarray("flux").dtype)
    dst.set_statistical_error("flux", 0.5 * (beam.get_statistical_error("flux") + scatter.get_statistical_error("flux")))

    handled = {"flux"}

    # --- spectrum: exactly RadField3D-NN ChannelsJoin (flux-weighted mix, renormalized per voxel) ---
    if "spectrum" in beam_layers and "spectrum" in scatter_layers:
        eps = 1e-8
        ratio_beam = (beam_flux + eps) / (total_flux + eps)
        ratio_scatter = (scatter_flux + eps) / (total_flux + eps)
        spectrum = (
            ratio_scatter * np.asarray(scatter.get_layer_as_ndarray("spectrum"))
            + ratio_beam * np.asarray(beam.get_layer_as_ndarray("spectrum"))
        )
        spectrum = spectrum / np.clip(spectrum.sum(axis=-1, keepdims=True), eps, None)
        spectrum = np.where(total_flux <= 0, 0.0, spectrum)
        _add_layer_like(beam, dst, "spectrum")
        dst.get_layer_as_ndarray("spectrum")[...] = spectrum.astype(dst.get_layer_as_ndarray("spectrum").dtype)
        dst.set_statistical_error("spectrum", 0.5 * (beam.get_statistical_error("spectrum") + scatter.get_statistical_error("spectrum")))
        handled.add("spectrum")

    # --- error: mean of the two channels' statistical-error estimates ---
    if "error" in beam_layers and "error" in scatter_layers:
        error = 0.5 * (np.asarray(beam.get_layer_as_ndarray("error")) + np.asarray(scatter.get_layer_as_ndarray("error")))
        _add_layer_like(beam, dst, "error")
        dst.get_layer_as_ndarray("error")[...] = error.astype(dst.get_layer_as_ndarray("error").dtype)
        handled.add("error")

    # --- any other per-voxel layer common to both channels: additive per primary (e.g. angular flux) ---
    for layer_name in beam.get_layers():
        if layer_name in handled or layer_name not in scatter_layers:
            continue
        _add_layer_like(beam, dst, layer_name)
        joined = np.asarray(beam.get_layer_as_ndarray(layer_name)) + np.asarray(scatter.get_layer_as_ndarray(layer_name))
        dst.get_layer_as_ndarray(layer_name)[...] = joined.astype(dst.get_layer_as_ndarray(layer_name).dtype)

    # --- preserve any non-beam channels (e.g. geometry) unchanged ---
    for channel_name in channel_names:
        if channel_name in (direct_channel, scatter_channel):
            continue
        src = field.get_channel(channel_name)
        _copy_channel(src, out.add_channel(channel_name), src.get_layers())

    metadata = FieldStore.load_metadata(path)
    FieldStore.store(out, metadata, path)
