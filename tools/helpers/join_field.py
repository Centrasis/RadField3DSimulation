from typing import Iterable

import numpy as np

from radfiled3d import CartesianRadiationField, DType, VMFMixtureVoxel
from radfiled3d.glm import vec3
from radfiled3d.store import FieldStore

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


VMF_JOIN_MODES = ("beam_lobe", "scatter_only")
VMF_LAYER = "vmf_lobes"
_VALUES_PER_LOBE = 5


def beam_lobe_kappa(distance, voxel_edge: float):
    """Concentration of the direct beam's directions in a voxel at ``distance`` from the focal spot.

    Every primary photon comes from the (point) focal spot, so the directions through one voxel spread over its
    apparent size, about ``voxel_edge / distance`` per axis, uniformly: variance ``(voxel_edge / distance)^2 / 12``.
    A narrow vMF has variance ``1 / kappa`` per axis, hence ``kappa = 12 distance^2 / voxel_edge^2``.
    """
    return 12.0 * np.square(distance) / (voxel_edge * voxel_edge)


def _sort_lobes(lobes: np.ndarray) -> np.ndarray:
    """Canonical slots, as the simulation stores them: strongest lobe first, unused (weight 0) slots all zero."""
    order = np.argsort(-lobes[..., 0], axis=-1, kind="stable")
    lobes = np.take_along_axis(lobes, order[..., None], axis=-2)
    lobes[lobes[..., 0] <= 0] = 0.0
    return lobes


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
    joined_channel: str = "radiation",
    vmf_join_mode: str = "beam_lobe",
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

    A ``vmf_lobes`` layer of the scatter channel (the direct beam has none) is joined by ``vmf_join_mode``:
    ``"beam_lobe"`` adds a lobe for the direct beam, pointing from the tube's focal spot through each voxel
    (``beam_lobe_kappa`` wide), and merges it with the scatter lobes by ``VMFMixtureVoxel.merge``, each weighted by its
    share of the voxel flux. The joined layer has one lobe more, so no lobes are merged into each other;
    ``"scatter_only"`` keeps the scatter lobes unchanged, so they describe only the scattered part.

    :param path: File path to the stored radiation field.
    :param direct_channel: Name of the direct-beam channel to consume.
    :param scatter_channel: Name of the scatter channel to consume.
    :param joined_channel: Name of the combined channel to create.
    :param vmf_join_mode: How a ``vmf_lobes`` layer is joined, one of ``VMF_JOIN_MODES``.
    """
    if vmf_join_mode not in VMF_JOIN_MODES:
        raise ValueError(f"vmf_join_mode must be one of {VMF_JOIN_MODES}, got '{vmf_join_mode}'")
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

    # --- vMF lobes: a mixture is not additive, the scatter lobes are merged with the beam's direction ---
    if VMF_LAYER in beam_layers:
        raise ValueError(f"{path}: the direct-beam channel carries a '{VMF_LAYER}' layer; only scatter lobes can be joined")
    if VMF_LAYER in scatter_layers:
        scatter_lobes = np.asarray(scatter.get_layer_as_ndarray(VMF_LAYER))
        n_lobes = scatter_lobes.shape[-2]
        unit = scatter.get_layer_unit(VMF_LAYER)
        if vmf_join_mode == "beam_lobe":
            origin = FieldStore.load_metadata(path).simulation.tube.radiation_origin
            focal_spot = np.array([origin.x, origin.y, origin.z])
            counts = (vc.x, vc.y, vc.z)
            size = np.array([vc.x * vd.x, vc.y * vd.y, vc.z * vd.z])
            edge = np.array([vd.x, vd.y, vd.z])
            # flat voxel indices of radfiled3d run x fastest (Fortran order)
            flat_beam = beam_flux.reshape(counts).ravel(order="F")
            flat_scatter = scatter_flux.reshape(counts).ravel(order="F")
            lit = np.flatnonzero(flat_beam + flat_scatter > 0)
            centres = (np.stack(np.unravel_index(lit, counts, order="F"), axis=-1) + 0.5) * edge - size / 2.0
            offsets = centres - focal_spot
            distances = np.linalg.norm(offsets, axis=-1)
            means = offsets / np.where(distances > 0, distances, 1.0)[:, None]
            kappas = beam_lobe_kappa(distances, float(edge.max()))

            # VMFMixtureVoxel has no Python constructor: the beam's one-lobe mixture lives in a one-voxel scratch field
            beam_mixture_channel = CartesianRadiationField(vec3(vd.x, vd.y, vd.z), vec3(vd.x, vd.y, vd.z)).add_channel("beam")
            beam_mixture_channel.add_vmf_layer(VMF_LAYER, 1, unit)
            beam_mixture = beam_mixture_channel.get_voxel_flat(VMF_LAYER, 0)
            n_lobes += 1
            dst.add_vmf_layer(VMF_LAYER, n_lobes, unit)
            # the scatter and the beam mixture, each weighted by its share of the voxel flux; the joined layer has room
            # for all lobes, so none are merged into each other
            for i, mean, kappa in zip(lit, means, kappas):
                i = int(i)
                beam_mixture.set_lobe(0, 1.0, mean.tolist(), float(kappa))
                VMFMixtureVoxel.merge(
                    scatter.get_voxel_flat(VMF_LAYER, i), float(flat_scatter[i]),
                    beam_mixture, float(flat_beam[i]),
                    dst.get_voxel_flat(VMF_LAYER, i),
                )
            joined = dst.get_layer_as_ndarray(VMF_LAYER)
            joined[...] = _sort_lobes(np.array(joined, dtype=np.float64)).astype(joined.dtype)
        else:
            dst.add_vmf_layer(VMF_LAYER, n_lobes, unit)
            dst.get_layer_as_ndarray(VMF_LAYER)[...] = scatter_lobes
        dst.set_statistical_error(VMF_LAYER, scatter.get_statistical_error(VMF_LAYER))
        handled.add(VMF_LAYER)

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
