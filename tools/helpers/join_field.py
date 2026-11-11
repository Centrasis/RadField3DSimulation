import os
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


VMF_JOIN_MODES = ("beam_lobe", "beam_lobe_merged", "scatter_only")
VMF_LAYER = "vmf_lobes"
_VALUES_PER_LOBE = 5
# A beam lobe fused with a forward-scatter lobe by the store's merge is stored this wide at most (the store's cap)
# and up to about 0.8 degrees off the focal-spot ray; see reduce_joined_vmf_lobes.
BEAM_KAPPA_STORED_MAX = 1.0e4
BEAM_COS_TOLERANCE = 1.0e-4


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


def _voxel_geometry(counts, voxel_dims, flat_indices: np.ndarray, focal_spot: np.ndarray):
    """Directions from the focal spot through the centres of the voxels at ``flat_indices`` (radfiled3d's flat order,
    x fastest) and the beam-lobe concentration there."""
    size = np.array(counts) * voxel_dims
    centres = (np.stack(np.unravel_index(flat_indices, counts, order="F"), axis=-1) + 0.5) * voxel_dims - size / 2.0
    offsets = centres - focal_spot
    distances = np.linalg.norm(offsets, axis=-1)
    means = offsets / np.where(distances > 0, distances, 1.0)[:, None]
    return means, beam_lobe_kappa(distances, float(np.max(voxel_dims)))


def _scratch_channel(counts, voxel_dims):
    """A channel to hold VMFMixtureVoxel buffers, which have no Python constructor."""
    size = np.array(counts) * voxel_dims
    return CartesianRadiationField(vec3(*size), vec3(*voxel_dims)).add_channel("vmf_scratch")


def _add_beam_lobes(scatter_channel, scatter_layer: str, flat_beam: np.ndarray, flat_scatter: np.ndarray,
                    counts, voxel_dims: np.ndarray, focal_spot: np.ndarray, dst_channel, unit: str,
                    scatter_lobes: int) -> None:
    """Fills the ``VMF_LAYER`` of ``dst_channel`` (``scatter_lobes + 1`` slots): per voxel, the scatter mixture of
    ``scatter_channel``/``scatter_layer``, reduced to ``scatter_lobes`` lobes if it has more (``VMFMixtureVoxel.merge``
    merges the two lobes with the smallest angle between them, repeatedly), merged with one narrow lobe for the direct
    beam from the focal spot through the voxel centre, each weighted by its share of the voxel flux. Sorted by weight."""
    lit = np.flatnonzero(flat_beam + flat_scatter > 0)
    means, kappas = _voxel_geometry(counts, voxel_dims, lit, focal_spot)
    one_voxel = CartesianRadiationField(vec3(*voxel_dims), vec3(*voxel_dims)).add_channel("beam")
    one_voxel.add_vmf_layer(VMF_LAYER, 1, unit)
    beam_mixture = one_voxel.get_voxel_flat(VMF_LAYER, 0)
    source_lobes = scatter_channel.get_voxel_flat(scatter_layer, 0).get_lobes()
    if source_lobes > scatter_lobes:
        reduced_channel = _scratch_channel(counts, voxel_dims)
        reduced_channel.add_vmf_layer("reduced", scatter_lobes, unit)
    dst_channel.add_vmf_layer(VMF_LAYER, scatter_lobes + 1, unit)
    for i, mean, kappa in zip(lit, means, kappas):
        i = int(i)
        scatter_mixture = scatter_channel.get_voxel_flat(scatter_layer, i)
        if source_lobes > scatter_lobes:
            reduced = reduced_channel.get_voxel_flat("reduced", i)
            VMFMixtureVoxel.merge(scatter_mixture, 1.0, scatter_mixture, 0.0, reduced)
            scatter_mixture = reduced
        beam_mixture.set_lobe(0, 1.0, mean.tolist(), float(kappa))
        VMFMixtureVoxel.merge(scatter_mixture, float(flat_scatter[i]), beam_mixture, float(flat_beam[i]), dst_channel.get_voxel_flat(VMF_LAYER, i))
    joined = dst_channel.get_layer_as_ndarray(VMF_LAYER)
    joined[...] = _sort_lobes(np.array(joined, dtype=np.float64)).astype(joined.dtype)


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
    error is the relative error of the summed flux (absolute errors added in quadrature). Any further per-voxel layer common to both channels (e.g. an
    angular layer) is summed per primary. Channels other than the two beam channels (e.g. a geometry
    channel) are copied unchanged. The result is stored back to ``path``, overriding it. If either
    beam channel is missing, nothing is done.

    A ``vmf_lobes`` layer of the scatter channel (the direct beam has none) is joined by ``vmf_join_mode``:
    ``"beam_lobe"`` adds a lobe for the direct beam, pointing from the tube's focal spot through each voxel
    (``beam_lobe_kappa`` wide), and merges it with the scatter lobes by ``VMFMixtureVoxel.merge``, each weighted by its
    share of the voxel flux. The joined layer has one lobe more, so no lobes are merged into each other;
    ``"beam_lobe_merged"`` does the same but first reduces the scatter lobes by one (merging the two with the smallest
    angle between them), so the joined layer keeps the scatter layer's lobe count and the beam lobe stays intact;
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

    # --- error: relative standard errors of the flux; the absolute errors of the two independent parts add in
    # quadrature, sigma = sqrt((R_beam * flux_beam)^2 + (R_scatter * flux_scatter)^2), R = sigma / total flux ---
    if "error" in beam_layers and "error" in scatter_layers:
        sigma = np.hypot(np.asarray(beam.get_layer_as_ndarray("error")) * beam_flux, np.asarray(scatter.get_layer_as_ndarray("error")) * scatter_flux)
        with np.errstate(divide="ignore", invalid="ignore"):
            error = np.where(total_flux > 0, sigma / total_flux, 1.0)
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
        if vmf_join_mode == "scatter_only":
            dst.add_vmf_layer(VMF_LAYER, n_lobes, unit)
            dst.get_layer_as_ndarray(VMF_LAYER)[...] = scatter_lobes
        else:
            scatter_slots = n_lobes - 1 if vmf_join_mode == "beam_lobe_merged" else n_lobes
            if scatter_slots < 1:
                raise ValueError(f"{path}: '{vmf_join_mode}' needs at least 2 scatter lobes, the field has {n_lobes}")
            origin = FieldStore.load_metadata(path).simulation.tube.radiation_origin
            counts = (vc.x, vc.y, vc.z)
            _add_beam_lobes(
                scatter, VMF_LAYER,
                beam_flux.reshape(counts).ravel(order="F"), scatter_flux.reshape(counts).ravel(order="F"),
                counts, np.array([vd.x, vd.y, vd.z]), np.array([origin.x, origin.y, origin.z]), dst, unit, scatter_slots,
            )
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


def reduce_joined_vmf_lobes(path: str, lobes: int, channel: str = "radiation") -> bool:
    """Reduces the ``vmf_lobes`` layer of an already joined field (``"beam_lobe"``, scatter lobes + 1 beam lobe) to
    ``lobes`` slots, giving the same result as joining with ``"beam_lobe_merged"`` would have.

    The beam lobe of each voxel is recognized by its concentration (``beam_lobe_kappa``) and direction (from the focal
    spot through the voxel centre) and kept; the other lobes are reduced to ``lobes - 1`` by ``VMFMixtureVoxel.merge``
    and merged back with the beam lobe by their weights. The field is replaced through a temporary file.

    :return: False if the layer has no more than ``lobes`` slots already (nothing done), True otherwise.
    """
    field = FieldStore.load(path)
    metadata = FieldStore.load_metadata(path)
    src = field.get_channel(channel)
    current = np.array(src.get_layer_as_ndarray(VMF_LAYER), dtype=np.float64)
    slots = current.shape[-2]
    if slots <= lobes:
        return False
    if lobes < 2:
        raise ValueError("at least 2 lobes are needed: one for the beam and one for the scatter")
    unit = src.get_layer_unit(VMF_LAYER)
    vd = field.get_voxel_dimensions(); vc = field.get_voxel_counts()
    counts = (vc.x, vc.y, vc.z)
    voxel_dims = np.array([vd.x, vd.y, vd.z])
    origin = metadata.simulation.tube.radiation_origin
    focal_spot = np.array([origin.x, origin.y, origin.z])

    flat = current.reshape(counts + current.shape[-2:]).reshape((-1,) + current.shape[-2:], order="F")
    lit = np.flatnonzero(flat[:, :, 0].sum(axis=-1) > 0)
    means, kappas = _voxel_geometry(counts, voxel_dims, lit, focal_spot)
    lobes_lit = flat[lit]
    # The join wrote the beam lobe with beam_lobe_kappa along the focal-spot ray, but the store fuses lobes of almost
    # the same direction when it merges mixtures: where the scatter mixture holds a forward-scatter lobe along the
    # ray, the stored beam lobe is that fusion, a few 1e-3 off the ray and at most BEAM_KAPPA_STORED_MAX wide. The
    # beam lobe is therefore the lobe within BEAM_COS_TOLERANCE of the ray that still has the join's concentration,
    # or else the most concentrated one there, provided it is at least half as narrow as the join wrote it (a
    # forward-scatter lobe beside an intact beam lobe can be narrower than the beam lobe, so exactness ranks first).
    # A voxel whose lobe was fused below that keeps its lobes as scatter.
    along_ray = (lobes_lit[:, :, 0] > 0) & (np.einsum("vkc,vc->vk", lobes_lit[:, :, 1:4], means) > 1.0 - BEAM_COS_TOLERANCE)
    exact = np.isclose(lobes_lit[:, :, 4], kappas[:, None], rtol=1e-4)
    rank = np.where(along_ray, np.where(exact, np.inf, lobes_lit[:, :, 4]), -np.inf)
    expected = np.minimum(kappas, BEAM_KAPPA_STORED_MAX)
    is_beam = np.zeros(along_ray.shape, dtype=bool)
    is_beam[np.arange(len(lit)), rank.argmax(axis=-1)] = along_ray.any(axis=-1) & (rank.max(axis=-1) >= 0.5 * expected)
    beam_weight = np.where(is_beam, lobes_lit[:, :, 0], 0.0).sum(axis=-1)

    # the lobes without the beam lobe, in a scratch layer the merge can read; a voxel without a recognized beam lobe
    # hands over all its slots (the empty beam slot is then an all-zero lobe, which the merge ignores)
    scratch = _scratch_channel(counts, voxel_dims)
    scratch.add_vmf_layer("scatter", slots, unit)
    scatter_only = np.zeros((len(lit), slots, 5))
    for v in range(len(lit)):
        kept = lobes_lit[v][~is_beam[v]]
        scatter_only[v, : len(kept)] = kept
    scatter_flat = np.zeros((flat.shape[0], slots, 5))
    scatter_flat[lit] = scatter_only
    scratch.get_layer_as_ndarray("scatter")[...] = scatter_flat.reshape(counts + (slots, 5), order="F").astype(np.float32)
    scatter_weight = np.zeros(flat.shape[0])
    scatter_weight[lit] = scatter_only[:, :, 0].sum(axis=-1)
    beam_weight_flat = np.zeros(flat.shape[0])
    beam_weight_flat[lit] = beam_weight

    def fill(dst):
        # the beam lobes carry exactly the beam's flux share, so these weights reproduce the original join
        _add_beam_lobes(scratch, "scatter", beam_weight_flat, scatter_weight, counts, voxel_dims, focal_spot, dst, unit, lobes - 1)
    _replace_vmf_layer(path, field, metadata, channel, fill)
    return True


def _replace_vmf_layer(path: str, field, metadata, channel: str, fill) -> None:
    """Stores ``field`` with its ``channel``/``vmf_lobes`` layer rebuilt by ``fill(dst_channel)`` (which must add the
    layer), all other layers copied, replacing ``path`` through a temporary file."""
    vd = field.get_voxel_dimensions(); vc = field.get_voxel_counts()
    out = CartesianRadiationField(vec3(vc.x * vd.x, vc.y * vd.y, vc.z * vd.z), vec3(vd.x, vd.y, vd.z))
    for channel_name in field.get_channel_names():
        source = field.get_channel(channel_name)
        target = out.add_channel(channel_name)
        _copy_channel(source, target, [l for l in source.get_layers() if not (channel_name == channel and l == VMF_LAYER)])
    src = field.get_channel(channel)
    dst = out.get_channel(channel)
    fill(dst)
    dst.set_statistical_error(VMF_LAYER, src.get_statistical_error(VMF_LAYER))
    tmp = path + ".reducing"
    FieldStore.store(out, metadata, tmp)
    os.replace(tmp, path)


def reduce_scatter_vmf_lobes(path: str, lobes: int, channel: str = "scatter_field") -> bool:
    """Reduces a ``vmf_lobes`` layer holding scatter lobes only (a field whose channels were not joined) to ``lobes``
    slots: per voxel ``VMFMixtureVoxel.merge`` merges the two lobes with the smallest angle between them until
    ``lobes`` are left; sorted by weight. The field is replaced through a temporary file.

    :return: False if the layer has no more than ``lobes`` slots already (nothing done), True otherwise.
    """
    if lobes < 1:
        raise ValueError("at least 1 lobe is needed")
    field = FieldStore.load(path)
    metadata = FieldStore.load_metadata(path)
    src = field.get_channel(channel)
    slots = src.get_voxel_flat(VMF_LAYER, 0).get_lobes()
    if slots <= lobes:
        return False
    unit = src.get_layer_unit(VMF_LAYER)
    weights = np.asarray(src.get_layer_as_ndarray(VMF_LAYER))[..., 0].sum(axis=-1)
    lit = np.flatnonzero(weights.ravel(order="F") > 0)

    def fill(dst):
        dst.add_vmf_layer(VMF_LAYER, lobes, unit)
        for i in lit:
            mixture = src.get_voxel_flat(VMF_LAYER, int(i))
            VMFMixtureVoxel.merge(mixture, 1.0, mixture, 0.0, dst.get_voxel_flat(VMF_LAYER, int(i)))
        reduced = dst.get_layer_as_ndarray(VMF_LAYER)
        reduced[...] = _sort_lobes(np.array(reduced, dtype=np.float64)).astype(reduced.dtype)
    _replace_vmf_layer(path, field, metadata, channel, fill)
    return True


def reduce_vmf_lobes(path: str, lobes: int, joined_channel: str = "radiation", scatter_channel: str = "scatter_field") -> str:
    """Reduces the ``vmf_lobes`` of a field to ``lobes`` slots if it has more: a joined field by
    ``reduce_joined_vmf_lobes`` (keeps the beam lobe), an unjoined one by ``reduce_scatter_vmf_lobes``.

    :return: ``"reduced"``, ``"unchanged"`` (no more than ``lobes`` slots) or ``"no_lobes"`` (no ``vmf_lobes`` layer).
    """
    channels = FieldStore.load(path).get_channel_names()
    for channel, reduce in ((joined_channel, reduce_joined_vmf_lobes), (scatter_channel, reduce_scatter_vmf_lobes)):
        if channel in channels and VMF_LAYER in FieldStore.load(path).get_channel(channel).get_layers():
            return "reduced" if reduce(path, lobes, channel) else "unchanged"
    return "no_lobes"


UNIFORM_LOBE = (1.0, 0.0, 0.0, 1.0, 0.0)   # weight 1, any unit mean, kappa 0: the uniform distribution
RELATIVE_ERROR_UNIT = "relative std. error of flux"
# a valid relative-error layer of a 1e8-photon field has a median of about 0.02 over the lit voxels; the layers written
# before the error was fixed sit at about 0.5 everywhere
_MAX_PLAUSIBLE_MEDIAN_ERROR = 0.3


def fill_uniform_vmf_lobes(path: str, min_samples: float = None, joined_channel: str = "radiation",
                           scatter_channel: str = "scatter_field") -> dict:
    """Gives an existing field the lobes the current VMFTrainer stores for voxels without enough samples: the uniform
    distribution (weight 1, kappa 0).

    - Exact: a voxel with flux but no lobes gets the uniform lobe (in an unjoined field every voxel without lobes, as
      the trainer starts every voxel uniform; a joined field keeps voxels without flux empty, as the join does).
    - Approximate, only with ``min_samples``: a joined voxel without a beam lobe (all its lobes are scatter) whose
      relative flux error R gives fewer than ``min_samples`` effective samples (n = 1 / R^2) becomes uniform too. The
      error covers the whole run while the trainer counts its last two passes, so this keeps somewhat more fitted
      voxels than the trainer would. Needs a valid relative-error layer; a field without one raises ``ValueError``.

    The field is replaced through a temporary file only if something changed.

    :return: ``{"channel": ..., "empty_filled": n, "sparse_filled": n}``, or ``{"channel": None, ...}`` without lobes.
    """
    field = FieldStore.load(path)
    channels = field.get_channel_names()
    joined = joined_channel in channels and VMF_LAYER in field.get_channel(joined_channel).get_layers()
    channel = joined_channel if joined else scatter_channel
    result = {"channel": None, "empty_filled": 0, "sparse_filled": 0}
    if channel not in channels or VMF_LAYER not in field.get_channel(channel).get_layers():
        return result
    result["channel"] = channel
    src = field.get_channel(channel)
    lobes = np.asarray(src.get_layer_as_ndarray(VMF_LAYER))
    weights = lobes[..., 0].sum(axis=-1)
    uniform = np.zeros(lobes.shape[-2:], dtype=lobes.dtype)
    uniform[0] = UNIFORM_LOBE

    if joined:
        flux = np.asarray(src.get_layer_as_ndarray("flux"))[..., 0]
        empty = (flux > 0) & (weights <= 0)
    else:
        empty = weights <= 0
    sparse = np.zeros_like(empty)

    if min_samples is not None and joined:
        if "error" not in src.get_layers() or src.get_layer_unit("error") != RELATIVE_ERROR_UNIT:
            raise ValueError(f"{path}: --min-samples needs a '{RELATIVE_ERROR_UNIT}' error layer")
        error = np.asarray(src.get_layer_as_ndarray("error"), dtype=np.float64)[..., 0]
        lit = flux > 0
        if np.median(error[lit]) > _MAX_PLAUSIBLE_MEDIAN_ERROR:
            raise ValueError(f"{path}: the error layer is not a valid relative error (median {np.median(error[lit]):.2f} over the lit voxels)")
        metadata = FieldStore.load_metadata(path)
        vd = field.get_voxel_dimensions(); vc = field.get_voxel_counts()
        counts = (vc.x, vc.y, vc.z)
        origin = metadata.simulation.tube.radiation_origin
        flat_lobes = lobes.reshape((-1,) + lobes.shape[-2:], order="F")
        candidates = np.flatnonzero((weights > 0).ravel(order="F"))
        means, kappas = _voxel_geometry(counts, np.array([vd.x, vd.y, vd.z]), candidates, np.array([origin.x, origin.y, origin.z]))
        cand = flat_lobes[candidates].astype(np.float64)
        has_beam = (
            (cand[:, :, 0] > 0)
            & np.isclose(cand[:, :, 4], kappas[:, None], rtol=1e-4)
            & (np.einsum("vkc,vc->vk", cand[:, :, 1:4], means) > 1.0 - 1e-6)
        ).any(axis=-1)
        flat_error = error.ravel(order="F")[candidates]
        with np.errstate(divide="ignore"):
            effective = np.where(flat_error > 0, 1.0 / np.square(flat_error), np.inf)
        already_uniform = (cand[:, 0, 0] == 1.0) & (cand[:, 0, 4] == 0.0)
        chosen = candidates[~has_beam & (effective < min_samples) & ~already_uniform]
        flat_sparse = np.zeros(int(np.prod(counts)), dtype=bool)
        flat_sparse[chosen] = True
        sparse = flat_sparse.reshape(counts, order="F")

    change = empty | sparse
    result["empty_filled"] = int(empty.sum())
    result["sparse_filled"] = int(sparse.sum())
    if not change.any():
        return result
    lobes[change] = uniform
    tmp = path + ".filling"
    FieldStore.store(field, FieldStore.load_metadata(path), tmp)
    os.replace(tmp, path)
    return result
