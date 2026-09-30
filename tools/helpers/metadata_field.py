from typing import Sequence

import numpy as np

from radfiled3d import DType
from radfiled3d.glm import vec3
from radfiled3d.store import FieldStore


def add_patient_translation(rf3_path: str, translation: Sequence[float]) -> None:
    """Store the patient translation as a vec3 dynamic-metadata entry on the field at ``rf3_path``.

    The translation ``(x, y, z)`` is the value the patient mesh was moved by (metres) and is written
    under the key ``patient_translation``. The field is loaded, the metadata is extended and the field
    is stored back, overriding the file.

    :param rf3_path: File path to the stored radiation field.
    :param translation: The patient translation as (x, y, z) in metres.
    """
    field = FieldStore.load(rf3_path)
    metadata = FieldStore.load_metadata(rf3_path)
    if "patient_translation" in metadata.get_dynamic_metadata_keys():
        vx = metadata.get_dynamic_metadata("patient_translation")
    else:
        vx = metadata.add_dynamic_metadata("patient_translation", DType.VEC3)
    vx.set_data(vec3(float(translation[0]), float(translation[1]), float(translation[2])))
    FieldStore.store(field, metadata, rf3_path)


PATIENT_OVERLAP_KEY = "patient_overlap"
C_ARM_LAYERS = ("imagedetector", "xraytube")


def mark_patient_overlap(rf3_path: str) -> bool:
    """Flag a field whose C-arm collides with the patient: an ``imagedetector`` or ``xraytube`` geometry layer shares a
    voxel with the ``patient`` layer, or the tube's focal spot lies in a patient voxel.

    Only overlapping fields get the byte dynamic metadata ``patient_overlap`` (value 1); others stay unchanged, so the
    key's presence is the filter. The check works on the stored geometry channel, so a detector within one voxel of
    the patient counts as overlapping too.

    :param rf3_path: File path to the stored radiation field.
    :return: Whether the field overlaps.
    """
    field = FieldStore.load(rf3_path)
    if not field.has_channel("geometry"):
        return False
    geometry = field.get_channel("geometry")
    layers = geometry.get_layers()
    if "patient" not in layers:
        return False
    patient = np.squeeze(geometry.get_layer_as_ndarray("patient"), axis=-1) > 0

    overlap = any(
        np.any(patient & (np.squeeze(geometry.get_layer_as_ndarray(layer), axis=-1) > 0))
        for layer in C_ARM_LAYERS if layer in layers
    )

    metadata = FieldStore.load_metadata(rf3_path)
    if not overlap:
        # voxel centres lie at (index + 0.5) * voxel size - field size / 2
        origin = metadata.simulation.tube.radiation_origin
        counts = geometry.get_voxel_counts()
        voxel = field.get_voxel_dimensions()
        size = field.get_field_dimensions()
        index = [
            int(np.floor((coord + extent / 2.0) / edge))
            for coord, extent, edge in zip((origin.x, origin.y, origin.z), (size.x, size.y, size.z), (voxel.x, voxel.y, voxel.z))
        ]
        if all(0 <= i < n for i, n in zip(index, (counts.x, counts.y, counts.z))):
            overlap = bool(patient[index[0], index[1], index[2]])

    if overlap:
        if PATIENT_OVERLAP_KEY in metadata.get_dynamic_metadata_keys():
            vx = metadata.get_dynamic_metadata(PATIENT_OVERLAP_KEY)
        else:
            vx = metadata.add_dynamic_metadata(PATIENT_OVERLAP_KEY, DType.BYTE)
        vx.set_data(1)
        FieldStore.store(field, metadata, rf3_path)
    return overlap
