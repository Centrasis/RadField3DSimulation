from typing import Sequence

from RadFiled3D.RadFiled3D import FieldStore, vec3, DType


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
