import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

from create_dataset import GeometrySampler
import numpy as np


desc_file = """
{
    "patient": {
        "MaterialName": "G4_TISSUE_SOFT_ICRU-4",
        "Type": "Patient",
        "Transform": {
            "Rotation": {
                "X": 0,
                "Y": 0.0,
                "Z": 0.0
            },
            "Translation": {
                "X": 0.0,
                "Y": 0.5,
                "Z": 0.0
            },
            "Scale": {
                "X": 1.0,
                "Y": 1.0,
                "Z": 1.0
            }
        },
        "Children": {
            "lung": {
                "Type": "Organ",
                "MaterialName": "G4_LUNG_ICRP",
                "Transform": {
                    "Rotation": {
                        "X": 0,
                        "Y": 0.0,
                        "Z": 0.0
                    },
                    "Translation": {
                        "X": 0.0,
                        "Y": 0.0,
                        "Z": 0.0
                    },
                    "Scale": {
                        "X": 1,
                        "Y": 1,
                        "Z": 1
                    }
                }
            }
        }
    }
}
"""

with open("temp_geom_description.json", "w") as f:
    f.write(desc_file)

def test_geometrysampling():
    sampler = GeometrySampler("temp_geom_description.json", {
        "patient": {
            "Translation": {
                "X": [-0.2, 0.2],
                "Y": [-0.0, 0.0],
                "Z": [-0.5, 0.5]
            }
        }
    })
    assert sampler is not None

    xs = []
    ys = []
    zs = []
    for i in range(2000):
        for obj_name, transformation in sampler._sample_transformation_per_object().items():
            x = transformation["Translation"]["X"]
            y = transformation["Translation"]["Y"]
            z = transformation["Translation"]["Z"]
            xs.append(x)
            ys.append(y)
            zs.append(z)

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    # check if the distribution was uniform
    assert np.isclose(np.mean(xs), 0, atol=0.1), f"Mean X {np.mean(xs)} not close to 0"
    assert np.isclose(np.mean(ys), 0, atol=0.1), f"Mean Y {np.mean(ys)} not close to 0"
    assert np.isclose(np.mean(zs), 0, atol=0.1), f"Mean Z {np.mean(zs)} not close to 0"
    assert np.isclose(np.std(xs), 0.4 / np.sqrt(12), rtol=0.1, atol=0.01), f"Std X {np.std(xs)} not close to {0.4 / np.sqrt(12)}"
    assert np.isclose(np.std(ys), 0.0, rtol=0.1, atol=0.1), f"Std Y {np.std(ys)} not close to 0.0"
    assert np.isclose(np.std(zs), 1.0 / np.sqrt(12), rtol=0.1, atol=0.02), f"Std Z {np.std(zs)} not close to {1.0 / np.sqrt(12)}"


GROUP_DESC = {
    "patient": {"Type": "patient", "Transform": {"Translation": {"X": 0.0, "Y": 0.5, "Z": 0.0}},
                "Children": {"lung": {"Transform": {"Translation": {"X": 0.0, "Y": 0.0, "Z": 0.0}}}}},
    "HardDesk": {"Transform": {"Translation": {"X": 0.1, "Y": -0.2, "Z": 0.0}}},
    "DetectorPlane": {"Transform": {"Translation": {"X": 0.0, "Y": 0.0, "Z": 0.0}}},
}


def _write_group_desc(tmp_path):
    import json
    path = tmp_path / "group.desc"
    path.write_text(json.dumps(GROUP_DESC))
    return str(path)


def _translation(desc, name):
    return GeometrySampler._find_mesh(desc, name)["Transform"]["Translation"]


def test_group_moves_its_meshes_rigidly(tmp_path):
    import json
    sampler = GeometrySampler(_write_group_desc(tmp_path), {
        "patient_on_table": {"Meshes": ["patient", "HardDesk"], "Translation": {"X": [-0.25, 0.25], "Z": [-0.5, 0.5]}}
    })
    out = tmp_path / "out.desc"
    for _ in range(50):
        sampler.sample_transformations(str(out))
        desc = json.loads(out.read_text())
        patient, desk = _translation(desc, "patient"), _translation(desc, "HardDesk")
        # one offset for both, added to each mesh's own translation: the arrangement is kept
        dx, dz = patient["X"] - 0.0, patient["Z"] - 0.0
        assert -0.25 <= dx <= 0.25 and -0.5 <= dz <= 0.5
        assert np.isclose(desk["X"], 0.1 + dx) and np.isclose(desk["Z"], 0.0 + dz)
        assert np.isclose(patient["Y"], 0.5) and np.isclose(desk["Y"], -0.2)
        # the child keeps its placement relative to the moving parent; untouched meshes stay put
        assert _translation(desc, "lung") == {"X": 0.0, "Y": 0.0, "Z": 0.0}
        assert _translation(desc, "DetectorPlane") == {"X": 0.0, "Y": 0.0, "Z": 0.0}


def test_mesh_entry_still_replaces(tmp_path):
    import json
    sampler = GeometrySampler(_write_group_desc(tmp_path), {"patient": {"Translation": {"Y": 0.1}}})
    out = tmp_path / "out.desc"
    sampler.sample_transformations(str(out))
    assert np.isclose(_translation(json.loads(out.read_text()), "patient")["Y"], 0.1)


def test_group_definition_errors(tmp_path):
    import pytest
    desc = _write_group_desc(tmp_path)
    with pytest.raises(ValueError, match="only define a 'Translation'"):
        GeometrySampler(desc, {"g": {"Meshes": ["patient"], "Rotation": {"X": [0.0, 1.0]}}})
    with pytest.raises(ValueError, match="more than one"):
        GeometrySampler(desc, {"g": {"Meshes": ["patient", "HardDesk"], "Translation": {"X": 0.1}}, "patient": {"Translation": {"X": 0.0}}})
    with pytest.raises(ValueError, match="non-empty list"):
        GeometrySampler(desc, {"g": {"Meshes": [], "Translation": {"X": 0.1}}})
    with pytest.raises(ValueError, match="not in"):
        GeometrySampler(desc, {"g": {"Meshes": ["patient", "Tabel"], "Translation": {"X": 0.1}}}).sample_transformations(str(tmp_path / "o.desc"))
