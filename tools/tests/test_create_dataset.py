import sys
import os

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

from create_dataset import GeometrySampler
import numpy as np


desc_file = """
{
    "patient": {
        "MaterialName": "G4_TISSUE_SOFT_ICRU-4",
        "Patient": true,
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
                "Patient": false,
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
