from dataset.spectra import generate_sprectrum, TubeDefinition, LayerDefinition, MATERIALS
import random
import torch
import json
import os
from rich.progress import track
import logging
import argparse
from rich import print


def draw_chaotic_definition() -> TubeDefinition:
    filters = []
    for _ in range(random.randint(2, 8)):
        filters.append(
            LayerDefinition(
                material=MATERIALS[random.randint(0, len(MATERIALS) - 1)],
                thickness_mm=random.randint(int(1e+3), int(1e+6)) / 1e+3
            )
        )
    return TubeDefinition(
        energy=random.randint(70, 250) * 1000.0,
        anode_angle_deg=random.randint(10, 45),
        target_material='W',
        layers=filters
    )


def draw_realistic_definition() -> TubeDefinition:
    filters = []
    series = random.randint(0, 2)
    energy_kV = random.randint(10, 300)
    be_thickness = random.randint(1, 7) * 1.0
    al_thickness = random.randint(1, 3) * 1.0
    cu_thickness = random.randint(0, 3) * 1.0
    sn_thickness = random.randint(0, 3) * 1.0
    pb_thickness = random.randint(0, 3) * 1.0

    if series == 0: # PTB-Series A like setup
        filters.append(LayerDefinition(material='Be', thickness_mm=be_thickness))
        filters.append(LayerDefinition(material='Al', thickness_mm=al_thickness))
        if energy_kV >= 40 and cu_thickness > 0.0:
            filters.append(LayerDefinition(material='Cu', thickness_mm=cu_thickness))
        if energy_kV >= 120 and sn_thickness > 0.0:
            filters.append(LayerDefinition(material='Sn', thickness_mm=sn_thickness))
        if energy_kV >= 200 and pb_thickness > 0.0:
            filters.append(LayerDefinition(material='Pb', thickness_mm=pb_thickness))
    elif series == 1: # PTB-Series B like setup
        filters.append(LayerDefinition(material='Be', thickness_mm=be_thickness))
        filters.append(LayerDefinition(material='Al', thickness_mm=al_thickness))
        if energy_kV >= 60 and energy_kV < 150 and cu_thickness > 0.0:
            filters.append(LayerDefinition(material='Cu', thickness_mm=cu_thickness))
        elif energy_kV >= 150 and sn_thickness > 0.0:
            filters.append(LayerDefinition(material='Sn', thickness_mm=sn_thickness))
    elif series == 2: # PTB-Series C like setup
        filters.append(LayerDefinition(material='Be', thickness_mm=be_thickness))
        if energy_kV >= 20:
            filters.append(LayerDefinition(material='Al', thickness_mm=al_thickness))
        if energy_kV >= 100 and cu_thickness > 0.0:
            filters.append(LayerDefinition(material='Cu', thickness_mm=cu_thickness))

    filters.append(LayerDefinition(material='Air', thickness_mm=random.randint(1, 25) * 100.0))
    return TubeDefinition(
        energy=energy_kV * 1000.0,
        anode_angle_deg=random.randint(10, 15),
        target_material='W',
        layers=filters
    )


def draw_120kV_definition(energy_kV: float = 100.0) -> TubeDefinition:
    filters = []
    filters.append(LayerDefinition(material='Be', thickness_mm=1))
    filters.append(LayerDefinition(material='Cu', thickness_mm=4.028))
    filters.append(LayerDefinition(material='Al', thickness_mm=0.151))
    return TubeDefinition(
        energy=energy_kV * 1000.0,
        anode_angle_deg=35,
        target_material='W',
        layers=filters
    )

def draw_N120_definition(energy_kV: float = 100.0) -> TubeDefinition:
    filters = []
    filters.append(LayerDefinition(material='Be', thickness_mm=1))
    filters.append(LayerDefinition(material='Al', thickness_mm=4.0))
    filters.append(LayerDefinition(material='Cu', thickness_mm=5.0))
    filters.append(LayerDefinition(material='Sn', thickness_mm=1.0))
    return TubeDefinition(
        energy=energy_kV * 1000.0,
        anode_angle_deg=35,
        target_material='W',
        layers=filters
    )

def draw_miniX2_definition(energy_kV: float = 50.0) -> TubeDefinition:
    filters = []
    filters.append(LayerDefinition(material='Be', thickness_mm=0.125))
    return TubeDefinition(
        energy=energy_kV * 1000.0,
        anode_angle_deg=0,
        target_material='Au',
        layers=filters
    )


def draw_carm_definition() -> TubeDefinition:
    energy_kV = random.randint(40, 125)
    filters = []
    filters.append(LayerDefinition(material='Al', thickness_mm=random.randint(2, 7) + 0.5)) # Min 2.5mm 7.5mm in 1mm steps
    cu_mm = random.randint(0, 9) * 0.1
    if cu_mm != 0:
        filters.append(LayerDefinition(material='Cu', thickness_mm=cu_mm))
    return TubeDefinition(
        energy=energy_kV * 1000.0,
        anode_angle_deg=random.randint(8, 12),
        target_material='W',
        layers=filters
    )


def draw_artisQ_carm_definition() -> TubeDefinition:
    energy_kV = random.randint(40, 125)
    filters = []
    filters.append(LayerDefinition(material='Al', thickness_mm=random.randint(2, 7) + 0.5)) # Min 2.5mm 7.5mm in 1mm steps
    cu_mm = random.randint(0, 9) * 0.1
    if cu_mm != 0:
        filters.append(LayerDefinition(material='Cu', thickness_mm=cu_mm))
    return TubeDefinition(
        energy=energy_kV * 1000.0,
        anode_angle_deg=9.5,
        target_material='W',
        layers=filters
    )


MODES: dict = {
    "Chaos": draw_chaotic_definition,
    "Realistic": draw_realistic_definition,
    "C-Arm": draw_carm_definition,
    "ArtisQ": draw_artisQ_carm_definition,
    "120kV_C100": lambda: draw_120kV_definition(100.0),
    "N120": lambda: draw_N120_definition(120),
    "MiniX2_50keV": lambda: draw_miniX2_definition(50.0)
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate dataset of X-ray spectra")
    parser.add_argument("--size", type=int, default=100000, help="Size of the dataset")
    parser.add_argument("--output", type=str, default="datasets/Spectra-1", help="Output directory")
    parser.add_argument("--mode", type=str, default="Chaos", help=f"Mode of the dataset ({MODES.keys()})")
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING)
    dataset_size = args.size
    dataset_path = args.output
    generator_fn = MODES[args.mode]

    print(f"Generating dataset of size {dataset_size} in {dataset_path} using {args.mode} mode")

    if not os.path.exists(dataset_path):
        os.makedirs(dataset_path)

    for _ in track(range(int(dataset_size)), description="Generating Dataset"):
        definition = generator_fn()
        spectrum = generate_sprectrum(definition)
        if spectrum[0,:].min() == spectrum[0,:].max() or spectrum[1,:].min() == spectrum[1,:].max():
            logging.warning("Spectrum is empty, skipping...")
            continue
        id_str = str(random.getrandbits(64))
        torch.save(spectrum, f"{dataset_path}/{id_str}.spectrum")
        with open(f"{dataset_path}/{id_str}.info", "w") as f:
            f.write(json.dumps(definition.__json__()))
