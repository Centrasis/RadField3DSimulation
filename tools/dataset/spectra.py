import spekpy as sp
from typing import List, NamedTuple, Dict
from torch import Tensor
import torch
from torch import nn



#MATERIALS = [
#    "W",
#    "Be",
#    "Sn",
#    "Cu",
#    "Pb",
#    "Al",
#    "Air"
#]

## PTB Materials
MATERIALS = [
    "W",
    "Au",
    "Be",
    "Al",
    "Cu",
    "Sn",
    "Pb",
    "Air"
]

VOCABULARY = {
    '[PAD]': -100,
    '[UNK]': 1,
    '[CLS]': 2,
    '[SEP]': 3
}
VOCABULARY.update({mat: i + len(VOCABULARY.keys()) for i, mat in enumerate(sorted(set(MATERIALS)))})
VOCABULARY.update({str(i): i + len(VOCABULARY.keys()) for i in range(25001)})
VOCABULARY.update({
    "[Energy]": len(VOCABULARY),
    "[AnodeAngle]": len(VOCABULARY) + 1
})


class LayerDefinition(NamedTuple):
    material: str
    thickness_mm: float


class TubeDefinition(object):
    layers: List[LayerDefinition]
    target_material: str
    anode_angle_deg: float
    energy: float
    VOCABULARY: Dict[str, int] = None

    def __init__(self, energy: float, target_material: str, anode_angle_deg: float, layers: List[LayerDefinition], max_energy: float = 2.5e+5, max_layer_thickness: float = 1e+3, use_spekpy: bool = True) -> None:
        self.energy = energy
        self.anode_angle_deg = anode_angle_deg
        self.target_material = target_material
        self.layers = layers
        self.max_energy = max_energy
        self.max_layer_thickness = max_layer_thickness
        self._spekpy: sp.Spek = self.create_SPEK() if use_spekpy else None
        self.embedding_dim = 768
        self.embedding = nn.Embedding(len(VOCABULARY.keys()), self.embedding_dim)

    def to_tensor(self) -> Tensor:
        encode_layers = []

        #Ignore the constant Tungsten
        #encode_layers.append(self.target_material)

        e_steps = int(self.energy / 1000)
        encode_layers.append("[CLS]")
        encode_layers.append("[Energy]")
        encode_layers.append(str(e_steps))
        
        encode_layers.append("[AnodeAngle]")

        self.anode_angle_deg = int(self.anode_angle_deg)
        if self.anode_angle_deg < 0:
            self.anode_angle_deg = 360 - self.anode_angle_deg
        if self.anode_angle_deg >= 360:
            self.anode_angle_deg = self.anode_angle_deg % 360
        encode_layers.append(str(self.anode_angle_deg))

        encode_layers.append("[SEP]")

        for layer in self.layers:
            encode_layers.append(layer.material)
            count = int(layer.thickness_mm / 0.1)
            encode_layers.append(str(count))
        layer_ids = torch.tensor([(VOCABULARY[layer] if layer in VOCABULARY else VOCABULARY["[UNK]"]) for layer in encode_layers], dtype=torch.long)
        return layer_ids

    def create_SPEK(self) -> sp.Spek:
        spek_def = sp.Spek(kvp=self.energy/1000, th=self.anode_angle_deg, targ=self.target_material)
        for layer in self.layers:
            spek_def.filter(layer.material, layer.thickness_mm)
        return spek_def
    
    def __json__(self):
        return {
            "energy": self.energy,
            "anode_angle_deg": self.anode_angle_deg,
            "target_material": self.target_material,
            "layers": [
                {
                    "material": layer.material,
                    "thickness_mm": layer.thickness_mm
                } for layer in self.layers
            ]
        }
    
    @staticmethod
    def from_json(json: dict):
        return TubeDefinition(
            energy=json["energy"],
            anode_angle_deg=json["anode_angle_deg"],
            target_material=json["target_material"],
            layers=[
                LayerDefinition(
                    material=layer["material"],
                    thickness_mm=layer["thickness_mm"]
                ) for layer in json["layers"]
            ],
            use_spekpy=False
        )


def generate_sprectrum(definition: TubeDefinition) -> Tensor:
    spectrum = definition._spekpy.get_spectrum()
    return torch.stack([
        torch.from_numpy(spectrum[0]),
        torch.from_numpy(spectrum[1])
    ])


def write_spectrum_file(definition: TubeDefinition, out_path: str):
    # get spectrum from setup
    spectrum = definition._spekpy.get_spectrum()

    # build spectrum histrogram
    energies = spectrum[0]
    fluence = spectrum[1]

    # write histogram to file
    with open(out_path, "w") as f:
        f.write(f"Energy[eV]    Fluence[]")
        for i in range(len(energies)):
            f.write(f"\n{energies[i] * 1000}    {fluence[i]}")
