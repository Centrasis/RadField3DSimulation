import subprocess
import os
import re
import random
from typing import List
from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn
from rich import print
import argparse
import torch
import json
import shutil


pattern_energy = re.compile(r"Set X-Ray source energy to: ([\de\+\-]+)eV")
pattern_rotation = re.compile(r"Set X-Ray source roation \((\d+).\, (\d+).\)")
pattern_loc = re.compile(r"Set X-Ray source location to \(([\+\-\d\.]+)\, ([\+\-\d\.]+)\, ([\+\-\d\.]+)\)")
pattern_voxel_and_particles = re.compile(r"Start simulating with voxel dimension: ([\d\.\+]+)m and an particle count of ([\d\.]+)")


class Parameter:
    def __init__(self, energy: int, source_distance: float, source_angle_alpha: int, source_angle_beta: int, source_opening_angle: float):
        self.energy = energy
        self.source_distance = source_distance
        self.source_angle_alpha = source_angle_alpha
        self.source_angle_beta = source_angle_beta
        self.source_opening_angle = source_opening_angle

    @staticmethod
    def from_json(json: dict) -> "Parameter":
        return Parameter(
            energy=json["energy"],
            source_distance=json["source_distance"],
            source_angle_alpha=json["source_angle_alpha"],
            source_angle_beta=json["source_angle_beta"],
            source_opening_angle=json["source_opening_angle"] if "source_opening_angle" in json else None
        )


class ParameterSelector:
    def get_n_samples(self) -> int:
        pass

    def __next__(self) -> Parameter:
        pass

    def __iter__(self):
        return self
    
    def __len__(self):
        return self.get_n_samples()
    
    def __getitem__(self, item) -> Parameter:
        pass

    def get_max_energy(self) -> float:
        pass


class ParameterSampler(ParameterSelector):
    def __init__(self, energy_range: tuple[float, float], source_distance_range: tuple[float, float], source_angle_alpha_range: tuple[int, int], source_angle_beta_range: tuple[int, int], energy_resolution: float, source_opening_angle: tuple[float, float], n_samples: int):
        self.energy_range = energy_range
        self.source_distance_range = source_distance_range
        self.source_angle_alpha_range = source_angle_alpha_range
        self.source_angle_beta_range = source_angle_beta_range
        self.n_samples = n_samples
        self.energy_resolution = energy_resolution
        self.source_opening_angle = source_opening_angle
        self.idx = 0

    def get_n_samples(self) -> int:
        return self.n_samples
    
    def get_max_energy(self) -> float:
        return self.energy_range[1]
    
    def __iter__(self):
        self.idx = 0
        return self
    
    def __next__(self) -> Parameter:
        if self.idx >= self.n_samples:
            raise StopIteration
        self.idx += 1
        return Parameter(
            energy=self.energy_range[0] + ((random.randint(0, int(self.energy_resolution)) / self.energy_resolution) * (self.energy_range[1] - self.energy_range[0])),
            source_distance=random.uniform(self.source_distance_range[0], self.source_distance_range[1]),
            source_angle_alpha=random.randint(int(self.source_angle_alpha_range[0]), int(self.source_angle_alpha_range[1])),
            source_angle_beta=random.randint(int(self.source_angle_beta_range[0]), int(self.source_angle_beta_range[1])),
            source_opening_angle=random.uniform(self.source_opening_angle[0], self.source_opening_angle[1])
        )

    def __getitem__(self, item) -> Parameter:
        self.idx = item
        return next(self)
    

class ParameterSequence(ParameterSelector):
    def __init__(self, parameters_file: str, default_source_opening_angle: float):
        self.parameters: List[Parameter] = [
            Parameter.from_json(elem)
            for elem in json.load(open(parameters_file, "r"))
        ]
        self.idx = 0
        self.default_source_opening_angle = default_source_opening_angle
        self.max_energy = max([p.energy for p in self.parameters])
    
    def get_n_samples(self) -> int:
        return len(self.parameters)
    
    def get_max_energy(self) -> float:
        return self.max_energy
    
    def __iter__(self):
        self.idx = 0
        return self
    
    def __next__(self) -> Parameter:
        if self.idx >= len(self.parameters):
            raise StopIteration
        self.idx += 1
        element = self.parameters[self.idx - 1]
        if "source_opening_angle" not in element.__dict__ or element.source_opening_angle is None:
            element["source_opening_angle"] = self.default_source_opening_angle
        return element
    
    def __getitem__(self, item) -> Parameter:
        element: Parameter = self.parameters[item]
        if "source_opening_angle" not in element.__dict__ or element.source_opening_angle is None:
            element.source_opening_angle = self.default_source_opening_angle
        return element


def write_spectum_file(src_file: str, out_path: str):
    if not os.path.exists(os.path.dirname(out_path)):
        os.makedirs(os.path.dirname(out_path))

    spectrum = torch.load(src_file, weights_only=True)

    energies = spectrum[0]
    fluence = spectrum[1]

    # write histogram to file
    with open(out_path, "w") as f:
        f.write(f"Energy[eV]    Fluence[]")
        for i in range(len(energies)):
            f.write(f"\n{energies[i] * 1000}    {fluence[i]}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Calculate Radiation Field")
    parser.add_argument("--dest", default="Dataset", type=str, nargs=1, required=False, help="Destination for the field outputs")
    parser.add_argument("--fields", default=10, type=int, nargs=1, required=False, help="Number of fields to calculate")
    parser.add_argument("--Emin", default=1e+2, type=float, nargs=1, required=False, help="Emin in eV for the gamma particles radiation")
    parser.add_argument("--Emax", default=1e+5, type=float, nargs=1, required=False, help="Emax in eV for the gamma particles radiation")
    parser.add_argument("--particles", default=1.0, type=float, nargs=1, required=False, help="Maximum number of particles to simulate")
    parser.add_argument("--geometry", type=str, nargs=1, default='', required=False, help="Path to a geometry file (practically any format) where a material definition file is stored as well (*.vm)")
    parser.add_argument("--voxel_size", default=0.05, type=float, nargs=1, required=False, help="Dimension of the cubic voxels in m")
    parser.add_argument("--world_size", default=[1, 1, 1], type=float, nargs=3, required=False, help="Dimension of the rectangular world in m")
    parser.add_argument("--energy_res", default=1e+2, type=float, nargs=1, required=False, help="Energy resolution to use in eV for the sampling of new energies during dataset creation.")
    parser.add_argument("--binary", default="RadField3D.exe", type=str, nargs=1, required=False, help="Path to RadField3D Binary")
    parser.add_argument("--spectra", default=None, type=str, nargs=1, required=False, help="Path to a folder of a spectra dataset or a single spectrum that should be used. (Disables energy sampling)")
    parser.add_argument("--source_distance", default=1.0, type=float, nargs=1, required=False, help="Distance of the source to the geometry in m")
    parser.add_argument("--source_shape", default="cone", type=str, nargs=1, required=False, help="Shape of the source (cone, rectangle)")
    parser.add_argument("--source_angles", default=None, type=float, nargs=2, required=False, help="Source angles in degrees (alpha, beta). If not set the angles will be randomly sampled.")
    parser.add_argument("--source_opening_angle", default='20.0', type=str, nargs='+', required=False, help="Opening angle of the source in degrees (One to two values depending on the source shape)")
    parser.add_argument("--clean", default=False, action="store_true", required=False, help="Clean the output radiation field before calculating it new. Otherwise append, if the fields are representing the same.")
    parser.add_argument("--bin_count", default=None, required=False, type=float, help="Optional: Define the number of energy bins to store for each fluence in each voxel. Defaults to match a bin width of 1 eV.")
    parser.add_argument("--sequence_file", default=None, type=str, nargs=1, required=False, help="Path to a sequence file that should be used. (Disables energy and angle sampling)")
    parser.add_argument("--tracer_algorithm", default="sampling", type=str, nargs=1, required=False, help="Tracer algorithm to use (sampling, bresenham, linetracing)")
    parser.add_argument("--cluster_node_partition", default=None, type=int, nargs=2, required=False, help="Partition the node into a cluster of nodes. (Number of nodes, Node ID)")

    args = parser.parse_args()

    out_dir: str = args.dest[0] if isinstance(args.dest, list) else args.dest
    n_sample_fields: int = args.fields[0] if isinstance(args.fields, list) else args.fields
    Emin: float = args.Emin[0] if isinstance(args.Emin, list) else args.Emin
    Emax = None
    if "Emax" in args and args.Emax is not None:
        Emax: float = args.Emax[0] if isinstance(args.Emax, list) else args.Emax
    energy_range = [Emin, Emax]
    particles: float = args.particles[0] if isinstance(args.particles, list) else args.particles
    energy_resolution: float = args.energy_res[0] if isinstance(args.energy_res, list) else args.energy_res
    world_size: List[float] = args.world_size
    voxel_size: float = args.voxel_size[0] if isinstance(args.voxel_size, list) else args.voxel_size
    geometry_file: str = args.geometry[0] if isinstance(args.geometry, list) else args.geometry
    binary_path: str = args.binary[0] if isinstance(args.binary, list) else args.binary
    spectra_path: str = args.spectra[0] if isinstance(args.spectra, list) else args.spectra
    source_distance: float = args.source_distance[0] if isinstance(args.source_distance, list) else args.source_distance
    source_shape: str = args.source_shape[0] if isinstance(args.source_shape, list) else args.source_shape
    source_angles: List[float] = args.source_angles
    should_sample_angles = source_angles is None
    source_opening_angle = None
    tracer_algorithm: str = args.tracer_algorithm[0] if isinstance(args.tracer_algorithm, list) else args.tracer_algorithm
    if "source_opening_angle" in args and args.source_opening_angle is not None:
        source_opening_angle: str = ' '.join(args.source_opening_angle) if isinstance(args.source_opening_angle, list) else args.source_opening_angle
    simulation_energy_resolution: float = int(Emax / args.bin_count) if args.bin_count is not None else 1e+3
    sequence_file: str = None
    if "sequence_file" in args and args.sequence_file is not None:
        sequence_file: str = args.sequence_file[0] if isinstance(args.sequence_file, list) else args.sequence_file

    cluster_node_partition = None
    if "cluster_node_partition" in args and args.cluster_node_partition is not None:
        cluster_node_partition = args.cluster_node_partition

    opening_angle = source_opening_angle
    if opening_angle is None:
        opening_angle = (1, 30)
    elif isinstance(opening_angle, str):
        opening_angle = opening_angle.split(" ")
        if len(opening_angle) == 1:
            opening_angle = (float(opening_angle[0]), float(opening_angle[0]))
        elif len(opening_angle) == 2:
            opening_angle = (float(opening_angle[0]), float(opening_angle[1]))
        else:
            raise Exception("Opening angle must have one or two values")
    elif isinstance(opening_angle, list) or isinstance(opening_angle, tuple):
        if len(opening_angle) == 1:
            opening_angle = (opening_angle[0], opening_angle[0])
        elif len(opening_angle) == 2:
            opening_angle = (opening_angle[0], opening_angle[1])
        else:
            raise Exception("Opening angle must have one or two values")
    elif isinstance(opening_angle, float):
        opening_angle = (opening_angle, opening_angle)
    else:
        raise Exception("Opening angle must be a string, list, tuple or float")

    parameters: ParameterSelector = ParameterSampler(
        energy_range=energy_range,
        source_distance_range=(source_distance, source_distance),
        source_angle_alpha_range=(-90, 90) if should_sample_angles else (source_angles[0], source_angles[0]),
        source_angle_beta_range=(-45, 45) if should_sample_angles else (source_angles[1], source_angles[1]),
        n_samples=n_sample_fields,
        energy_resolution=energy_resolution,
        source_opening_angle=opening_angle
    ) if sequence_file is None else ParameterSequence(sequence_file, default_source_opening_angle=source_opening_angle)

    n_sample_fields = len(parameters)


    if not os.path.isabs(out_dir):
        out_dir = os.path.join(os.getcwd(), out_dir)
    if geometry_file != '' and not os.path.isabs(geometry_file):
        geometry_file = os.path.join(os.getcwd(), geometry_file)
    if not os.path.isabs(binary_path):
        binary_path = os.path.join(os.getcwd(), binary_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    print(f"Start calculating {n_sample_fields} fields for dataset!")
    if cluster_node_partition is not None:
        print(f"Partitioning in cluster mode as node {cluster_node_partition[1]} of {cluster_node_partition[0]}")

    if Emax is None:
        Emax = parameters.get_max_energy()

    parameters = [p for p in parameters]

    if cluster_node_partition is not None:
        min_idx = cluster_node_partition[1] * (len(parameters) // cluster_node_partition[0])
        max_idx = (cluster_node_partition[1] + 1) * (len(parameters) // cluster_node_partition[0])
        max_idx = min(max_idx, len(parameters))
        if cluster_node_partition[1] + 1 >= cluster_node_partition[0]:
            max_idx = len(parameters)
        parameters = parameters[min_idx:max_idx]

    with Progress(SpinnerColumn(), *Progress.get_default_columns(), TimeElapsedColumn()) as progress:
        sampling_task = progress.add_task("Calculating field...", total=len(parameters))
        for sample in parameters:
            energy = sample.energy
            alpha = sample.source_angle_alpha
            beta = sample.source_angle_beta

            spec_args = []
            if spectra_path is not None:
                spec_file = ""
                spec_csv_file = None
                if os.path.isdir(spectra_path):
                    spec_files = [os.path.join(spectra_path, f) for f in os.listdir(spectra_path) if f.endswith(".spectrum")]
                    spec_file = spec_files[random.randint(0, len(spec_files) - 1)]
                elif spectra_path.endswith(".spectrum"):
                    spec_file = spectra_path
                elif spectra_path.endswith(".csv"):
                    spec_csv_file = spectra_path
                    energy = 0.0
                else:
                    raise Exception(f"Spectra file was not valid: {spec_file}")

                if spec_csv_file is None:
                    info_file = spec_file.removesuffix(".spectrum") + ".info"
                    if not os.path.exists(info_file):
                        raise Exception(f"Spec Info file not found: {info_file}")
                    info = json.loads(open(info_file, "r").read())

                    energy = float(info["energy"])

                    spec_csv_file = os.path.join(out_dir, "spectra", os.path.basename(spec_file).removesuffix(".spectrum") + ".csv")
                    write_spectum_file(
                        src_file=spec_file,
                        out_path=spec_csv_file
                    )
                    if not os.path.exists(os.path.join(out_dir, os.path.basename(info_file))):
                        shutil.copyfile(info_file, os.path.join(out_dir, "spectra", os.path.basename(info_file)))
                spec_args = [
                    "--spectrum", os.path.normpath(spec_csv_file)
                ]

            out_name = f"RF_{'{:1.1f}'.format(energy/1000)}keV_{alpha}_{beta}_{os.path.basename(os.path.splitext(geometry_file)[0])}.rf3"
            progress.update(sampling_task, description=f"Calculating field: {out_name}...")
            out_path = os.path.normpath(os.path.join(out_dir, "fields", out_name))

            if not os.path.exists(os.path.dirname(out_path)):
                os.makedirs(os.path.dirname(out_path))

            if not args.clean:
                spec_args.append("--append")

            out = subprocess.run([
                    binary_path,
                    "--out", out_path,
                    "--max-energy", str(Emax),
                    "--source-alpha", str(alpha),
                    "--source-beta", str(beta),
                    "--source-distance", str(sample.source_distance),
                    "--world-dim", f"{world_size[0]} {world_size[1]} {world_size[2]}",
                    "--particles", str(particles),
                    "--source-shape", source_shape,
                    "--voxel-dim", str(voxel_size),
                    "--energy-resolution", str(simulation_energy_resolution),
                    "--source-opening-angle", f"{sample.source_opening_angle}",
                    "--tracing-algorithm", f"{tracer_algorithm}"
                ] + spec_args + (["--geom", geometry_file] if geometry_file != '' else []),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=out_dir
            )
            err = out.stderr
            stdout = out.stdout.decode()

            progress.update(sampling_task, advance=1)

            if isinstance(err, str) or isinstance(err, bytes):
                if isinstance(err, bytes):
                    err = err.decode()
                if len(err) > 0:
                    print(stdout)
                    print(f"\n\nError while calculating {out_name}: '{err}'")
            print(f"[white]Field was written to -> [green]{out_path}" if os.path.exists(out_path) else f"[white]Field was [red]not [white]written to -> [red]{out_path}")
