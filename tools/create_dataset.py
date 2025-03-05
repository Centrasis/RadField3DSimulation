import subprocess
import os
import re
import random
from typing import List, Union
from rich.progress import Progress, SpinnerColumn, TimeElapsedColumn
from rich import print
import argparse
import torch
import json
import shutil
from typing import NamedTuple
from logging import getLogger


pattern_energy = re.compile(r"Set X-Ray source energy to: ([\de\+\-]+)eV")
pattern_rotation = re.compile(r"Set X-Ray source roation \((\d+).\, (\d+).\)")
pattern_loc = re.compile(r"Set X-Ray source location to \(([\+\-\d\.]+)\, ([\+\-\d\.]+)\, ([\+\-\d\.]+)\)")
pattern_voxel_and_particles = re.compile(r"Start simulating with voxel dimension: ([\d\.\+]+)m and an particle count of ([\d\.]+)")


class ParameterValue(NamedTuple):
    name: str
    value: any

    @staticmethod
    def from_json(json: dict) -> "ParameterValue":
        name = json["name"]
        assert isinstance(name, str), "Name must be a string"
        name = name.lower()
        assert name in Parameter.VALID_NAMES, f"Parameter name must be one of {Parameter.VALID_NAMES} not {name}"
        return ParameterValue(
            name=name,
            value=json["value"]
        )


class Parameter(object):
    VALID_NAMES = ["energy", "source_distance", "source_angle_alpha", "source_angle_beta", "source_opening_angle", "source_shape", "source_spectra", "geometry", "tracer_algorithm", "bin_count", "voxel_size", "particles", "world_dim", "world_material"]

    def __init__(self, name: str, range: Union[tuple[int, int], tuple[float, float], list[any]], is_range: bool = None):
        self.name = name.lower()
        assert name in Parameter.VALID_NAMES, f"Parameter name must be one of {Parameter.VALID_NAMES} not {name}"
        self.range = range
        assert len(range) > 0, "Range must have at least one element"
        if is_range is None:
            self.is_range = len(range) == 2 and ((isinstance(range[0], int) and isinstance(range[1], int)) or (isinstance(range[0], float) and isinstance(range[1], float)))
        else:
            self.is_range = is_range

    def sample(self) -> any:
        if self.is_range:
            if isinstance(self.range[0], int):
                return random.randint(self.range[0], self.range[1])
            elif isinstance(self.range[0], float):
                return random.uniform(self.range[0], self.range[1])
            else:
                raise Exception("Range must be either int or float")
        else:
            idx = random.randint(0, len(self.range) - 1)
            return self.range[idx]

    @staticmethod
    def from_json(json: dict) -> "Parameter":
        range = json["range"] if "range" in json else None
        is_range = range is not None
        if range is None and "value" in json:
            range = [json["value"]]
        if range is None and "values" in json:
            range = json["values"]
        assert range is not None, "Range must be defined in the json object"

        return Parameter(
            name=json["name"],
            range=range,
            is_range=is_range
        )


class FixedParameterValueSet(object):
    def __init__(self, values: List[ParameterValue]):
        self.values = values


class ParameterSelector:
    def get_n_samples(self) -> int:
        pass

    def __next__(self) -> list[ParameterValue]:
        pass

    def __iter__(self):
        return self
    
    def __len__(self):
        return self.get_n_samples()
    
    def __getitem__(self, item) -> list[ParameterValue]:
        pass

    def get_max_energy(self) -> float:
        pass


class ParameterSequence(ParameterSelector):
    def __init__(self, parameters_file: str):
        file_content = json.load(open(parameters_file, "r"))
        assert "Metaparameters" in file_content, "Metaparameters must be defined in the sequence file"
        assert "ParameterSets" in file_content, "ParameterSets must be defined in the sequence file"
        self.parameter_sets: List[FixedParameterValueSet] = [
            FixedParameterValueSet(
                [
                    ParameterValue.from_json(val)
                    for val in pset
                ]
            )
            for pset in file_content["ParameterSets"]
        ]
        geometry_file = Parameter("geometry", [file_content["Metaparameters"]["GeometryFile"]]) if "GeometryFile" in file_content["Metaparameters"] else None
        particles = Parameter("particles", [file_content["Metaparameters"]["Particles"]]) if "Particles" in file_content["Metaparameters"] else None
        tracer_algo = Parameter("tracer_algorithm", [file_content["Metaparameters"]["TracerAlgorithm"]]) if "TracerAlgorithm" in file_content["Metaparameters"] else None
        bin_count = Parameter("bin_count", [file_content["Metaparameters"]["BinCount"]]) if "BinCount" in file_content["Metaparameters"] else None
        voxel_size = Parameter("voxel_size", [file_content["Metaparameters"]["VoxelSize"]]) if "VoxelSize" in file_content["Metaparameters"] else None
        world_dim = Parameter("world_dim", [file_content["Metaparameters"]["WorldDim"]]) if "WorldDim" in file_content["Metaparameters"] else None
        world_material = Parameter("world_material", [file_content["Metaparameters"]["WorldMaterial"]]) if "WorldMaterial" in file_content["Metaparameters"] else None

        for i, pset in enumerate(self.parameter_sets):
            all_names = [p.name for p in pset.values]
            if not any([p.name == "geometry" for p in pset.values]) and geometry_file is not None:
                pset.values.append(geometry_file)
            if not any([p.name == "particles" for p in pset.values]) and particles is not None:
                pset.values.append(particles)
            if not any([p.name == "tracer_algorithm" for p in pset.values]) and tracer_algo is not None:
                pset.values.append(tracer_algo)
            if not any([p.name == "bin_count" for p in pset.values]) and bin_count is not None:
                pset.values.append(bin_count)
            if not any([p.name == "voxel_size" for p in pset.values]) and voxel_size is not None:
                pset.values.append(voxel_size)
            if not any([p.name == "world_dim" for p in pset.values]) and world_dim is not None:
                pset.values.append(world_dim)
            if not any([p.name == "world_material" for p in pset.values]) and world_material is not None:
                pset.values.append(world_material)
            assert len(set(all_names)) == len(all_names), f"Parameter names must be unique, found duplicates in set {i}"

        self.idx = 0
        self.max_energy = file_content["Metaparameters"]["MaxEnergy"]

    def get_n_samples(self) -> int:
        return len(self.parameter_sets)
    
    def get_max_energy(self) -> float:
        return self.max_energy
    
    def __iter__(self):
        self.idx = 0
        return self
    
    def __next__(self) -> list[ParameterValue]:
        if self.idx >= len(self.parameter_sets):
            raise StopIteration
        self.idx += 1
        pset = self.parameter_sets[self.idx - 1]
        return pset.values
    
    def __getitem__(self, item) -> list[ParameterValue]:
        self.idx = item
        return next(self)


class ParameterizedSampler(ParameterSelector):
    def __init__(self, parameters: List[Parameter], n_samples: int, max_energy: float):
        self.parameters = parameters
        self.n_samples = n_samples
        self.max_energy = max_energy
        self.idx = 0

    @staticmethod
    def load_from(definition_file: str):
        file_content = json.load(open(definition_file, "r"))
        assert "Metaparameters" in file_content, "Metaparameters must be defined in the dataset definition file"
        max_energy = file_content["Metaparameters"]["MaxEnergy"]
        assert "Parameters" in file_content, "Parameters must be defined in the dataset definition file"
        parameters: List[Parameter] = [
            Parameter.from_json(p)
            for p in file_content["Parameters"]
        ]
        
        geometry_file = Parameter("geometry", [file_content["Metaparameters"]["GeometryFile"]]) if "GeometryFile" in file_content["Metaparameters"] else None
        particles = Parameter("particles", [file_content["Metaparameters"]["Particles"]]) if "Particles" in file_content["Metaparameters"] else None
        tracer_algo = Parameter("tracer_algorithm", [file_content["Metaparameters"]["TracerAlgorithm"]]) if "TracerAlgorithm" in file_content["Metaparameters"] else None
        bin_count = Parameter("bin_count", [file_content["Metaparameters"]["BinCount"]]) if "BinCount" in file_content["Metaparameters"] else None
        voxel_size = Parameter("voxel_size", [file_content["Metaparameters"]["VoxelSize"]]) if "VoxelSize" in file_content["Metaparameters"] else None
        world_dim = Parameter("world_dim", [file_content["Metaparameters"]["WorldDim"]]) if "WorldDim" in file_content["Metaparameters"] else None
        world_material = Parameter("world_material", [file_content["Metaparameters"]["WorldMaterial"]]) if "WorldMaterial" in file_content["Metaparameters"] else None
        if not any([p.name == "geometry" for p in parameters]) and geometry_file is not None:
            parameters.append(geometry_file)
        if not any([p.name == "particles" for p in parameters]) and particles is not None:
            parameters.append(particles)
        if not any([p.name == "tracer_algorithm" for p in parameters]) and tracer_algo is not None:
            parameters.append(tracer_algo)
        if not any([p.name == "bin_count" for p in parameters]) and bin_count is not None:
            parameters.append(bin_count)
        if not any([p.name == "voxel_size" for p in parameters]) and voxel_size is not None:
            parameters.append(voxel_size)
        if not any([p.name == "world_dim" for p in parameters]) and world_dim is not None:
            parameters.append(world_dim)
        if not any([p.name == "world_material" for p in parameters]) and world_material is not None:
            parameters.append(world_material)

        all_names = [p.name for p in parameters]
        assert len(set(all_names)) == len(all_names), "Parameter names must be unique"
        n_samples = file_content["Metaparameters"]["nSamples"]
        assert isinstance(n_samples, int), "nSamples must be an integer"
        return ParameterizedSampler(
            parameters=parameters,
            n_samples=n_samples,
            max_energy=max_energy
        )

    def get_n_samples(self) -> int:
        return self.n_samples
    
    def get_max_energy(self) -> float:
        return self.max_energy
    
    def __iter__(self):
        self.idx = 0
        return self
    
    def __next__(self) -> list[ParameterValue]:
        if self.idx >= self.n_samples:
            raise StopIteration
        self.idx += 1
        return [
            ParameterValue(
                name=p.name,
                value=p.sample()
            )
            for p in self.parameters
        ]
    
    def __getitem__(self, item) -> list[ParameterValue]:
        self.idx = item
        return next(self)


def write_spectum_file(src_file: str, out_path: str):
    if not os.path.exists(os.path.dirname(out_path)):
        try:
            os.makedirs(os.path.dirname(out_path))
        except FileExistsError:
            pass
        except Exception as e:
            getLogger().warning(f"Could not create directory {os.path.dirname(out_path)} for field {os.path.basename(src_file)} -> {e}")
            raise e

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
    parser.add_argument("--dataset_definition", default=None, type=str, nargs=1, required=False, help="Path to a dataset definition file that should be used. (Overrides energy/angle7/source_distance/source_shape/source_opening_angle sampling, energy_resolution, ...)")
    
    subcommands = parser.add_subparsers(title="Subcommands", description="valid subcommands", help="")
    cluster_parser = subcommands.add_parser("cluster", help="Cluster mode")
    cluster_parser.add_argument("--generate_batch", action="store_true", help="Generate cluster batch files")
    cluster_parser.add_argument("--type", type=str, required=False, default=None, help="Type of the cluster system ('slurm')")
    cluster_parser.add_argument("--node_partition", default=None, type=int, nargs=2, required=False, help="Partition the node into a cluster of nodes. (Number of nodes, Node ID)")
    cluster_parser.add_argument("--node_initialization_batch_path", default=None, type=str, required=False, help="Path to a batch file that should be executed on the cluster nodes before the calculation starts")
    cluster_parser.set_defaults(cluster=True, namespace="cluster")

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

    dataset_definition_file: str = None
    if "dataset_definition" in args and args.dataset_definition is not None:
        dataset_definition_file: str = args.dataset_definition[0] if isinstance(args.dataset_definition, list) else args.dataset_definition

    cluster_node_partition = None
    cluster_type = None
    cluster_should_generate_batch = False
    cluster_node_initialization_batch_path = None

    if "cluster" in args and args.cluster:
        if "type" in args and args.type is not None:
            cluster_type = args.type
            if cluster_type not in ["slurm"]:
                raise Exception("Cluster type must be one of ['slurm']")
        if "node_partition" in args and args.node_partition is not None:
            cluster_node_partition = args.node_partition
        if "generate_batch" in args and args.generate_batch is not None:
            cluster_should_generate_batch = args.generate_batch
        if "node_initialization_batch_path" in args and args.node_initialization_batch_path is not None:
            cluster_node_initialization_batch_path = args.node_initialization_batch_path
        
        if cluster_should_generate_batch and cluster_type is None:
            raise Exception("Cluster type must be defined to generate batch files")


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

    parameters: ParameterSelector = None
    
    if dataset_definition_file is not None and sequence_file is not None:
        raise Exception("Cannot use dataset definition and sequence file at the same time")
    
    if sequence_file is not None:
        parameters = ParameterSequence(sequence_file)
    
    if dataset_definition_file is not None:
        parameters = ParameterizedSampler.load_from(dataset_definition_file)
    
    if parameters is None:
        params = [
            Parameter("energy", energy_range),
            Parameter("source_distance", (source_distance, source_distance)),
            Parameter("source_angle_alpha", (-90, 90) if should_sample_angles else (source_angles[0], source_angles[0])),
            Parameter("source_angle_beta", (-45, 45) if should_sample_angles else (source_angles[1], source_angles[1])),
            Parameter("source_opening_angle", [opening_angle]),
            Parameter("source_shape", [source_shape]),
            Parameter("geometry", [geometry_file] if geometry_file != '' else []),
            Parameter("tracer_algorithm", [tracer_algorithm]),
            Parameter("bin_count", [simulation_energy_resolution]),
            Parameter("voxel_size", [voxel_size]),
            Parameter("particles", [particles]),
            Parameter("world_dim", [world_size]),
            Parameter("world_material", ["Air"])
        ]
        if spectra_path is not None:
            params.append(Parameter("source_spectra", [spectra_path]))
        parameters = ParameterizedSampler(params, n_sample_fields, Emax)
    n_sample_fields = len(parameters)


    if not os.path.isabs(out_dir):
        out_dir = os.path.join(os.getcwd(), out_dir)
    if geometry_file != '' and not os.path.isabs(geometry_file):
        geometry_file = os.path.join(os.getcwd(), geometry_file)
    if not os.path.isabs(binary_path):
        binary_path = os.path.join(os.getcwd(), binary_path)
    if not os.path.exists(out_dir):
        try:
            os.makedirs(out_dir)
        except FileExistsError:
            pass
        except Exception as e:
            getLogger().warning(f"Could not create directory {out_dir} -> {e}")
            raise e

    print(f"Start calculating {n_sample_fields} fields for dataset!")
    if cluster_node_partition is not None:
        print(f"Partitioning in cluster mode as node {cluster_node_partition[1]} of {cluster_node_partition[0]}")

    Emax = parameters.get_max_energy()

    preexisting_samples = [f for f in os.listdir(os.path.join(out_dir, "fields")) if f.endswith(".rf3")] if os.path.exists(os.path.join(out_dir, "fields")) else []
    try:
        preexisting_samples = [int(f.removesuffix(".rf3")) for f in preexisting_samples]
    except:
        getLogger().warning("Could not parse preexisting samples! Was not a standard dataset naming convention.")
    preexisting_samples_max_nb = (max(preexisting_samples) + 1) if len(preexisting_samples) > 0 else 0

    if cluster_node_partition is not None:
        parameters = [p for p in parameters]
        min_idx = cluster_node_partition[1] * (len(parameters) // cluster_node_partition[0])
        max_idx = (cluster_node_partition[1] + 1) * (len(parameters) // cluster_node_partition[0])
        max_idx = min(max_idx, len(parameters))
        if cluster_node_partition[1] + 1 >= cluster_node_partition[0]:
            max_idx = len(parameters)
        parameters = parameters[min_idx:max_idx]
        preexisting_samples_max_nb += min_idx

    cluster_batch_file_name = ""
    if cluster_should_generate_batch:
        cluster_batch_file_name = out_dir + ".sh"
        if os.path.exists(cluster_batch_file_name):
            os.remove(cluster_batch_file_name)

    if cluster_should_generate_batch:
        if cluster_type == "slurm":
            with open(cluster_batch_file_name, "w") as f:
                f.write("#!/bin/bash\n")
                f.write(f"#SBATCH --job-name={os.path.basename(out_dir)}\n")
                f.write(f"#SBATCH -o {os.path.basename(out_dir)}_%j.out\n")
                f.write(f"#SBATCH -e {os.path.basename(out_dir)}_%j.err\n")
                f.write("#SBATCH --ntasks-per-node=1\n")
                f.write(f"#SBATCH -N {cluster_node_partition[0]}\n")

                if cluster_node_initialization_batch_path is not None:
                    init_content = open(cluster_node_initialization_batch_path, "r").read()
                    f.write(init_content)
                    f.write("\n")
                f.write("\n")
        else:
            raise Exception(f"Cluster type {cluster_type} not supported")

    with Progress(SpinnerColumn(), *Progress.get_default_columns(), TimeElapsedColumn()) as progress:
        sampling_task = progress.add_task("Calculating field...", total=len(parameters))
        for nb_sample, sample_parameters in enumerate(parameters):
            spectra_path = None
            world_material = "Air"
            nb_sample += preexisting_samples_max_nb
            for param in sample_parameters:
                if param.name == "energy":
                    energy = param.value
                elif param.name == "source_distance":
                    source_distance = param.value
                elif param.name == "source_angle_alpha":
                    alpha = param.value
                elif param.name == "source_angle_beta":
                    beta = param.value
                elif param.name == "source_opening_angle":
                    source_opening_angle = param.value
                elif param.name == "source_shape":
                    source_shape = param.value
                elif param.name == "source_spectra":
                    spectra_path = param.value
                elif param.name == "geometry":
                    geometry_file = param.value
                    if not os.path.isabs(geometry_file):
                        geometry_file = os.path.join(os.getcwd(), geometry_file)
                elif param.name == "tracer_algorithm":
                    tracer_algorithm = param.value
                elif param.name == "bin_count":
                    simulation_energy_resolution = int(Emax / param.value)
                elif param.name == "voxel_size":
                    voxel_size = param.value
                elif param.name == "particles":
                    particles = param.value
                elif param.name == "world_dim":
                    world_size = param.value
                elif param.name == "world_material":
                    world_material = param.value
                else:
                    raise Exception(f"Parameter {param.name} was not recognized")   

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
                spec_csv_file = spec_csv_file.replace("\\", "/")
                spec_args = [
                    "--spectrum", spec_csv_file
                ]

            #out_name = f"RF_{'{:1.1f}'.format(energy/1000)}keV_{alpha}_{beta}_{os.path.basename(os.path.splitext(geometry_file)[0])}.rf3"
            out_name = f"{nb_sample:04d}.rf3"
            progress.update(sampling_task, description=f"Calculating field: {out_name}...")
            out_path = os.path.normpath(os.path.join(out_dir, "fields", out_name))

            if not os.path.exists(os.path.dirname(out_path)):
                try:
                    os.makedirs(os.path.dirname(out_path))
                except FileExistsError:
                    pass
                except Exception as e:
                    getLogger().warning(f"Could not create directory {os.path.dirname(out_path)} for field {out_name} -> {e}")
                    raise e

            if not args.clean:
                spec_args.append("--append")

            out_dir = out_dir.replace("\\", "/")
            out_path = out_path.replace("\\", "/")
            geometry_file = geometry_file.replace("\\", "/")
            binary_path = binary_path.replace("\\", "/")

            cmd_args = [
                        binary_path,
                        "--out", out_path,
                        "--max-energy", str(Emax),
                        "--source-alpha", str(alpha),
                        "--source-beta", str(beta),
                        "--source-distance", str(source_distance),
                        "--world-dim", f"{world_size[0]} {world_size[1]} {world_size[2]}",
                        "--particles", str(particles),
                        "--source-shape", source_shape,
                        "--voxel-dim", str(voxel_size),
                        "--energy-resolution", str(simulation_energy_resolution),
                        "--source-opening-angle", f"{source_opening_angle}",
                        "--tracing-algorithm", f"{tracer_algorithm}",
                        "--world-material", world_material
                    ] + spec_args + (["--geom", geometry_file] if geometry_file != '' else [])

            if cluster_should_generate_batch:
                err = None
                if cluster_type == "slurm":
                    with open(cluster_batch_file_name, "a") as f:
                        f.write("srun --exclusive --ntasks=1 --nodes=1 --no-requeue ")
                        f.write(" ".join(cmd_args))
                        f.write(" || true &\n")
                else:
                    raise Exception(f"Cluster type {cluster_type} not supported")
            else:
                out = subprocess.run(
                    cmd_args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    cwd=out_dir
                )
                err = out.stderr
                stdout = out.stdout.decode()
                ret_code = out.returncode

            progress.update(sampling_task, advance=1)

            if isinstance(err, str) or isinstance(err, bytes):
                if isinstance(err, bytes):
                    err = err.decode()
                if len(err) > 0:
                    print(stdout)
                    print(f"\n\nError while calculating {out_name}: '{err}'")
            if not cluster_should_generate_batch:
                print(f"[white]Field was written to -> [green]{out_path}" if os.path.exists(out_path) else f"[white]Field was [red]not [white]written to -> [red]{out_path}")
        if cluster_should_generate_batch:
            if cluster_type == "slurm":
                with open(cluster_batch_file_name, "a") as f:
                    f.write("wait")
