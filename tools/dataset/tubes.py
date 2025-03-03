from torch.utils.data import DataLoader
from .spectra import TubeDefinition, VOCABULARY
from typing import List, Tuple, Union
from torch import Tensor
from torch.utils.data import Dataset
import json
import torch
from torch.nn.utils.rnn import pad_sequence
from torch import nn
import os
from zipfile import ZipFile


class DatasetDescriptor(List[Tuple[str, str]]):
    def __init__(self, dataset_path: Union[str, None]) -> None:
        self._definitions = None
        self.dataset_path = ""
        self._is_zip = False
        if dataset_path is not None:
            self.dataset_path = dataset_path
            self._is_zip = os.path.isfile(self.dataset_path) and self.dataset_path.endswith(".zip")
            self._definitions = self.definitions

    @property
    def is_zip(self) -> bool:
        return self._is_zip
    
    @property
    def definitions(self) -> List[Tuple[str, str]]:
        if self._definitions is None:
            if self.is_zip:
                ds_zip = ZipFile(self.dataset_path)
                all_files = ds_zip.namelist()
                self._definitions = [(f.removesuffix(".spectrum") + ".info", f) for f in all_files if f.endswith('.spectrum')]
                ds_zip.close()
            else:
                self._definitions = [
                    (os.path.join(self.dataset_path, f), os.path.join(self.dataset_path, f.removesuffix(".info") + ".spectrum"))
                    for f in os.listdir(self.dataset_path) if f.endswith(".info")
                ]
        return self._definitions
    
    def __len__(self) -> int:
        return len(self.definitions)
    
    def __getitem__(self, idx: Union[int, slice]) -> Union[Tuple[str, str], List[Tuple[str, str]]]:
        return self.definitions[idx]

    def __iter__(self):
        return iter(self.definitions)
    
    def __contains__(self, item: Tuple[str, str]) -> bool:
        return item in self.definitions
    
    def __copy__(self) -> "DatasetDescriptor":
        new = DatasetDescriptor(None)
        new._definitions = self._definitions
        new._is_zip = self._is_zip
        new.dataset_path = self.dataset_path
        return new

    @staticmethod
    def from_list(desc: "DatasetDescriptor", lst: List[Tuple[str, str]]) -> "DatasetDescriptor":
        new_desc = desc.__copy__()
        new_desc._definitions = lst
        return new_desc


class TubesDataLoader(DataLoader):
    def __init__(self, dataset, batch_size=1, shuffle=False, num_workers=0, pin_memory=False, persistent_workers=False):
        super().__init__(
            dataset=dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            num_workers=num_workers,
            pin_memory=pin_memory,
            persistent_workers=persistent_workers,
            collate_fn=self.collate_fn
        )
        self.pad_value = VOCABULARY["[PAD]"]

    def collate_fn(self, batch):
        x, y = zip(*batch)
        padded_sequences = (
            pad_sequence(x, batch_first=True, padding_value=self.pad_value),
            pad_sequence(y, batch_first=True, padding_value=self.pad_value)
        )
        return padded_sequences


class TubesDataset(Dataset):
    def __init__(self, dataset: DatasetDescriptor, max_x: float = 1.0):
        self.dataset = dataset
        self.max_x = max_x
        self.ds_zip = None

    def __len__(self) -> int:
        return len(self.dataset)
    
    def normalize(self, spectrum: Tensor, norm_x: bool = True) -> Tensor:
        if norm_x:
            spectrum[:, 0] /= self.max_x
        spectrum[:, 1] /= spectrum[:, 1].sum(dim=-1)
        return spectrum

    def __getitem__(self, idx: int) -> Tuple[Tensor, Tensor]:
        info_file, spectrum_file = self.dataset[idx]
        if self.dataset.is_zip:
            if self.ds_zip is None:
                self.ds_zip = ZipFile(self.dataset.dataset_path)
            info_file = self.ds_zip.open(info_file)
            spectrum_file = self.ds_zip.open(spectrum_file, 'r')
        else:
            info_file = open(info_file)
            spectrum_file = open(spectrum_file, 'rb')
        try:
            definition = TubeDefinition.from_json(json.load(info_file))
            with torch.no_grad():
                spectrum: Tensor = torch.load(spectrum_file).permute(1, 0)
                spectrum = spectrum.to(torch.float32)
                spectrum = self.normalize(spectrum)

                assert spectrum[:, 1].min() != spectrum[:, 0].max(), f"{spectrum[:, 1].min()} != {spectrum[:, 0].max()}"
                assert spectrum[:, 0].min() != spectrum[:, 1].max(), f"{spectrum[:, 0].min()} != {spectrum[:, 1].max()}"

                return (definition.to_tensor(), spectrum)
        finally:
            info_file.close()
            spectrum_file.close()

    def __del__(self):
        if self.ds_zip is not None:
            self.ds_zip.close()


class AugmentatedTubeDataset(TubesDataset):
    def __init__(self, dataset: DatasetDescriptor, augmentations: List[nn.Module], max_x: float = 1.0, samling_rate: int = 2):
        super().__init__(dataset=dataset, max_x=max_x)
        self.samling_rate = samling_rate
        self.augmentations = augmentations

    def __len__(self) -> int:
        length = super().__len__()
        if len(self.augmentations) > 0:
            length *= self.samling_rate
        return length

    def __getitem__(self, idx: int) -> Tuple[Tensor, Tensor]:
        tube, spectrum = super().__getitem__(idx % len(self.dataset))
        for aug in self.augmentations:
            spectrum = aug(spectrum)
        if len(self.augmentations) > 0:
            spectrum = self.normalize(spectrum, False)
        return (tube, spectrum)
