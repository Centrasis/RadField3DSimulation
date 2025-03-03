from NeuralDashboard.dataset import DatasetLoader
from .tubes import TubesDataLoader, AugmentatedTubeDataset, TubesDataset, DatasetDescriptor
from torch.utils.data import Dataset
from torch.utils.data import random_split
from lightning.pytorch import LightningDataModule
from NeuralDashboard.datacontainers import Points2D
from torch import nn
from typing import List
import torch
import os
import multiprocessing


class TubesDataModule(LightningDataModule, DatasetLoader):
    N_WORKERS = max(0, int((multiprocessing.cpu_count() - 1) / 4))

    def __init__(self, data_dir: str, batch_size: int = 32, augmentations: List[nn.Module] = [], samling_rate=2) -> None:
        LightningDataModule.__init__(self)
        self.full_dataset = None
        self.samling_rate = samling_rate
        self.batch_size = batch_size
        self.data_dir = data_dir
        self.augmentations = augmentations
        self.init_info(
            name=os.path.basename(data_dir),
            description="Radiation Field Dataset",
            path=data_dir,
            data_item_info=Points2D(representation="bars")
        )
        self.max_energy = 300.0
        self.dataset_desc = DatasetDescriptor(data_dir)

    def guess_max_bin_count(self, sampling: float = 0.01) -> int:
        max_size = 0
        self.max_energy = 0.0
        files = DatasetDescriptor.from_list(self.dataset_desc, self.dataset_desc[0:int(len(self.dataset_desc) * sampling)])
        for _, spectrum in TubesDataset(files):
            if spectrum.shape[0] > 0:
                max_size = max(max_size, spectrum.shape[0])
                self.max_energy = max(self.max_energy, spectrum[:, 0].max())
        self.max_energy = self.max_energy * 1.1
        return int(max_size * 1.1)

    def setup(self, stage: str):
        if stage is not None:
            self.dataset_info.session_info.item_count = len(self.dataset_desc)
            self.dataset_info.session_info.loaded = False
            self.dataset_info.session_info.loading_inprogress = True
            val_size = int(len(self.dataset_desc) * 0.1)
            test_size = int(len(self.dataset_desc) * 0.1)
            self.train_set, self.val_set, self.test_set = random_split(
                self.dataset_desc.definitions, [len(self.dataset_desc) - val_size - test_size, val_size, test_size], generator=torch.Generator().manual_seed(42)
            )
            self.full_dataset = AugmentatedTubeDataset(self.dataset_desc, self.augmentations, self.max_energy, self.samling_rate)
            self.set_active_dataset(self)
            self.train_set = DatasetDescriptor.from_list(self.dataset_desc, self.train_set)
            self.val_set = DatasetDescriptor.from_list(self.dataset_desc, self.val_set)
            self.test_set = DatasetDescriptor.from_list(self.dataset_desc, self.test_set)
            self.dataset_info.session_info.loaded = True
            self.dataset_info.session_info.loading_inprogress = False

    def __len__(self) -> int:
        if self.full_dataset is not None:
            length = int(len(self.full_dataset) * self.samling_rate)
            self.dataset_info.session_info.item_count = length
            return length
        return -1
    
    def __getitem__(self, idx: int):
        return self.full_dataset[idx][1]
    
    def __data_set_info__(self):
        if self.full_dataset is not None:
            self.dataset_info.session_info.item_count = len(self.full_dataset)
        return self.dataset_info
    
    def __data_set_statistics__(self) -> dict:
        if self.full_dataset is None:
            return {}

        return {
            "n_train": len(self.train_set),
            "n_val": len(self.val_set),
            "n_test": len(self.test_set)
        }

    def train_dataloader(self):
        ds = AugmentatedTubeDataset(self.train_set, self.augmentations, self.max_energy, self.samling_rate)
        return TubesDataLoader(
            ds,
            batch_size=self.batch_size,
            shuffle=True,
            num_workers=TubesDataModule.N_WORKERS,
            persistent_workers=TubesDataModule.N_WORKERS > 0
        )
    
    def val_dataloader(self):
        ds = AugmentatedTubeDataset(self.val_set, self.augmentations, self.max_energy, self.samling_rate)
        return TubesDataLoader(
            ds,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=TubesDataModule.N_WORKERS,
            persistent_workers=TubesDataModule.N_WORKERS > 0
        )
    
    def test_dataloader(self):
        ds = TubesDataset(self.test_set, max_x=self.max_energy)
        return TubesDataLoader(
            ds,
            batch_size=self.batch_size,
            shuffle=False,
            num_workers=TubesDataModule.N_WORKERS,
            persistent_workers=TubesDataModule.N_WORKERS > 0
        )
