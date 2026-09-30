#!/usr/bin/env python3
"""Extract tube spectra from a RadField3D dataset's field metadata into .spectrum/.info files
usable as a ``source_spectra`` directory by create_dataset.py — e.g. to re-simulate with the
exact spectra family of an existing Monte-Carlo dataset (analytical pre-flight runs).
"""
import argparse
import glob
import hashlib
import json
import os
import random

import torch

from radfiled3d.metadata.v1 import Metadata
from radfiled3d.store import FieldStore


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True, help="Dataset root or fields directory with .rf3 files")
    parser.add_argument("--out", required=True, help="Output directory for .spectrum/.info files")
    parser.add_argument("--max_spectra", type=int, default=512, help="Number of files to sample spectra from")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    fields_dir = os.path.join(args.dataset, "fields") if os.path.isdir(os.path.join(args.dataset, "fields")) else args.dataset
    files = sorted(glob.glob(os.path.join(fields_dir, "*.rf3")))
    if not files:
        raise SystemExit(f"No .rf3 files found in {fields_dir}")
    random.Random(args.seed).shuffle(files)
    os.makedirs(args.out, exist_ok=True)

    written = 0
    seen = set()
    for path in files:
        if written >= args.max_spectra:
            break
        try:
            meta = Metadata.from_raw_metadata(FieldStore.load_metadata(path))
            table = meta.simulation.tube.spectrum
        except Exception:
            continue
        if table is None or (table[:, 1] <= 0).all():
            continue
        digest = hashlib.sha1(table.tobytes()).hexdigest()[:16]
        if digest in seen:
            continue
        seen.add(digest)
        # Trim trailing zero-fluence bins: SpectrumLoader's PDF/CDF otherwise smears the last real
        # bin's probability mass linearly up to the highest listed energy and reports it as max().
        last = int((table[:, 1] > 0).nonzero()[0].max())
        table = table[:last + 2] if last + 2 <= len(table) else table
        energies_kev = torch.tensor(table[:, 0] / 1000.0, dtype=torch.float64)
        fluence = torch.tensor(table[:, 1], dtype=torch.float64)
        torch.save(torch.stack([energies_kev, fluence]), os.path.join(args.out, f"{digest}.spectrum"))
        max_energy_ev = float(table[table[:, 1] > 0][:, 0].max())
        with open(os.path.join(args.out, f"{digest}.info"), "w") as f:
            json.dump({"energy": max_energy_ev, "source_field": os.path.basename(path)}, f)
        written += 1

    print(f"Extracted {written} unique spectra to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
