#!/usr/bin/env python3
"""Trim a ProteinDF pdf-archive-h5.py reference.h5 to the minimal regression set.

Keeps:
  - all top-level attributes
  - TEs, basisset/*, control/*, molecule/* (copied as-is)
  - occ/* (copied as-is)
  - energy_level/<runtype>_<N> and P/<runtype>_<N> for the converged
    iteration N (top-level attribute "iterations")

Drops everything else (per-iteration C, F, P, energy_level of non-converged
iterations, h, h2, s, ...), so the reference.h5 does not keep intermediate
SCF history that is not meant to be compared against.
"""

import argparse
import sys

import h5py

KEEP_AS_IS = ("TEs", "basisset", "control", "molecule", "occ")
KEEP_FINAL_ITERATION = ("energy_level", "P")


def trim(src_path, dst_path):
    with h5py.File(src_path, "r") as src, h5py.File(dst_path, "w") as dst:
        for key, value in src.attrs.items():
            dst.attrs[key] = value

        if "iterations" not in src.attrs:
            raise SystemExit(f"ERROR: {src_path}: top-level attribute 'iterations' is missing")
        iterations = int(src.attrs["iterations"])
        suffix = f"_{iterations}"

        for name in KEEP_AS_IS:
            if name not in src:
                raise SystemExit(f"ERROR: {src_path}: '{name}' is missing")
            src.copy(name, dst)

        for name in KEEP_FINAL_ITERATION:
            if name not in src:
                raise SystemExit(f"ERROR: {src_path}: '{name}' is missing")
            group = src[name]
            final_keys = [
                key for key in group.keys() if key.rsplit("_", 1)[-1] == str(iterations)
            ]
            if not final_keys:
                raise SystemExit(
                    f"ERROR: {src_path}: '{name}' has no dataset for the converged "
                    f"iteration ({suffix})"
                )
            dst_group = dst.create_group(name)
            for key in final_keys:
                src.copy(f"{name}/{key}", dst_group)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="input reference.h5 (full pdf-archive-h5.py output)")
    parser.add_argument("output", help="output reference.h5 (trimmed)")
    args = parser.parse_args()

    trim(args.input, args.output)


if __name__ == "__main__":
    sys.exit(main())
