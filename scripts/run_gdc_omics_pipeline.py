#!/usr/bin/env python3
"""Run the GDC TCGA-OV data pull, survival join, and label harmonization."""

from __future__ import annotations

import argparse
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the TCGA-OV GDC omics harmonization pipeline."
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download mRNA and miRNA expression files with gdc-client.",
    )
    parser.add_argument(
        "--out-dir",
        default="data/gdc_tcga_ov_omics",
        help="Pipeline output directory.",
    )
    parser.add_argument(
        "--seed-metadata",
        default="data/tcga-ov-metadata/metadata.cohort.2025-06-23.json",
        help="Seed metadata JSON defining the starting TCGA-OV cohort.",
    )
    return parser.parse_args()


def run(command: list[str]) -> None:
    print("Running:", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> int:
    args = parse_args()

    pull_command = [
        sys.executable,
        "scripts/pull_gdc_tcga_ov_omics.py",
        "--seed-metadata",
        args.seed_metadata,
        "--out-dir",
        args.out_dir,
    ]
    if args.download:
        pull_command.append("--download")

    run(pull_command)
    run([
        sys.executable,
        "scripts/link_gdc_omics_samples_to_survival.py",
        "--gdc-dir",
        args.out_dir,
    ])
    run([
        sys.executable,
        "scripts/harmonize_gdc_omics_labels.py",
        "--gdc-dir",
        args.out_dir,
    ])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
