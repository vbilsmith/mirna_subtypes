#!/usr/bin/env python3
"""Run the GDC TCGA-OV data acquisition pipeline."""

from __future__ import annotations

import argparse
import subprocess
import sys


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run TCGA-OV GDC omics data acquisition and optional downstream steps."
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="Download mRNA and miRNA expression files with gdc-client.",
    )
    parser.add_argument(
        "--skip-pull",
        action="store_true",
        help="Skip GDC queries and use existing files in --out-dir.",
    )
    parser.add_argument(
        "--harmonize-labels",
        action="store_true",
        help=(
            "After cluster labels have been generated, harmonize existing "
            "mRNA/miRNA labels to the current GDC sample tables."
        ),
    )
    parser.add_argument(
        "--mrna-label-dir",
        default="mRNA_clusters/output",
        help="Directory containing mRNA subtype score exports.",
    )
    parser.add_argument(
        "--mirna-label-dir",
        default=None,
        help=(
            "Optional directory containing miRNA ConsensusOV_labels.csv and "
            "ConsensusOV_probs.csv. Omit until miRNA clusters are available."
        ),
    )
    parser.add_argument(
        "--survival",
        action="store_true",
        help=(
            "Calculate survival fields and append them to harmonized labels. "
            "Use after --harmonize-labels, or after label outputs already exist."
        ),
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

    if not args.skip_pull:
        run(pull_command)

    if args.harmonize_labels:
        harmonize_command = [
            sys.executable,
            "scripts/harmonize_gdc_omics_labels.py",
            "--gdc-dir",
            args.out_dir,
            "--mrna-label-dir",
            args.mrna_label_dir,
        ]
        if args.mirna_label_dir:
            harmonize_command.extend(["--mirna-label-dir", args.mirna_label_dir])
        run(harmonize_command)

    if args.survival:
        run([
            sys.executable,
            "scripts/link_gdc_omics_samples_to_survival.py",
            "--gdc-dir",
            args.out_dir,
        ])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
