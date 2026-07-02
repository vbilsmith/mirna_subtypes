#!/usr/bin/env python3
"""Harmonize subtype labels with GDC mRNA/miRNA sample-survival tables."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


DEFAULT_GDC_DIR = "data/gdc_tcga_ov_omics"
DEFAULT_MRNA_LABEL_DIR = "mRNA_clusters/output"
DEFAULT_MIRNA_LABEL_DIR = "miRNA_clusters"


MRNA_SUBTYPE_FILES = {
    "consensusOV": "consensusOV_subtypes_scores.csv",
    "konecny": "konecny_subtypes_scores.csv",
    "helland": "helland_subtypes_scores.csv",
    "verhaak": "verhaak_subtypes_scores.csv",
    "bentink": "bentink_subtypes_scores.csv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Join existing mRNA/miRNA subtype labels onto the current GDC "
            "sample-survival tables produced by the pipeline."
        )
    )
    parser.add_argument("--gdc-dir", default=DEFAULT_GDC_DIR)
    parser.add_argument("--mrna-label-dir", default=DEFAULT_MRNA_LABEL_DIR)
    parser.add_argument("--mirna-label-dir", default=DEFAULT_MIRNA_LABEL_DIR)
    parser.add_argument(
        "--out-dir",
        default=None,
        help="Output directory. Defaults to <gdc-dir>/harmonized_labels.",
    )
    return parser.parse_args()


def read_delimited(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_csv(path: Path) -> list[dict[str, str]]:
    return read_delimited(path, ",")


def read_tsv(path: Path) -> list[dict[str, str]]:
    return read_delimited(path, "\t")


def write_delimited(path: Path, rows: list[dict[str, Any]], delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_both(out_dir: Path, stem: str, rows: list[dict[str, Any]]) -> None:
    write_delimited(out_dir / f"{stem}.tsv", rows, "\t")
    write_delimited(out_dir / f"{stem}.csv", rows, ",")


def prefix_mrna_subtype_rows(rows: list[dict[str, str]], method: str) -> dict[str, dict[str, str]]:
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        sample = row.get("sample", "")
        if not sample:
            continue

        prefixed: dict[str, str] = {}
        for key, value in row.items():
            if key == "sample":
                continue
            if key == "assigned_subtype":
                prefixed[f"{method}_subtype"] = value
            elif key == "max_score":
                prefixed[f"{method}_max_score"] = value
            else:
                prefixed[f"{method}_{key}"] = value
        indexed[sample] = prefixed
    return indexed


def load_mrna_labels(label_dir: Path) -> dict[str, dict[str, str]]:
    merged: dict[str, dict[str, str]] = defaultdict(dict)
    for method, file_name in MRNA_SUBTYPE_FILES.items():
        path = label_dir / file_name
        if not path.exists():
            continue
        for sample, labels in prefix_mrna_subtype_rows(read_csv(path), method).items():
            merged[sample].update(labels)
    return dict(merged)


def load_mirna_consensus_labels(label_dir: Path) -> dict[str, dict[str, str]]:
    labels_path = label_dir / "ConsensusOV_labels.csv"
    probs_path = label_dir / "ConsensusOV_probs.csv"
    if not labels_path.exists():
        return {}

    labels: dict[str, dict[str, str]] = {}
    for row in read_csv(labels_path):
        sample = row.get("sample", "")
        if not sample:
            continue
        labels[sample] = {
            "miRNA_consensusOV_subtype": row.get("type", ""),
            "miRNA_consensusOV_max_score": row.get("score", ""),
        }

    if probs_path.exists():
        for row in read_csv(probs_path):
            sample = row.get("") or row.get("sample") or ""
            if not sample:
                continue
            labels.setdefault(sample, {})
            for key, value in row.items():
                if key in {"", "sample"}:
                    continue
                labels[sample][f"miRNA_consensusOV_{key}"] = value

    return labels


def add_labels_by_file_name(
    rows: list[dict[str, str]],
    labels_by_file: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], int]:
    linked: list[dict[str, str]] = []
    matched = 0
    for row in rows:
        labels = labels_by_file.get(row.get("file_name", ""), {})
        if labels:
            matched += 1
        linked.append({**row, **labels})
    return linked, matched


def resolve_mirna_label_samples(
    mirna_labels_by_sample: dict[str, dict[str, str]],
    mrna_rows: list[dict[str, str]],
    mirna_rows: list[dict[str, str]],
) -> tuple[dict[str, dict[str, str]], dict[str, Any]]:
    mrna_by_file = {row["file_name"]: row for row in mrna_rows}
    mirna_by_file = {row["file_name"]: row for row in mirna_rows}

    labels_by_sample_submitter: dict[str, dict[str, str]] = {}
    source_counts: Counter[str] = Counter()
    unresolved: list[str] = []

    for label_sample, labels in mirna_labels_by_sample.items():
        labels_with_source = dict(labels)
        if label_sample in mirna_by_file:
            sample_id = mirna_by_file[label_sample]["sample_submitter_id"]
            labels_with_source["miRNA_label_source"] = "miRNA_file_name"
            labels_by_sample_submitter[sample_id] = labels_with_source
            source_counts["miRNA_file_name"] += 1
        elif label_sample in mrna_by_file:
            sample_id = mrna_by_file[label_sample]["sample_submitter_id"]
            labels_with_source["miRNA_label_source"] = "mRNA_file_name"
            labels_with_source["miRNA_label_source_file_name"] = label_sample
            labels_by_sample_submitter[sample_id] = labels_with_source
            source_counts["mRNA_file_name"] += 1
        else:
            unresolved.append(label_sample)
            source_counts["unresolved"] += 1

    diagnostics = {
        "source_counts": dict(source_counts),
        "unresolved": unresolved,
    }
    return labels_by_sample_submitter, diagnostics


def add_labels_by_sample_submitter(
    rows: list[dict[str, str]],
    labels_by_sample: dict[str, dict[str, str]],
) -> tuple[list[dict[str, str]], int]:
    linked: list[dict[str, str]] = []
    matched = 0
    for row in rows:
        labels = labels_by_sample.get(row.get("sample_submitter_id", ""), {})
        if labels:
            matched += 1
        linked.append({**row, **labels})
    return linked, matched


def collapse_case_level(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    preferred = sorted(
        rows,
        key=lambda row: (
            row.get("case_submitter_id", ""),
            0 if row.get("sample_type_code") == "01" else 1,
            row.get("sample_submitter_id", ""),
            row.get("file_name", ""),
        )
    )
    seen: set[str] = set()
    collapsed: list[dict[str, str]] = []
    for row in preferred:
        case_id = row.get("case_id", "")
        if case_id in seen:
            continue
        seen.add(case_id)
        collapsed.append(row)
    return collapsed


def make_pairwise_confusion(rows: list[dict[str, str]], subtype_cols: dict[str, str]) -> list[dict[str, Any]]:
    methods = list(subtype_cols)
    out: list[dict[str, Any]] = []
    for i, row_method in enumerate(methods):
        for column_method in methods[i + 1:]:
            row_col = subtype_cols[row_method]
            column_col = subtype_cols[column_method]
            counts: Counter[tuple[str, str]] = Counter()
            for row in rows:
                row_subtype = row.get(row_col, "")
                column_subtype = row.get(column_col, "")
                if row_subtype and column_subtype:
                    counts[(row_subtype, column_subtype)] += 1
            for (row_subtype, column_subtype), n in sorted(counts.items()):
                out.append({
                    "row_method": row_method,
                    "column_method": column_method,
                    "row_subtype": row_subtype,
                    "column_subtype": column_subtype,
                    "n": n,
                })
    return out


def availability_with_labels(
    availability_rows: list[dict[str, str]],
    mrna_by_sample: dict[str, dict[str, str]],
    mirna_by_sample: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in availability_rows:
        sample = row.get("sample_submitter_id", "")
        mrna_labels = {
            key: value
            for key, value in mrna_by_sample.get(sample, {}).items()
            if key.endswith("_subtype") or key.endswith("_max_score")
        }
        mirna_labels = mirna_by_sample.get(sample, {})
        rows.append({**row, **mrna_labels, **mirna_labels})
    return rows


def main() -> int:
    args = parse_args()
    gdc_dir = Path(args.gdc_dir)
    out_dir = Path(args.out_dir) if args.out_dir else gdc_dir / "harmonized_labels"

    mrna_rows = read_tsv(gdc_dir / "mrna_samples_survival.tsv")
    mirna_rows = read_tsv(gdc_dir / "mirna_samples_survival.tsv")
    availability_rows = read_tsv(gdc_dir / "omics_sample_assay_availability.tsv")

    mrna_labels_by_file = load_mrna_labels(Path(args.mrna_label_dir))
    mrna_labeled, mrna_label_matches = add_labels_by_file_name(mrna_rows, mrna_labels_by_file)
    mrna_labels_by_sample = {
        row["sample_submitter_id"]: {
            key: value for key, value in row.items()
            if key.endswith("_subtype") or key.endswith("_max_score")
        }
        for row in mrna_labeled
        if row.get("sample_submitter_id")
    }

    mirna_labels_raw = load_mirna_consensus_labels(Path(args.mirna_label_dir))
    mirna_labels_by_sample, mirna_diagnostics = resolve_mirna_label_samples(
        mirna_labels_raw,
        mrna_rows,
        mirna_rows,
    )
    mirna_labeled, mirna_label_matches = add_labels_by_sample_submitter(
        mirna_rows,
        mirna_labels_by_sample,
    )

    long_labeled = sorted(
        mrna_labeled + mirna_labeled,
        key=lambda row: (
            row.get("case_submitter_id", ""),
            row.get("sample_type_code", ""),
            row.get("sample_submitter_id", ""),
            row.get("assay", ""),
            row.get("file_name", ""),
        )
    )
    availability_labeled = availability_with_labels(
        availability_rows,
        mrna_labels_by_sample,
        mirna_labels_by_sample,
    )
    mrna_case_level = collapse_case_level(mrna_labeled)
    mirna_case_level = collapse_case_level(mirna_labeled)

    subtype_cols = {
        "mRNA_consensusOV": "consensusOV_subtype",
        "mRNA_konecny": "konecny_subtype",
        "mRNA_helland": "helland_subtype",
        "mRNA_verhaak": "verhaak_subtype",
        "mRNA_bentink": "bentink_subtype",
        "miRNA_consensusOV": "miRNA_consensusOV_subtype",
    }
    sample_confusion = make_pairwise_confusion(availability_labeled, subtype_cols)

    write_both(out_dir, "mrna_samples_survival_labels", mrna_labeled)
    write_both(out_dir, "mirna_samples_survival_labels", mirna_labeled)
    write_both(out_dir, "omics_samples_survival_labels_long", long_labeled)
    write_both(out_dir, "omics_sample_assay_availability_labels", availability_labeled)
    write_both(out_dir, "mrna_case_level_survival_labels", mrna_case_level)
    write_both(out_dir, "mirna_case_level_survival_labels", mirna_case_level)
    write_both(out_dir, "subtype_confusion_long", sample_confusion)

    diagnostics = [
        {"metric": "mRNA_sample_rows", "value": len(mrna_rows)},
        {"metric": "mRNA_rows_with_any_label", "value": mrna_label_matches},
        {"metric": "miRNA_sample_rows", "value": len(mirna_rows)},
        {"metric": "miRNA_rows_with_consensusOV_label", "value": mirna_label_matches},
        {"metric": "miRNA_label_input_rows", "value": len(mirna_labels_raw)},
        {
            "metric": "miRNA_label_resolved_from_mRNA_file_name",
            "value": mirna_diagnostics["source_counts"].get("mRNA_file_name", 0),
        },
        {
            "metric": "miRNA_label_resolved_from_miRNA_file_name",
            "value": mirna_diagnostics["source_counts"].get("miRNA_file_name", 0),
        },
        {
            "metric": "miRNA_label_unresolved",
            "value": mirna_diagnostics["source_counts"].get("unresolved", 0),
        },
        {"metric": "availability_rows", "value": len(availability_labeled)},
    ]
    write_both(out_dir, "label_harmonization_diagnostics", diagnostics)

    if mirna_diagnostics["unresolved"]:
        write_both(
            out_dir,
            "unresolved_miRNA_label_samples",
            [{"sample": sample} for sample in mirna_diagnostics["unresolved"]],
        )

    print(f"mRNA rows: {len(mrna_rows)}")
    print(f"mRNA rows with any subtype labels: {mrna_label_matches}")
    print(f"miRNA rows: {len(mirna_rows)}")
    print(f"miRNA rows with ConsensusOV labels: {mirna_label_matches}")
    print(f"miRNA label resolution: {mirna_diagnostics['source_counts']}")
    print(f"Outputs written to: {out_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
