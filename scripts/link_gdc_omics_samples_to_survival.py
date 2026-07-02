#!/usr/bin/env python3
"""Link GDC mRNA/miRNA sample manifests to clinical survival fields."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


DEFAULT_GDC_DIR = "data/gdc_tcga_ov_omics"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Join mRNA and miRNA GDC file/sample link tables to case-level "
            "clinical survival information from the GDC clinical JSON."
        )
    )
    parser.add_argument(
        "--gdc-dir",
        default=DEFAULT_GDC_DIR,
        help="Directory containing outputs from pull_gdc_tcga_ov_omics.py.",
    )
    return parser.parse_args()


def numeric_or_none(value: Any) -> float | None:
    if value in (None, "", "'--", "--", "not reported", "Not Reported"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def integer_string(value: float | None) -> str:
    if value is None:
        return ""
    if float(value).is_integer():
        return str(int(value))
    return str(value)


def first_non_missing(values: list[Any]) -> str:
    for value in values:
        if value not in (None, "", "'--", "--", "not reported", "Not Reported"):
            return str(value)
    return ""


def max_numeric(values: list[Any]) -> float | None:
    numeric_values = [numeric_or_none(value) for value in values]
    numeric_values = [value for value in numeric_values if value is not None]
    return max(numeric_values) if numeric_values else None


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def case_survival_rows(clinical_json: Path) -> list[dict[str, str]]:
    with clinical_json.open() as handle:
        cases = json.load(handle)

    rows: list[dict[str, str]] = []
    for case in cases:
        demographic = case.get("demographic") or {}
        diagnoses = case.get("diagnoses") or []
        follow_ups = case.get("follow_ups") or []

        diagnosis_follow_up_days = [
            diagnosis.get("days_to_last_follow_up") for diagnosis in diagnoses
        ]
        follow_up_days = [follow_up.get("days_to_follow_up") for follow_up in follow_ups]
        days_to_last_follow_up = max_numeric(diagnosis_follow_up_days + follow_up_days)

        days_to_death = numeric_or_none(demographic.get("days_to_death"))
        vital_status = str(demographic.get("vital_status") or "")
        os_event = "1" if vital_status.lower() == "dead" else "0"

        if os_event == "1" and days_to_death is not None:
            os_days = days_to_death
            os_time_source = "days_to_death"
        elif days_to_last_follow_up is not None:
            os_days = days_to_last_follow_up
            os_time_source = "days_to_last_follow_up"
        else:
            os_days = None
            os_time_source = ""

        primary_diagnosis = first_non_missing([
            diagnosis.get("primary_diagnosis") for diagnosis in diagnoses
        ])
        figo_stage = first_non_missing([
            diagnosis.get("figo_stage") for diagnosis in diagnoses
        ])
        tumor_grade = first_non_missing([
            diagnosis.get("tumor_grade") for diagnosis in diagnoses
        ])
        year_of_diagnosis = first_non_missing([
            diagnosis.get("year_of_diagnosis") for diagnosis in diagnoses
        ])
        days_to_diagnosis = first_non_missing([
            diagnosis.get("days_to_diagnosis") for diagnosis in diagnoses
        ])
        progression_or_recurrence = first_non_missing([
            diagnosis.get("progression_or_recurrence") for diagnosis in diagnoses
        ])
        prior_malignancy = first_non_missing([
            diagnosis.get("prior_malignancy") for diagnosis in diagnoses
        ])
        prior_treatment = first_non_missing([
            diagnosis.get("prior_treatment") for diagnosis in diagnoses
        ])

        rows.append({
            "case_id": case.get("case_id", ""),
            "case_submitter_id": case.get("submitter_id", ""),
            "vital_status": vital_status,
            "os_event": os_event,
            "os_days": integer_string(os_days),
            "os_time_source": os_time_source,
            "days_to_death": integer_string(days_to_death),
            "days_to_last_follow_up": integer_string(days_to_last_follow_up),
            "age_at_index": str(demographic.get("age_at_index") or ""),
            "days_to_birth": integer_string(numeric_or_none(demographic.get("days_to_birth"))),
            "primary_diagnosis": primary_diagnosis,
            "figo_stage": figo_stage,
            "tumor_grade": tumor_grade,
            "year_of_diagnosis": year_of_diagnosis,
            "days_to_diagnosis": days_to_diagnosis,
            "progression_or_recurrence": progression_or_recurrence,
            "prior_malignancy": prior_malignancy,
            "prior_treatment": prior_treatment,
            "n_diagnosis_records": str(len(diagnoses)),
            "n_follow_up_records": str(len(follow_ups)),
        })

    return sorted(rows, key=lambda row: row["case_submitter_id"])


def link_assay_samples(
    assay: str,
    sample_rows: list[dict[str, str]],
    survival_by_case: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    linked: list[dict[str, str]] = []
    for row in sample_rows:
        survival = survival_by_case.get(row["case_id"], {})
        linked_row = {
            "assay": assay,
            **row,
            **{key: value for key, value in survival.items() if key not in {"case_id", "case_submitter_id"}},
            "has_survival": "TRUE" if survival else "FALSE",
            "has_os_time": "TRUE" if survival.get("os_days") else "FALSE",
        }
        linked.append(linked_row)

    return sorted(linked, key=lambda row: (
        row.get("case_submitter_id", ""),
        row.get("sample_type_code", ""),
        row.get("sample_submitter_id", ""),
        row.get("assay", ""),
        row.get("file_name", ""),
    ))


def sample_availability_rows(long_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], dict[str, Any]] = {}
    for row in long_rows:
        key = (row["case_id"], row["sample_submitter_id"])
        current = grouped.setdefault(key, {
            "case_id": row["case_id"],
            "case_submitter_id": row["case_submitter_id"],
            "sample_submitter_id": row["sample_submitter_id"],
            "sample_type": row["sample_type"],
            "sample_type_code": row["sample_type_code"],
            "mRNA_files": [],
            "miRNA_files": [],
            "vital_status": row.get("vital_status", ""),
            "os_event": row.get("os_event", ""),
            "os_days": row.get("os_days", ""),
            "days_to_last_follow_up": row.get("days_to_last_follow_up", ""),
            "figo_stage": row.get("figo_stage", ""),
            "tumor_grade": row.get("tumor_grade", ""),
        })
        current[f"{row['assay']}_files"].append(row["file_name"])

    availability: list[dict[str, str]] = []
    for row in grouped.values():
        mrna_files = sorted(row["mRNA_files"])
        mirna_files = sorted(row["miRNA_files"])
        availability.append({
            "case_id": row["case_id"],
            "case_submitter_id": row["case_submitter_id"],
            "sample_submitter_id": row["sample_submitter_id"],
            "sample_type": row["sample_type"],
            "sample_type_code": row["sample_type_code"],
            "has_mRNA": "TRUE" if mrna_files else "FALSE",
            "has_miRNA": "TRUE" if mirna_files else "FALSE",
            "n_mRNA_files": str(len(mrna_files)),
            "n_miRNA_files": str(len(mirna_files)),
            "mRNA_files": ";".join(mrna_files),
            "miRNA_files": ";".join(mirna_files),
            "vital_status": row["vital_status"],
            "os_event": row["os_event"],
            "os_days": row["os_days"],
            "days_to_last_follow_up": row["days_to_last_follow_up"],
            "figo_stage": row["figo_stage"],
            "tumor_grade": row["tumor_grade"],
        })

    return sorted(availability, key=lambda row: (
        row["case_submitter_id"],
        row["sample_type_code"],
        row["sample_submitter_id"],
    ))


def main() -> int:
    args = parse_args()
    gdc_dir = Path(args.gdc_dir)

    survival_rows = case_survival_rows(gdc_dir / "clinical_union_current_gdc.json")
    survival_by_case = {row["case_id"]: row for row in survival_rows}

    mrna_rows = read_tsv(gdc_dir / "mrna_file_sample_link.tsv")
    mirna_rows = read_tsv(gdc_dir / "mirna_file_sample_link.tsv")

    linked_mrna = link_assay_samples("mRNA", mrna_rows, survival_by_case)
    linked_mirna = link_assay_samples("miRNA", mirna_rows, survival_by_case)
    linked_long = sorted(
        linked_mrna + linked_mirna,
        key=lambda row: (
            row.get("case_submitter_id", ""),
            row.get("sample_type_code", ""),
            row.get("sample_submitter_id", ""),
            row.get("assay", ""),
            row.get("file_name", ""),
        )
    )
    availability = sample_availability_rows(linked_long)

    sample_fields = [
        "assay",
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "state",
        "data_type",
        "experimental_strategy",
        "data_format",
        "case_id",
        "case_submitter_id",
        "sample_id",
        "sample_submitter_id",
        "sample_type",
        "sample_type_code",
        "aliquot_barcode",
        "vital_status",
        "os_event",
        "os_days",
        "os_time_source",
        "days_to_death",
        "days_to_last_follow_up",
        "age_at_index",
        "days_to_birth",
        "primary_diagnosis",
        "figo_stage",
        "tumor_grade",
        "year_of_diagnosis",
        "days_to_diagnosis",
        "progression_or_recurrence",
        "prior_malignancy",
        "prior_treatment",
        "n_diagnosis_records",
        "n_follow_up_records",
        "has_survival",
        "has_os_time",
    ]
    survival_fields = [
        "case_id",
        "case_submitter_id",
        "vital_status",
        "os_event",
        "os_days",
        "os_time_source",
        "days_to_death",
        "days_to_last_follow_up",
        "age_at_index",
        "days_to_birth",
        "primary_diagnosis",
        "figo_stage",
        "tumor_grade",
        "year_of_diagnosis",
        "days_to_diagnosis",
        "progression_or_recurrence",
        "prior_malignancy",
        "prior_treatment",
        "n_diagnosis_records",
        "n_follow_up_records",
    ]
    availability_fields = [
        "case_id",
        "case_submitter_id",
        "sample_submitter_id",
        "sample_type",
        "sample_type_code",
        "has_mRNA",
        "has_miRNA",
        "n_mRNA_files",
        "n_miRNA_files",
        "mRNA_files",
        "miRNA_files",
        "vital_status",
        "os_event",
        "os_days",
        "days_to_last_follow_up",
        "figo_stage",
        "tumor_grade",
    ]

    write_tsv(gdc_dir / "clinical_survival_case_level.tsv", survival_rows, survival_fields)
    write_tsv(gdc_dir / "mrna_samples_survival.tsv", linked_mrna, sample_fields)
    write_tsv(gdc_dir / "mirna_samples_survival.tsv", linked_mirna, sample_fields)
    write_tsv(gdc_dir / "omics_samples_survival_long.tsv", linked_long, sample_fields)
    write_tsv(gdc_dir / "omics_sample_assay_availability.tsv", availability, availability_fields)

    write_csv(gdc_dir / "clinical_survival_case_level.csv", survival_rows, survival_fields)
    write_csv(gdc_dir / "mrna_samples_survival.csv", linked_mrna, sample_fields)
    write_csv(gdc_dir / "mirna_samples_survival.csv", linked_mirna, sample_fields)
    write_csv(gdc_dir / "omics_samples_survival_long.csv", linked_long, sample_fields)
    write_csv(gdc_dir / "omics_sample_assay_availability.csv", availability, availability_fields)

    missing_mrna = sum(1 for row in linked_mrna if row["has_survival"] == "FALSE")
    missing_mirna = sum(1 for row in linked_mirna if row["has_survival"] == "FALSE")
    missing_mrna_os_time = sum(1 for row in linked_mrna if row["has_os_time"] == "FALSE")
    missing_mirna_os_time = sum(1 for row in linked_mirna if row["has_os_time"] == "FALSE")
    paired_samples = sum(1 for row in availability if row["has_mRNA"] == "TRUE" and row["has_miRNA"] == "TRUE")

    print(f"Case-level survival rows: {len(survival_rows)}")
    print(f"mRNA sample/file rows linked to survival: {len(linked_mrna)}")
    print(f"miRNA sample/file rows linked to survival: {len(linked_mirna)}")
    print(f"mRNA rows missing survival: {missing_mrna}")
    print(f"miRNA rows missing survival: {missing_mirna}")
    print(f"mRNA rows missing OS time: {missing_mrna_os_time}")
    print(f"miRNA rows missing OS time: {missing_mirna_os_time}")
    print(f"Unique assay/sample availability rows: {len(availability)}")
    print(f"Samples with both mRNA and miRNA files: {paired_samples}")
    print(f"Outputs written to: {gdc_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
