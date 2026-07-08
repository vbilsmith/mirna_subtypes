#!/usr/bin/env python3
"""Pull current TCGA-OV mRNA, miRNA, and clinical data from GDC.

The only local input is a GDC metadata JSON used to define the TCGA-OV case
cohort. The script then queries GDC for current open-access mRNA and miRNA
quantification files for those cases, writes separate manifests, fetches
clinical data for the union of mRNA/miRNA cases, and optionally downloads the
files with gdc-client.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any


GDC_API = "https://api.gdc.cancer.gov"
DEFAULT_SEED_METADATA = "data/tcga-ov-metadata/metadata.cohort.2025-06-23.json"
DEFAULT_OUT_DIR = "data/gdc_tcga_ov_omics"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Query GDC for TCGA-OV mRNA STAR count files, miRNA mature-miRNA "
            "quantification files, and clinical records for the union of cases."
        )
    )
    parser.add_argument(
        "--seed-metadata",
        default=DEFAULT_SEED_METADATA,
        help="Local GDC metadata JSON used only to define the seed case cohort.",
    )
    parser.add_argument(
        "--out-dir",
        default=DEFAULT_OUT_DIR,
        help="Output directory for manifests, metadata, clinical JSON, and downloads.",
    )
    parser.add_argument(
        "--download",
        action="store_true",
        help="After writing manifests, run gdc-client download for mRNA and miRNA files.",
    )
    parser.add_argument(
        "--gdc-client",
        default="gdc-client",
        help="Path to gdc-client executable used when --download is set.",
    )
    parser.add_argument(
        "--api-url",
        default=GDC_API,
        help="GDC API base URL.",
    )
    parser.add_argument(
        "--page-size",
        type=int,
        default=1000,
        help="Number of GDC records requested per API page.",
    )
    return parser.parse_args()


def read_seed_case_ids(seed_metadata: Path) -> list[str]:
    with seed_metadata.open() as handle:
        records = json.load(handle)

    case_ids: set[str] = set()
    for record in records:
        for entity in record.get("associated_entities") or []:
            case_id = entity.get("case_id")
            if case_id:
                case_ids.add(case_id)

    if not case_ids:
        raise ValueError(f"No case IDs found in {seed_metadata}")

    return sorted(case_ids)


def op_in(field: str, values: list[str]) -> dict[str, Any]:
    return {"op": "in", "content": {"field": field, "value": values}}


def op_eq(field: str, value: str) -> dict[str, Any]:
    return op_in(field, [value])


def op_and(*filters: dict[str, Any]) -> dict[str, Any]:
    return {"op": "and", "content": list(filters)}


def gdc_post(
    api_url: str,
    endpoint: str,
    params: dict[str, str | int],
    payload: dict[str, Any],
    retries: int = 3,
) -> dict[str, Any]:
    url = f"{api_url.rstrip('/')}/{endpoint.lstrip('/')}?{urllib.parse.urlencode(params)}"
    data = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=data,
        headers={"Content-Type": "application/json"},
        method="POST",
    )

    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            if attempt == retries:
                raise RuntimeError(f"GDC HTTP {exc.code} for {url}\n{body}") from exc
        except urllib.error.URLError as exc:
            if attempt == retries:
                raise RuntimeError(f"GDC request failed for {url}: {exc}") from exc
        time.sleep(2 * attempt)

    raise RuntimeError(f"GDC request failed for {url}")


def fetch_all(
    api_url: str,
    endpoint: str,
    filters: dict[str, Any],
    fields: list[str],
    page_size: int,
    expand: list[str] | None = None,
) -> list[dict[str, Any]]:
    params: dict[str, str | int] = {
        "format": "JSON",
        "size": page_size,
        "from": 0,
    }
    if fields:
        params["fields"] = ",".join(fields)
    if expand:
        params["expand"] = ",".join(expand)

    records: list[dict[str, Any]] = []
    total: int | None = None

    while total is None or len(records) < total:
        params["from"] = len(records)
        response = gdc_post(api_url, endpoint, params, {"filters": filters})
        data = response.get("data") or {}
        hits = data.get("hits") or []
        records.extend(hits)
        total = int((data.get("pagination") or {}).get("total") or len(records))
        if not hits:
            break

    return records


def file_filter(case_ids: list[str], data_type: str, strategy: str) -> dict[str, Any]:
    return op_and(
        op_eq("cases.project.project_id", "TCGA-OV"),
        op_in("cases.case_id", case_ids),
        op_eq("access", "open"),
        op_eq("state", "released"),
        op_eq("data_category", "Transcriptome Profiling"),
        op_eq("data_type", data_type),
        op_eq("experimental_strategy", strategy),
    )


def first_case(record: dict[str, Any]) -> dict[str, Any]:
    cases = record.get("cases") or []
    return cases[0] if cases else {}


def first_sample(record: dict[str, Any]) -> dict[str, Any]:
    samples = first_case(record).get("samples") or []
    return samples[0] if samples else {}


def collect_aliquot_barcodes(node: Any) -> list[str]:
    barcodes: list[str] = []
    if isinstance(node, dict):
        submitter_id = node.get("submitter_id")
        if isinstance(submitter_id, str) and submitter_id.startswith("TCGA-"):
            if len(submitter_id.split("-")) >= 7:
                barcodes.append(submitter_id)
        for value in node.values():
            barcodes.extend(collect_aliquot_barcodes(value))
    elif isinstance(node, list):
        for value in node:
            barcodes.extend(collect_aliquot_barcodes(value))
    return sorted(set(barcodes))


def join_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return ";".join(str(item) for item in value if item is not None)
    return str(value)


def tcga_sample_type_code(sample_submitter_id: str) -> str:
    if len(sample_submitter_id) >= 15 and sample_submitter_id.startswith("TCGA-"):
        return sample_submitter_id[13:15]
    return ""


def simplify_file_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for record in records:
        case = first_case(record)
        sample = first_sample(record)
        aliquot_barcodes = collect_aliquot_barcodes(case)
        aliquot_barcode = aliquot_barcodes[0] if aliquot_barcodes else ""
        sample_submitter_id = sample.get("submitter_id", "")
        rows.append(
            {
                "file_id": record.get("file_id", ""),
                "file_name": record.get("file_name", ""),
                "md5sum": record.get("md5sum", ""),
                "file_size": record.get("file_size", ""),
                "state": record.get("state", ""),
                "data_type": record.get("data_type", ""),
                "experimental_strategy": record.get("experimental_strategy", ""),
                "data_format": record.get("data_format", ""),
                "case_id": case.get("case_id", ""),
                "case_submitter_id": case.get("submitter_id", ""),
                "sample_id": sample.get("sample_id", ""),
                "sample_submitter_id": sample_submitter_id,
                "sample_type": sample.get("sample_type", ""),
                "sample_type_code": tcga_sample_type_code(sample_submitter_id),
                "aliquot_barcode": aliquot_barcode,
            }
        )
    return sorted(rows, key=lambda row: (row["case_submitter_id"], row["file_name"]))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(path: Path, records: list[dict[str, Any]]) -> None:
    rows = [
        {
            "id": record.get("file_id", ""),
            "filename": record.get("file_name", ""),
            "md5": record.get("md5sum", ""),
            "size": record.get("file_size", ""),
            "state": record.get("state", ""),
        }
        for record in sorted(records, key=lambda row: row.get("file_name", ""))
    ]
    write_tsv(path, rows, ["id", "filename", "md5", "size", "state"])


def case_summary_rows(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case in sorted(cases, key=lambda item: item.get("submitter_id", "")):
        demographic = case.get("demographic") or {}
        diagnoses = case.get("diagnoses") or []
        diagnosis = diagnoses[0] if diagnoses else {}
        rows.append(
            {
                "case_id": case.get("case_id", ""),
                "case_submitter_id": case.get("submitter_id", ""),
                "project_id": (case.get("project") or {}).get("project_id", ""),
                "primary_site": join_value(case.get("primary_site")),
                "disease_type": join_value(case.get("disease_type")),
                "vital_status": demographic.get("vital_status", ""),
                "days_to_death": demographic.get("days_to_death", ""),
                "days_to_birth": demographic.get("days_to_birth", ""),
                "age_at_index": demographic.get("age_at_index", ""),
                "year_of_birth": demographic.get("year_of_birth", ""),
                "year_of_death": demographic.get("year_of_death", ""),
                "primary_diagnosis": diagnosis.get("primary_diagnosis", ""),
                "figo_stage": diagnosis.get("figo_stage", ""),
                "tumor_grade": diagnosis.get("tumor_grade", ""),
                "year_of_diagnosis": diagnosis.get("year_of_diagnosis", ""),
                "days_to_diagnosis": diagnosis.get("days_to_diagnosis", ""),
                "progression_or_recurrence": diagnosis.get("progression_or_recurrence", ""),
                "prior_malignancy": diagnosis.get("prior_malignancy", ""),
                "prior_treatment": diagnosis.get("prior_treatment", ""),
            }
        )
    return rows


def download_manifest(gdc_client: str, manifest: Path, destination: Path) -> None:
    if not shutil.which(gdc_client) and not Path(gdc_client).exists():
        raise FileNotFoundError(f"Could not find gdc-client executable: {gdc_client}")
    destination.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [gdc_client, "download", "-m", str(manifest), "-d", str(destination)],
        check=True,
    )


def main() -> int:
    args = parse_args()
    seed_metadata = Path(args.seed_metadata)
    out_dir = Path(args.out_dir)

    seed_case_ids = read_seed_case_ids(seed_metadata)
    print(f"Seed cases from {seed_metadata}: {len(seed_case_ids)}")

    file_fields = [
        "file_id",
        "file_name",
        "md5sum",
        "file_size",
        "state",
        "data_type",
        "experimental_strategy",
        "data_format",
        "cases.case_id",
        "cases.submitter_id",
        "cases.samples.sample_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.tumor_descriptor",
        "cases.samples.portions.analytes.aliquots.submitter_id",
    ]
    file_expand = [
        "cases",
        "cases.samples",
        "cases.samples.portions",
        "cases.samples.portions.analytes",
        "cases.samples.portions.analytes.aliquots",
    ]

    mrna_records = fetch_all(
        args.api_url,
        "files",
        file_filter(seed_case_ids, "Gene Expression Quantification", "RNA-Seq"),
        file_fields,
        args.page_size,
        file_expand,
    )
    mrna_records = [
        record
        for record in mrna_records
        if str(record.get("file_name", "")).endswith(".rna_seq.augmented_star_gene_counts.tsv")
    ]

    mirna_records = fetch_all(
        args.api_url,
        "files",
        file_filter(seed_case_ids, "miRNA Expression Quantification", "miRNA-Seq"),
        file_fields,
        args.page_size,
        file_expand,
    )
    mirna_records = [
        record
        for record in mirna_records
        if str(record.get("file_name", "")).endswith(".mirbase21.mirnas.quantification.txt")
    ]

    mrna_cases = {first_case(record).get("case_id") for record in mrna_records}
    mirna_cases = {first_case(record).get("case_id") for record in mirna_records}
    union_case_ids = sorted(case_id for case_id in (mrna_cases | mirna_cases) if case_id)

    clinical_fields = [
        "case_id",
        "submitter_id",
        "project.project_id",
        "primary_site",
        "disease_type",
        "demographic",
        "diagnoses",
        "treatments",
        "exposures",
        "follow_ups",
    ]
    clinical_expand = [
        "project",
        "demographic",
        "diagnoses",
        "diagnoses.treatments",
        "exposures",
        "follow_ups",
    ]
    clinical_cases = fetch_all(
        args.api_url,
        "cases",
        op_and(op_eq("project.project_id", "TCGA-OV"), op_in("case_id", union_case_ids)),
        clinical_fields,
        args.page_size,
        clinical_expand,
    )

    write_json(out_dir / "metadata_mrna_current_gdc.json", mrna_records)
    write_json(out_dir / "metadata_mirna_current_gdc.json", mirna_records)
    write_json(out_dir / "clinical_union_current_gdc.json", clinical_cases)
    write_manifest(out_dir / "gdc_manifest_mrna_star_counts.tsv", mrna_records)
    write_manifest(out_dir / "gdc_manifest_mirna_mirbase21_mature.tsv", mirna_records)
    write_tsv(out_dir / "mrna_file_sample_link.tsv", simplify_file_records(mrna_records))
    write_tsv(out_dir / "mirna_file_sample_link.tsv", simplify_file_records(mirna_records))
    write_tsv(
        out_dir / "clinical_union_case_summary.tsv",
        case_summary_rows(clinical_cases),
    )
    write_tsv(
        out_dir / "union_case_ids.tsv",
        [{"case_id": case_id} for case_id in union_case_ids],
        ["case_id"],
    )

    print(f"mRNA STAR count files: {len(mrna_records)}")
    print(f"miRNA mature quantification files: {len(mirna_records)}")
    print(f"mRNA cases: {len(mrna_cases)}")
    print(f"miRNA cases: {len(mirna_cases)}")
    print(f"Union cases for clinical pull: {len(union_case_ids)}")
    print(f"Clinical cases returned: {len(clinical_cases)}")
    print(f"Outputs written to: {out_dir}")

    if args.download:
        download_manifest(
            args.gdc_client,
            out_dir / "gdc_manifest_mrna_star_counts.tsv",
            out_dir / "downloads" / "mRNA",
        )
        download_manifest(
            args.gdc_client,
            out_dir / "gdc_manifest_mirna_mirbase21_mature.tsv",
            out_dir / "downloads" / "miRNA",
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
