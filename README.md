# Ovarian Cancer Subtyping Analysis

The manifest of TCGA-OV samples is located in: 
```text
data/tcga-ov-metadata/metadata.cohort.2025-06-23.json
```

That JSON is used to define the TCGA-OV case cohort. 
The pipeline then queries GDC for current mRNA, miRNA, and clinical data.

## Requirements

- [gdc client](https://anaconda.org/channels/bioconda/packages/gdc-client/overview)

## Pipeline Steps

Here are the steps of the full workflow:

1. `scripts/pull_gdc_tcga_ov_omics.py`
   - queries GDC for current TCGA-OV mRNA STAR-count files
   - queries GDC for current TCGA-OV mature miRNA quantification files
   - pulls clinical records for the union of mRNA/miRNA cases
   - writes manifests and file-to-sample link tables

2. Characterize mRNA and miRNA clusters
   - run the cluster assignment workflows
   - write cluster/subtype assignments before survival is calculated

3. `scripts/harmonize_gdc_omics_labels.py`
   - joins existing mRNA subtype labels onto current GDC mRNA rows
   - joins miRNA NMF subtype labels onto current GDC miRNA rows
   - writes harmonized label tables and subtype agreement summaries

4. `scripts/link_gdc_omics_samples_to_survival.py`
   - computes case-level overall survival fields
   - joins survival fields onto mRNA and miRNA sample/file rows
   - appends survival fields to the harmonized label tables


## Downloading and Harmonizing Data and Metadata

To identify the mRNA and miRNA data associated with the TCGA-OV cases and download the corresponding mRNA and miRNA data, run:

```bash
python3 scripts/run_gdc_omics_pipeline.py --download
```

This uses `gdc-client` with the generated manifests and writes outputs under:

```text
data/gdc_tcga_ov_omics/
```

### Key Outputs

```text
data/gdc_tcga_ov_omics/gdc_manifest_mrna_star_counts.tsv
data/gdc_tcga_ov_omics/gdc_manifest_mirna_mirbase21_mature.tsv
data/gdc_tcga_ov_omics/mrna_file_sample_link.tsv
data/gdc_tcga_ov_omics/mirna_file_sample_link.tsv
data/gdc_tcga_ov_omics/harmonized_labels/mrna_samples_labels.tsv
data/gdc_tcga_ov_omics/harmonized_labels/mirna_samples_labels.tsv
data/gdc_tcga_ov_omics/harmonized_labels/omics_sample_assay_availability_labels.tsv
data/gdc_tcga_ov_omics/harmonized_labels/omics_samples_survival_labels_long.tsv
```

## Identify mRNA Clusters

The first analysis generates mRNA clusters from the TCGA-OV data using the package ConsensusOV.
This analysis runs out of the `mRNA_clusters` directory.
Run the `AssignClusters.Rmd` notebook or `AssignClusters.R` script to generate cluster assignments and scores for each sample.
These outputs will be written to the `output` directory within `mRNA_clusters`

## Identify miRNA Clusters

The miRNA clustering workflow uses the miRNA files downloaded by the shared GDC
pipeline:

```text
data/gdc_tcga_ov_omics/downloads/miRNA/
data/gdc_tcga_ov_omics/mirna_file_sample_link.tsv
```

Run the miRNA notebooks in order:

```bash
quarto render miRNA_clusters/01_mirna_preprocess.qmd
quarto render miRNA_clusters/02a_mirna_pca.qmd
quarto render miRNA_clusters/02b_mirna_kmeans.qmd
quarto render miRNA_clusters/02c_mirna_nmf.qmd
```

The preprocessing notebook writes the normalized handoff matrix and GDC sample
metadata to:

```text
data/mirna_data/mirna_lcpm_matrix.rds
data/mirna_data/mirna_sample_metadata.csv
```

The NMF notebook writes the miRNA subtype assignments to:

```text
miRNA_clusters/mirna_nmf_subtypes_k4.csv
```

That CSV is the default miRNA subtype input for harmonization. It should contain
one row per clustered miRNA file/sample and a `cluster` column; the current
pipeline also carries through GDC identifiers such as `file_name`,
`case_submitter_id`, `sample_submitter_id`, and `aliquot_barcode`.

## Harmonize Labels With Clinical Data

After mRNA and miRNA cluster assignments have been generated, run:

```bash
python3 scripts/run_gdc_omics_pipeline.py --skip-pull --harmonize-labels --survival
```

This harmonizes mRNA subtype labels and miRNA NMF subtype labels, then appends
clinical/survival fields. By default, miRNA labels are read from:

```text
miRNA_clusters/mirna_nmf_subtypes_k4.csv
```

To use a different miRNA subtype file, pass it explicitly:

```bash
python3 scripts/run_gdc_omics_pipeline.py --skip-pull --harmonize-labels --survival --mirna-subtype-file path/to/mirna_nmf_subtypes_k4.csv
```

For survival analyses, filter to:

```text
has_os_time == TRUE
```
