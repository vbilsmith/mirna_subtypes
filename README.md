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
   - maps existing miRNA ConsensusOV labels by sample barcode
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

After cluster assignments have been generated, run:

```bash
python3 scripts/run_gdc_omics_pipeline.py --skip-pull --harmonize-labels --survival
```


