# Ovarian Cancer Subtyping Analysis

## mRNA Clusters

The first analysis generates mRNA clusters from the TCGA-OV data using the package ConsensusOV.
Here we use a manifest of TCGA-OV samples filtered for those that contain RNA-Seq data, corresponding to 424 samples.
The manifest is located in `gdc_manifest_2025-03-30.tsv`.

To download the data associated with this manifest, use the command:
```
mkdir -p gdc_download_20250330_115711.237016

gdc-client download \
  -m gdc_manifest_reconstructed_2025-03-30.tsv \
  -d gdc_download_20250330_115711.237016
```

Then run the `AssignClusters.Rmd` notebook to generate cluster assignments and scores for each sample following ConsensusOV.
These outputs will be written to the `output` directory within `mRNA_clusters`

