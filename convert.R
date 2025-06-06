library(dplyr)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://www.ensembl.org")

genes <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
  values = rna_expression_combined$ensembl_gene_id,
  mart = mart
)

# Check for duplicated values in the entrez or gene symbol columns
genes %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    n_entrez = n_distinct(entrezgene_id, na.rm = TRUE),
    n_symbol = n_distinct(hgnc_symbol, na.rm = TRUE),
    .groups = "drop"
    ) %>%
  filter(n_entrez > 1 | n_symbol > 1)

# Remove duplicated ENSEMBL --> entrez mappings by retaining the shortest
# entrez id (most likely to not be a LOC)
genes_clean <- genes %>%
  mutate(entrez_length = nchar(as.character(entrezgene_id))) %>%
  group_by(ensembl_gene_id) %>%
  slice_min(entrez_length, with_ties = FALSE) %>%
  ungroup() %>%
  select(-entrez_length)  # drop helper column

# Check for duplicated values in the entrez or gene symbol columns
# (should hopefully be empty now)
genes_clean %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    n_entrez = n_distinct(entrezgene_id, na.rm = TRUE),
    n_symbol = n_distinct(hgnc_symbol, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_entrez > 1 | n_symbol > 1)

  # Join to add 'name' from lookup table into df1
rna_expression_combined <- rna_expression_combined %>%
  left_join(genes_clean, by = "ensembl_gene_id")
