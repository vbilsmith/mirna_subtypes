## ----Install-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#source("https://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("biomaRt", "consensusOV")
optional_bioc_packages <- c("AnnotationDbi", "org.Hs.eg.db")
missing_bioc_packages <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]
missing_optional_bioc_packages <- optional_bioc_packages[
  !vapply(optional_bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc_packages) > 0) {
  BiocManager::install(missing_bioc_packages)
}

if (length(missing_optional_bioc_packages) > 0) {
  try(
    BiocManager::install(missing_optional_bioc_packages, ask = FALSE, update = FALSE),
    silent = TRUE
  )
}


## ----Imports-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(biomaRt)
library(consensusOV)
library(dplyr)
library(genefu)
library(readxl)

## -- Paths --
subset_dir <- "data/gdc_tcga_ov_omics/downloads/mRNA"
konecny_file_path <- "data/jnci_JNCI_14_0249_s05.xls"


## ----Identify Files----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Find all files ending in .tsv recursively
tsv_files <- list.files(
  path = subset_dir,
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Filter to make sure we're only keeping files (not directories or symlinks)
tsv_files <- tsv_files[file.info(tsv_files)$isdir == FALSE]

print(length(tsv_files))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Select only files with rna_seq in name
rna_files <- tsv_files[grepl("rna_seq", tsv_files)]
print(length(rna_files))


## ----Create Table------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Loop through each file and extract the gene_id and tpm, which we can use to track all of the genes across all of the individuals

rna_list <- list()
for (file in rna_files) {
  temp_data <- read.delim(file, header = TRUE, sep = "\t", comment.char = "#")
  temp_data <- temp_data[, c("gene_id", "tpm_unstranded")] #gene_name, gene_type, raw read counts (for unstranded and stranded, maybe useful in normalization?) -- are less relevant in classification using consensusOV
  sample_name <- basename(file) 
  colnames(temp_data)[2] <- sample_name 
  rna_list[[sample_name]] <- temp_data
}
# This will create a table with gene_id as the first column and each sample
# in columns 2:
rna_expression_combined <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), rna_list)

# Create a column of cleaned ENSEMBL ids (remove version)
rna_expression_combined$ensembl_gene_id <- sub("\\.\\d+$", "",rna_expression_combined$gene_id)



## ----Map Ensembl IDs onto Entrez and Gene Symbol-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
map_ensembl_ids <- function(
  ensembl_ids,
  cache_file = "output/ensembl_gene_mapping.csv",
  use_biomart_fallback = TRUE
) {
  ensembl_ids <- unique(na.omit(ensembl_ids))
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  normalize_gene_map <- function(gene_map) {
    gene_map %>%
      mutate(
        ensembl_gene_id = as.character(ensembl_gene_id),
        entrezgene_id = as.character(entrezgene_id),
        hgnc_symbol = as.character(hgnc_symbol)
      )
  }

  if (file.exists(cache_file)) {
    cached_map <- normalize_gene_map(read.csv(cache_file, check.names = FALSE))
    cached_map <- cached_map[cached_map$ensembl_gene_id %in% ensembl_ids, ]
    cached_mapped_ids <- cached_map$ensembl_gene_id[!is.na(cached_map$entrezgene_id)]
    if (length(setdiff(ensembl_ids, cached_mapped_ids)) == 0) {
      message("Using cached Ensembl mapping: ", cache_file)
      return(cached_map)
    }
  } else {
    cached_map <- NULL
  }

  org_map <- data.frame(
    ensembl_gene_id = character(),
    entrezgene_id = character(),
    hgnc_symbol = character(),
    stringsAsFactors = FALSE
  )

  if (
    requireNamespace("AnnotationDbi", quietly = TRUE) &&
      requireNamespace("org.Hs.eg.db", quietly = TRUE)
  ) {
    message("Mapping Ensembl IDs with org.Hs.eg.db...")
    org_db <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))
    org_map <- AnnotationDbi::select(
      org_db,
      keys = ensembl_ids,
      keytype = "ENSEMBL",
      columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
    ) %>%
      rename(
        ensembl_gene_id = ENSEMBL,
        entrezgene_id = ENTREZID,
        hgnc_symbol = SYMBOL
      ) %>%
      normalize_gene_map()
  } else {
    message("org.Hs.eg.db is not available; using biomaRt for Ensembl mapping.")
  }

  missing_ids <- setdiff(ensembl_ids, org_map$ensembl_gene_id[!is.na(org_map$entrezgene_id)])

  if (use_biomart_fallback && length(missing_ids) > 0) {
    message("Mapping ", length(missing_ids), " IDs with biomaRt...")
    biomart_map <- tryCatch({
      mart <- biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = "hsapiens_gene_ensembl",
        host = "https://www.ensembl.org"
      )

      biomaRt::getBM(
        filters = "ensembl_gene_id",
        attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
        values = missing_ids,
        mart = mart
      )
    }, error = function(e) {
      warning("biomaRt fallback failed: ", conditionMessage(e))
      NULL
    })

    if (!is.null(biomart_map) && nrow(biomart_map) > 0) {
      biomart_map <- normalize_gene_map(biomart_map)
      org_map <- bind_rows(
        org_map,
        biomart_map %>%
          filter(!ensembl_gene_id %in% org_map$ensembl_gene_id[!is.na(org_map$entrezgene_id)])
      )
    }
  }

  if (!is.null(cached_map)) {
    org_map <- bind_rows(
      cached_map,
      org_map %>% filter(!ensembl_gene_id %in% cached_map$ensembl_gene_id)
    )
  }

  if (nrow(org_map) == 0) {
    stop("Could not map Ensembl IDs with org.Hs.eg.db or biomaRt.")
  }

  write.csv(org_map, cache_file, row.names = FALSE)
  org_map
}

genes <- map_ensembl_ids(rna_expression_combined$ensembl_gene_id)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Check for duplicated values in the entrez or gene symbol columns
genes %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    n_entrez = n_distinct(entrezgene_id, na.rm = TRUE),
    n_symbol = n_distinct(hgnc_symbol, na.rm = TRUE),
    .groups = "drop"
    ) %>%
  filter(n_entrez > 1 | n_symbol > 1)

bad <- genes %>%
  group_by(entrezgene_id) %>%
  summarise(
    n_entrez = n_distinct(ensembl_gene_id, na.rm = TRUE),
    .groups = "drop"
    ) %>%
  filter(n_entrez > 1)

# Remove duplicated ENSEMBL --> entrez mappings by retaining the shortest
# entrez id (most likely to be a real gene rather than a LOC)
genes_clean <- genes %>%
  # Step 1: Keep shortest entrez ID per ensembl ID
  mutate(entrez_length = nchar(as.character(entrezgene_id))) %>%
  group_by(ensembl_gene_id) %>%
  slice_min(entrez_length, with_ties = FALSE) %>%
  ungroup() %>%
  select(-entrez_length) %>%

  # Step 2: Remove duplicate entrez IDs if the symbol is blank or NA
  group_by(entrezgene_id) %>%
  filter(!(n() > 1 & (hgnc_symbol == "" | is.na(hgnc_symbol)))) %>%
  ungroup() %>%

  # Step 3: For remaining entrezgene_id duplicates, keep the shortest symbol, or first alphabetically if tied
  mutate(symbol_length = nchar(hgnc_symbol)) %>%
  group_by(entrezgene_id) %>%
  arrange(symbol_length, hgnc_symbol) %>%
  slice(1) %>%
  ungroup() %>%
  select(-symbol_length)

  # Join to add 'name' from lookup table into rna_expression_combined
rna_expression_harmonized <- rna_expression_combined %>%
  left_join(genes_clean, by = "ensembl_gene_id")

# Remove rows where entrez_ids are NA or match the duplicated IDs to remove
sample_columns <- names(rna_list)
first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  x[[1]]
}

rna_expression_filtered <- rna_expression_harmonized %>%
  filter(!is.na(entrezgene_id)) %>%
  mutate(entrezgene_id = as.character(entrezgene_id)) %>%
  group_by(entrezgene_id) %>%
  summarise(
    across(all_of(sample_columns), ~ mean(.x, na.rm = TRUE)),
    hgnc_symbol = first_non_missing(hgnc_symbol),
    .groups = "drop"
  )



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# filtering the data. we only need geneid and expression for each patients in wide format
data_matrix <- as.matrix(rna_expression_filtered[, sample_columns, drop = FALSE])


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Extract entrez_ids and assign to row names
rownames(data_matrix) <- rna_expression_filtered$entrezgene_id



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# The Konecky subtyping from ConsensusOV was not working, so we pull the 
# xpression values associated with each cluster from the original paper's 
# supplementary materials and use them instead



# Read sheet 4 (UCLA_Mayo_PAMlist)
konecny.supplementary.data <- read_excel(konecny_file_path, sheet = 4)

# Extract relevant columns (keep EntrezGeneID + 4 centroid columns)
konecny.centroids.raw <- konecny.supplementary.data[, c(2, 4:7)]

# Rename for clarity (optional)
colnames(konecny.centroids.raw)[1] <- "EntrezID"

# Convert EntrezID to character
konecny.centroids.raw$EntrezID <- as.character(konecny.centroids.raw$EntrezID)

# Average duplicates: one row per Entrez ID
konecny.centroids <- konecny.centroids.raw %>%
  group_by(EntrezID) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()

# Set rownames and remove EntrezID column
rownames(konecny.centroids) <- konecny.centroids$EntrezID
konecny.centroids <- unique(konecny.centroids)
konecny.centroids$EntrezID <- NULL

shared_genes <- intersect(rownames(konecny.centroids), rownames(data_matrix))



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
get.konecny.subtypes.fixed <- function(expression.matrix, entrez.ids) {
  
  # Z-score normalize each gene (row)
  expression.matrix <- t(scale(t(expression.matrix)))
  
  # Ensure entrez.ids is character, in case it's being used later
  entrez.ids <- as.character(entrez.ids)
  
  # Identify shared genes
  intersecting.entrez.ids <- intersect(rownames(expression.matrix), rownames(konecny.centroids))
  
  # Subset both matrices to shared genes
  expression.matrix <- expression.matrix[intersecting.entrez.ids, , drop = FALSE]
  konecny.centroids <- konecny.centroids[intersecting.entrez.ids, , drop = FALSE]
  
  # Sanity check
  if (!all(rownames(expression.matrix) == rownames(konecny.centroids))) {
    stop("There is a mismatch between the Entrez IDs in the reference and input datasets.")
  }

    # Coerce both to matrices
  konecny.centroids <- as.matrix(konecny.centroids)
  expression.matrix  <- as.matrix(expression.matrix)

  # Step 1: Find common genes
  shared_genes <- intersect(rownames(expression.matrix), rownames(konecny.centroids))
  
  # Drop genes with NA or constant rows across either matrix
  valid_genes <- which(
    complete.cases(konecny.centroids) &
    complete.cases(expression.matrix) &
    apply(konecny.centroids, 1, sd) > 0 &
    apply(expression.matrix, 1, sd) > 0
  )
  
  konecny.centroids <- konecny.centroids[valid_genes, , drop = FALSE]
  expression.matrix  <- expression.matrix[valid_genes, , drop = FALSE]
  
  # Now run cor()
  spearman.cc.vals <- cor(konecny.centroids, expression.matrix, method = "spearman")

  # Each row is a subtype, each column is a sample
  # So we want the subtype (row) with the max correlation for each column (sample)
  max.idx <- apply(spearman.cc.vals, 2, which.max)
  subtypes <- rownames(spearman.cc.vals)[max.idx]
  subtypes <- factor(subtypes, levels = rownames(spearman.cc.vals))

  return(list(Konecny.subtypes = subtypes, spearman.cc.vals = t(spearman.cc.vals)))
}



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bentink.subtypes <- get.subtypes(data_matrix, rownames(data_matrix), method = "Bentink")

konecny.subtypes <- get.konecny.subtypes.fixed(data_matrix, rownames(data_matrix))

helland.subtypes <- get.subtypes(data_matrix, rownames(data_matrix), method = "Helland")

verhaak.subtypes <- get.subtypes(data_matrix, rownames(data_matrix), method = "Verhaak")

conc.subtypes <- get.subtypes(data_matrix, rownames(data_matrix), "consensusOV")


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
table(bentink.subtypes$Bentink.subtypes)
table(helland.subtypes$Helland.subtypes)
table(konecny.subtypes$Konecny.subtypes)
table(verhaak.subtypes$Verhaak.subtypes)
table(conc.subtypes$consensusOV.subtypes)


## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dir.create("output", showWarnings = FALSE)

score_object_to_df <- function(score_object, assignments) {
  n_samples <- length(assignments)
  
  if (is.matrix(score_object) || is.data.frame(score_object)) {
    return(as.data.frame(score_object, check.names = FALSE))
  }
  
  if (is.atomic(score_object) && length(score_object) == n_samples) {
    return(data.frame(score = score_object, check.names = FALSE))
  }
  
  if (!is.list(score_object)) {
    stop("Score object is not a matrix, data frame, vector, or list.")
  }
  
  score_parts <- list()
  
  for (part_name in names(score_object)) {
    part <- score_object[[part_name]]
    
    if (is.matrix(part) || is.data.frame(part)) {
      part_df <- as.data.frame(part, check.names = FALSE)
      
      if (nrow(part_df) == n_samples) {
        score_parts[[part_name]] <- part_df
      } else if (ncol(part_df) == n_samples) {
        score_parts[[part_name]] <- as.data.frame(t(as.matrix(part_df)), check.names = FALSE)
      }
    } else if (is.atomic(part) && length(part) == n_samples) {
      score_parts[[part_name]] <- data.frame(part, check.names = FALSE)
      names(score_parts[[part_name]]) <- part_name
    }
  }
  
  if (length(score_parts) == 0) {
    stop("No score table elements matched the number of samples.")
  }
  
  do.call(cbind, score_parts)
}

make_subtype_export <- function(score_matrix, assignments) {
  score_df <- score_object_to_df(score_matrix, assignments)
  assignment_sample_ids <- names(assignments)
  
  if (
    !is.null(assignment_sample_ids) &&
    !is.null(colnames(score_df)) &&
    all(assignment_sample_ids %in% colnames(score_df)) &&
    !all(assignment_sample_ids %in% rownames(score_df))
  ) {
    score_df <- as.data.frame(t(as.matrix(score_df)), check.names = FALSE)
  }
  
  numeric_score_cols <- vapply(score_df, is.numeric, logical(1))
  
  sample_ids <- rownames(score_df)
  
  if (is.null(sample_ids)) {
    sample_ids <- seq_len(nrow(score_df))
  }
  
  assignment_values <- as.character(assignments)
  if (!is.null(assignment_sample_ids)) {
    names(assignment_values) <- assignment_sample_ids
  }
  
  if (!is.null(assignment_sample_ids) && all(sample_ids %in% assignment_sample_ids)) {
    assignment_values <- assignment_values[sample_ids]
  } else if (length(assignment_values) != length(sample_ids)) {
    stop("Assignments and score matrix have different sample counts.")
  }
  
  data.frame(
    sample = sample_ids,
    assigned_subtype = assignment_values,
    max_score = apply(score_df[, numeric_score_cols, drop = FALSE], 1, max, na.rm = TRUE),
    score_df,
    row.names = NULL,
    check.names = FALSE
  )
}

write_assignment_export <- function(assignments, file) {
  assignment_values <- as.character(assignments)
  sample_ids <- names(assignments)
  
  if (is.null(sample_ids)) {
    sample_ids <- seq_along(assignment_values)
  }
  
  write.csv(
    data.frame(
      sample = sample_ids,
      assigned_subtype = assignment_values,
      row.names = NULL,
      check.names = FALSE
    ),
    file,
    row.names = FALSE
  )
}

write_subtype_export <- function(score_matrix, assignments, scores_file, assignment_file = NULL) {
  if (!is.null(score_matrix)) {
    export_df <- make_subtype_export(score_matrix, assignments)
    write.csv(export_df, scores_file, row.names = FALSE)
    return(export_df)
  }
  
  if (is.null(assignment_file)) {
    stop("No score table available for ", scores_file)
  }
  
  write_assignment_export(assignments, assignment_file)
  invisible(NULL)
}

konecny_export <- make_subtype_export(
  konecny.subtypes$spearman.cc.vals,
  konecny.subtypes$Konecny.subtypes
)
write.csv(konecny_export, "output/konecny_subtypes_scores.csv", row.names = FALSE)

helland_export <- write_subtype_export(
  helland.subtypes$subtype.scores,
  helland.subtypes$Helland.subtypes,
  "output/helland_subtypes_scores.csv",
  "output/helland_subtypes.csv"
)

verhaak_export <- write_subtype_export(
  verhaak.subtypes$gsva.out,
  verhaak.subtypes$Verhaak.subtypes,
  "output/verhaak_subtypes_scores.csv"
)

if (!is.null(verhaak.subtypes$gsva.out)) {
  write.csv(
    as.data.frame(verhaak.subtypes$gsva.out, check.names = FALSE),
    "output/verhaak_gsva_scores_raw.csv",
    row.names = TRUE
  )
}

bentink_export <- write_subtype_export(
  bentink.subtypes$angio,
  bentink.subtypes$Bentink.subtypes,
  "output/bentink_subtypes_scores.csv"
)



## ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
consensus_export <- make_subtype_export(
  conc.subtypes$rf.probs,
  conc.subtypes$consensusOV.subtypes
)

write.csv(consensus_export, "output/consensusOV_subtypes_scores.csv", row.names = FALSE)


hist(
  consensus_export$max_score,
  breaks = 20,
  col = "steelblue",
  main = "Histogram of Max Consensus Scores per Row",
  xlab = "Max Consensus Score",
  ylab = "Count"
)
