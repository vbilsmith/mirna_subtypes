clinical <- read.delim("data/clinical.cart.2025-04-08/clinical.tsv", header = TRUE, sep = "\t", comment.char = "#")
clinical <- clinical[, !sapply(clinical, function(col) all(col == "'--"))]

# Identify which columns are fixed across all rows
is_fixed <- sapply(clinical, function(col) length(unique(col)) == 1)

# Split the dataframe
clinical_fixed <- clinical[, is_fixed, drop = FALSE]
clinical_var   <- clinical[, !is_fixed, drop = FALSE]


library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)

### Step 1: Extract sample IDs from rownames / names for all subtype systems

# Helland
helland_ids <- rownames(helland.subtypes$subtype.scores) %>%
  str_extract("^[^\\.]+")
helland_df <- data.frame(
  case_id = helland_ids,
  Helland = helland.subtypes$Helland.subtypes,
  Helland_score = helland.subtypes$subtype.scores
)

# Konecny
konecny_ids <- rownames(konecny.subtypes$spearman.cc.vals) %>%
  str_extract("^[^\\.]+")
konecny_df <- data.frame(
  case_id = konecny_ids,
  Konecny = konecny.subtypes$Konecny.subtypes,
  Konecny_corr = konecny.subtypes$spearman.cc.vals
)

# Bentink
bentink_ids <- names(bentink.subtypes$angio$score) %>%
  str_extract("^[^\\.]+")
bentink_df <- data.frame(
  case_id = bentink_ids,
  Bentink = bentink.subtypes$angio$subtype
)

# Verhaak
verhaak_ids <- row.names(verhaak.subtypes$gsva.out) %>%
  str_extract("^[^\\.]+")
verhaak_df <- data.frame(
  case_id = verhaak_ids,
  Verhaak = verhaak.subtypes$Verhaak.subtypes,
  Verhaak_gsva = verhaak.subtypes$gsva.out
)

# Consensus
consensus_ids <- row.names(conc.subtypes$rf.probs) %>%
  str_extract("^[^\\.]+")
conc_df <- data.frame(
  case_id = consensus_ids,
  ConsensusOV = conc.subtypes$consensusOV.subtypes,
  ConsensusOV_probs = conc.subtypes$rf.probs
)


### Step 2: Combine all subtype assignments into one data frame

subtype_df <- reduce(
  list(helland_df, bentink_df, verhaak_df, konecny_df, conc_df),
  full_join,
  by = "case_id"
)

write.csv(subtype_df, "rna/subtypes.csv")


# Possible clinical variables of interest:
# diagnoses.primary_diagnosis
# diagnoses.tumor_grade
# diagnoses.year_of_diagnosis
# diagnoses.tissue_or_organ_of_origin
# demographic.days_to_birth
# demogragphipc.days_to_death



