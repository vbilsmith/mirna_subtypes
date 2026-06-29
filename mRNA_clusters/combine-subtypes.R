library(dplyr)
library(jsonlite)
library(stringr)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != "'--" & x != "--" & x != ""]
  if (length(x) == 0) NA else x[[1]]
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

annotation_summary <- function(annotations) {
  if (length(annotations) == 0) {
    return(data.frame(
      qc_flag = FALSE,
      qc_dnu = FALSE,
      qc_center_failed = FALSE,
      qc_low_coverage = FALSE,
      qc_high_mitochondrial = FALSE,
      noncanonical_flag = FALSE,
      ovarian_triplet_flag = FALSE,
      annotation_categories = NA_character_,
      annotation_notes = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  categories <- vapply(
    annotations,
    function(annotation) annotation$category %||% NA_character_,
    character(1)
  )
  notes <- vapply(
    annotations,
    function(annotation) annotation$notes %||% NA_character_,
    character(1)
  )
  annotation_text <- paste(categories, notes, collapse = " ")

  qc_dnu <- any(str_detect(categories, regex("Item flagged DNU", ignore_case = TRUE)))
  qc_center_failed <- any(str_detect(categories, regex("Center QC failed", ignore_case = TRUE)))
  qc_low_coverage <- any(str_detect(notes, regex("LOW 5/3 COVERAGE RATIO", ignore_case = TRUE)))
  qc_high_mitochondrial <- any(str_detect(notes, regex("HIGH MITOCHONDRIAL CONTENT", ignore_case = TRUE)))

  data.frame(
    qc_flag = qc_dnu || qc_center_failed || qc_low_coverage || qc_high_mitochondrial,
    qc_dnu = qc_dnu,
    qc_center_failed = qc_center_failed,
    qc_low_coverage = qc_low_coverage,
    qc_high_mitochondrial = qc_high_mitochondrial,
    noncanonical_flag = str_detect(annotation_text, regex("noncanonical", ignore_case = TRUE)),
    ovarian_triplet_flag = str_detect(annotation_text, regex("ovarian? triplet", ignore_case = TRUE)),
    annotation_categories = collapse_unique(categories),
    annotation_notes = collapse_unique(notes),
    stringsAsFactors = FALSE
  )
}

prefix_score_columns <- function(df, method) {
  keep <- c("sample", "assigned_subtype", "max_score")
  score_cols <- setdiff(names(df), keep)
  names(df)[names(df) == "assigned_subtype"] <- paste0(method, "_subtype")
  names(df)[names(df) == "max_score"] <- paste0(method, "_max_score")
  names(df)[names(df) %in% score_cols] <- paste0(method, "_", score_cols)
  df
}

read_subtype_export <- function(file, method) {
  read.csv(file, check.names = FALSE) |>
    prefix_score_columns(method)
}

make_pairwise_confusion <- function(df, subtype_cols) {
  method_pairs <- combn(names(subtype_cols), 2, simplify = FALSE)

  bind_rows(lapply(method_pairs, function(method_pair) {
    row_method <- method_pair[[1]]
    column_method <- method_pair[[2]]
    row_col <- subtype_cols[[row_method]]
    column_col <- subtype_cols[[column_method]]

    df |>
      filter(!is.na(.data[[row_col]]), !is.na(.data[[column_col]])) |>
      count(
        row_subtype = .data[[row_col]],
        column_subtype = .data[[column_col]],
        name = "n"
      ) |>
      mutate(
        row_method = row_method,
        column_method = column_method,
        .before = row_subtype
      )
  }))
}

write_pairwise_confusion_matrices <- function(confusion_df, level, output_dir) {
  split_keys <- paste(confusion_df$row_method, confusion_df$column_method, sep = "_vs_")

  invisible(lapply(split(confusion_df, split_keys), function(pair_df) {
    matrix_df <- as.data.frame.matrix(xtabs(n ~ row_subtype + column_subtype, data = pair_df))
    matrix_df <- cbind(row_subtype = row.names(matrix_df), matrix_df)
    row.names(matrix_df) <- NULL

    file_name <- paste0(
      "rna_",
      level,
      "_confusion_",
      pair_df$row_method[[1]],
      "_vs_",
      pair_df$column_method[[1]],
      ".csv"
    )

    write.csv(matrix_df, file.path(output_dir, file_name), row.names = FALSE)
  }))
}

metadata_file <- "data/metadata.cart.2025-03-30.json"
clinical_file <- "data/clinical.cart.2025-04-08/clinical.tsv"
output_dir <- "mRNA_clusters/output"

metadata_raw <- fromJSON(metadata_file, simplifyVector = FALSE)

rna_metadata <- bind_rows(lapply(metadata_raw, function(record) {
  entity <- record$associated_entities[[1]]
  aliquot_barcode <- entity$entity_submitter_id %||% NA_character_
  annotations <- annotation_summary(record$annotations %||% list())

  cbind(data.frame(
    file_id = record$file_id %||% NA_character_,
    file_name = record$file_name %||% NA_character_,
    case_id = entity$case_id %||% NA_character_,
    aliquot_id = entity$entity_id %||% NA_character_,
    aliquot_barcode = aliquot_barcode,
    sample_barcode = str_sub(aliquot_barcode, 1, 16),
    case_submitter_id = str_sub(aliquot_barcode, 1, 12),
    sample_type_code = str_sub(aliquot_barcode, 14, 15),
    stringsAsFactors = FALSE
  ), annotations)
}))

clinical <- read.delim(
  clinical_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  na.strings = c("'--", "--", "")
)

# The GDC clinical TSV has multiple rows per case when diagnoses/treatments repeat.
# Collapse to one row per case for downstream joins.
clinical_cases <- clinical |>
  group_by(`cases.case_id`, `cases.submitter_id`) |>
  summarise(
    clinical_rows = n(),
    vital_status = first_non_missing(`demographic.vital_status`),
    days_to_death = first_non_missing(`demographic.days_to_death`),
    days_to_last_follow_up = first_non_missing(`diagnoses.days_to_last_follow_up`),
    age_at_index = first_non_missing(`demographic.age_at_index`),
    days_to_birth = first_non_missing(`demographic.days_to_birth`),
    figo_stage = first_non_missing(`diagnoses.figo_stage`),
    tumor_grade = first_non_missing(`diagnoses.tumor_grade`),
    primary_diagnosis = first_non_missing(`diagnoses.primary_diagnosis`),
    year_of_diagnosis = first_non_missing(`diagnoses.year_of_diagnosis`),
    progression_or_recurrence = first_non_missing(`diagnoses.progression_or_recurrence`),
    prior_malignancy = first_non_missing(`diagnoses.prior_malignancy`),
    prior_treatment = first_non_missing(`diagnoses.prior_treatment`),
    .groups = "drop"
  )

subtype_files <- list(
  consensusOV = file.path(output_dir, "consensusOV_subtypes_scores.csv"),
  konecny = file.path(output_dir, "konecny_subtypes_scores.csv"),
  helland = file.path(output_dir, "helland_subtypes_scores.csv"),
  verhaak = file.path(output_dir, "verhaak_subtypes_scores.csv"),
  bentink = file.path(output_dir, "bentink_subtypes_scores.csv")
)

subtype_tables <- Map(read_subtype_export, subtype_files, names(subtype_files))
subtypes <- Reduce(function(x, y) full_join(x, y, by = "sample"), subtype_tables)

sample_level <- subtypes |>
  left_join(rna_metadata, by = c("sample" = "file_name")) |>
  left_join(
    clinical_cases,
    by = c(
      "case_id" = "cases.case_id",
      "case_submitter_id" = "cases.submitter_id"
    )
  ) |>
  relocate(
    sample,
    file_id,
    aliquot_barcode,
    sample_barcode,
    sample_type_code,
    case_submitter_id,
    case_id
  )

missing_metadata <- sample_level |> filter(is.na(case_id))
missing_clinical <- sample_level |>
  filter(!is.na(case_id), is.na(vital_status), is.na(primary_diagnosis))

if (nrow(missing_metadata) > 0) {
  warning(nrow(missing_metadata), " RNA subtype rows did not match metadata.")
}

if (nrow(missing_clinical) > 0) {
  warning(nrow(missing_clinical), " RNA subtype rows matched metadata but not clinical.")
}

# Keep all aliquots in the sample-level output. For case-level analyses, prefer
# primary tumor samples (01) and then use a deterministic aliquot/file ordering.
case_level <- sample_level |>
  mutate(sample_type_priority = if_else(sample_type_code == "01", 0L, 1L)) |>
  arrange(case_submitter_id, sample_type_priority, aliquot_barcode, sample) |>
  group_by(case_submitter_id) |>
  slice(1) |>
  ungroup() |>
  select(-sample_type_priority)

subtype_cols <- c(
  consensusOV = "consensusOV_subtype",
  konecny = "konecny_subtype",
  helland = "helland_subtype",
  verhaak = "verhaak_subtype",
  bentink = "bentink_subtype"
)

sample_confusion <- make_pairwise_confusion(sample_level, subtype_cols)
case_confusion <- make_pairwise_confusion(case_level, subtype_cols)

write.csv(
  sample_level,
  file.path(output_dir, "rna_subtypes_clinical_sample_level.csv"),
  row.names = FALSE
)

write.csv(
  case_level,
  file.path(output_dir, "rna_subtypes_clinical_case_level.csv"),
  row.names = FALSE
)

write.csv(
  sample_confusion,
  file.path(output_dir, "rna_sample_level_subtype_confusion_long.csv"),
  row.names = FALSE
)

write.csv(
  case_confusion,
  file.path(output_dir, "rna_case_level_subtype_confusion_long.csv"),
  row.names = FALSE
)

write_pairwise_confusion_matrices(sample_confusion, "sample_level", output_dir)
write_pairwise_confusion_matrices(case_confusion, "case_level", output_dir)

duplicate_cases <- sample_level |>
  count(case_submitter_id, sort = TRUE) |>
  filter(n > 1)

duplicate_case_details <- sample_level |>
  semi_join(duplicate_cases, by = "case_submitter_id") |>
  arrange(case_submitter_id, sample_type_code, aliquot_barcode, sample) |>
  select(
    case_submitter_id,
    sample_type_code,
    sample_barcode,
    aliquot_barcode,
    sample,
    consensusOV_subtype,
    konecny_subtype,
    helland_subtype,
    verhaak_subtype,
    bentink_subtype,
    qc_flag,
    qc_dnu,
    qc_center_failed,
    qc_low_coverage,
    qc_high_mitochondrial,
    noncanonical_flag,
    ovarian_triplet_flag,
    annotation_categories,
    annotation_notes,
    vital_status,
    figo_stage,
    tumor_grade
  )

low_qc_samples <- sample_level |>
  filter(qc_flag) |>
  arrange(case_submitter_id, sample_type_code, aliquot_barcode, sample) |>
  select(
    case_submitter_id,
    sample_type_code,
    sample_barcode,
    aliquot_barcode,
    sample,
    qc_flag,
    qc_dnu,
    qc_center_failed,
    qc_low_coverage,
    qc_high_mitochondrial,
    annotation_categories,
    annotation_notes,
    consensusOV_subtype,
    konecny_subtype,
    helland_subtype,
    verhaak_subtype,
    bentink_subtype
  )

duplicate_case_comparison <- duplicate_case_details |>
  group_by(case_submitter_id) |>
  summarise(
    n_samples = n(),
    sample_types = paste(sample_type_code, collapse = ";"),
    consensusOV_calls = paste(unique(consensusOV_subtype), collapse = ";"),
    consensusOV_concordant = n_distinct(consensusOV_subtype) == 1,
    konecny_calls = paste(unique(konecny_subtype), collapse = ";"),
    konecny_concordant = n_distinct(konecny_subtype) == 1,
    helland_calls = paste(unique(helland_subtype), collapse = ";"),
    helland_concordant = n_distinct(helland_subtype) == 1,
    verhaak_calls = paste(unique(verhaak_subtype), collapse = ";"),
    verhaak_concordant = n_distinct(verhaak_subtype) == 1,
    bentink_calls = paste(unique(bentink_subtype), collapse = ";"),
    bentink_concordant = n_distinct(bentink_subtype) == 1,
    all_methods_concordant = consensusOV_concordant &
      konecny_concordant &
      helland_concordant &
      verhaak_concordant &
      bentink_concordant,
    .groups = "drop"
  )

write.csv(
  duplicate_cases,
  file.path(output_dir, "rna_duplicate_cases.csv"),
  row.names = FALSE
)

write.csv(
  duplicate_case_details,
  file.path(output_dir, "rna_duplicate_cases_subtype_details.csv"),
  row.names = FALSE
)

write.csv(
  duplicate_case_comparison,
  file.path(output_dir, "rna_duplicate_cases_subtype_comparison.csv"),
  row.names = FALSE
)

write.csv(
  low_qc_samples,
  file.path(output_dir, "rna_low_qc_samples.csv"),
  row.names = FALSE
)

message("Sample-level rows: ", nrow(sample_level))
message("Case-level rows: ", nrow(case_level))
message("Rows missing metadata: ", nrow(missing_metadata))
message("Rows missing clinical: ", nrow(missing_clinical))
message("Rows flagged low QC: ", nrow(low_qc_samples))
message("Cases with multiple RNA files: ", nrow(duplicate_cases))
message(
  "Duplicate cases concordant across all subtype methods: ",
  sum(duplicate_case_comparison$all_methods_concordant)
)
