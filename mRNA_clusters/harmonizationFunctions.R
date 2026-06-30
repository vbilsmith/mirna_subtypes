`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

first_non_missing <- function(x) {
  x <- x[!is.na(x) & x != "'--" & x != "--" & x != ""]
  if (length(x) == 0) NA else x[[1]]
}

# Combines unique elements into a ;-separated string, or returns NA if
# no elements are non-empty
collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}

# Parse the annotations. First it sets a default blank records, then it
# extracts relevant annotations from the metadata and fills in as appropriate
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

# Streamline addition of multiple subtype columns by including the name of the
# subtyping paradigm used for each call
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