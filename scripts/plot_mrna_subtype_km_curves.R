#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(survival)
  library(survminer)
})

input_file <- "mRNA_clusters/output/rna_subtypes_clinical_sample_level.csv"
plot_dir <- "mRNA_clusters/mrna_plots"
summary_file <- file.path(plot_dir, "mrna_subtype_km_logrank_summary.csv")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c(
  consensusOV = "consensusOV_subtype",
  konecny     = "konecny_subtype",
  helland     = "helland_subtype",
  verhaak     = "verhaak_subtype",
  bentink     = "bentink_subtype"
)

method_titles <- c(
  consensusOV = "consensusOV",
  konecny     = "Konecny",
  helland     = "Helland",
  verhaak     = "Verhaak",
  bentink     = "Bentink"
)

poster_palette <- c(
  "#4F6FA8",
  "#7F9CCB",
  "#7A7F87",
  "#D1837D",
  "#A83F4D"
)

poster_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

make_safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

fmt_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  }
  if (p < 0.001) {
    return("<0.001")
  }
  formatC(p, format = "f", digits = 3)
}

mrna <- read_csv(input_file, show_col_types = FALSE)

surv_df <- mrna |>
  filter(!is.na(case_submitter_id)) |>
  mutate(is_primary = suppressWarnings(as.integer(sample_type_code)) == 1L) |>
  arrange(case_submitter_id, desc(is_primary)) |>
  group_by(case_submitter_id) |>
  slice(1) |>
  ungroup() |>
  mutate(
    os_event = if_else(vital_status == "Dead", 1L, 0L),
    os_time = if_else(
      vital_status == "Dead",
      as.numeric(days_to_death),
      as.numeric(days_to_last_follow_up)
    )
  ) |>
  filter(!is.na(os_time), os_time >= 0, !is.na(os_event))

summary_rows <- list()

for (method_name in names(methods)) {
  subtype_col <- methods[[method_name]]

  d <- surv_df |>
    filter(!is.na(.data[[subtype_col]]), .data[[subtype_col]] != "") |>
    mutate(subtype = factor(.data[[subtype_col]], levels = sort(unique(.data[[subtype_col]]))))

  if (n_distinct(d$subtype) < 2L) {
    warning("Skipping ", method_name, ": fewer than two subtype groups.")
    next
  }

  km_fit <- survfit(Surv(os_time, os_event) ~ subtype, data = d)
  logrank <- survdiff(Surv(os_time, os_event) ~ subtype, data = d)
  logrank_p <- pchisq(logrank$chisq, length(logrank$n) - 1L, lower.tail = FALSE)

  group_counts <- d |>
    count(subtype, name = "n") |>
    mutate(label = paste0(subtype, " (n=", n, ")"))

  palette <- poster_palette[seq_len(n_distinct(d$subtype))]

  km_plot <- ggsurvplot(
    km_fit,
    data = d,
    palette = palette,
    pval = paste0("Log-rank p = ", fmt_p(logrank_p)),
    risk.table = TRUE,
    conf.int = FALSE,
    censor.shape = "|",
    censor.size = 2.4,
    title = paste0("Overall survival by mRNA subtype: ", method_titles[[method_name]]),
    subtitle = paste0(nrow(d), " cases, ", sum(d$os_event), " deaths"),
    xlab = "Days",
    ylab = "Overall survival probability",
    legend.title = method_titles[[method_name]],
    legend.labs = as.character(levels(d$subtype)),
    ggtheme = poster_theme,
    tables.theme = theme_minimal(base_size = 10)
  )

  output_file <- file.path(
    plot_dir,
    paste0("poster_mrna_km_survival_", make_safe_name(method_name), ".png")
  )

  png(output_file, width = 1100, height = 900, res = 130)
  print(km_plot)
  dev.off()

  summary_rows[[method_name]] <- tibble(
    method = method_titles[[method_name]],
    subtype_column = subtype_col,
    n_cases = nrow(d),
    n_events = sum(d$os_event),
    n_groups = n_distinct(d$subtype),
    logrank_p = logrank_p,
    logrank_p_label = fmt_p(logrank_p),
    group_counts = paste(group_counts$label, collapse = "; "),
    plot_file = output_file
  )

  message("Saved ", output_file)
}

summary_tbl <- bind_rows(summary_rows) |>
  arrange(logrank_p)

write_csv(summary_tbl, summary_file)
message("Saved ", summary_file)
