#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(methods)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
})

project_root <- normalizePath(
  if (basename(getwd()) == "scripts") file.path(getwd(), "..") else getwd(),
  mustWork = TRUE
)
mirna_data_dir <- file.path(project_root, "data", "mirna_data")
mirna_plot_dir <- file.path(project_root, "miRNA_clusters", "mirna_plots")
dir.create(mirna_plot_dir, recursive = TRUE, showWarnings = FALSE)

k <- 5L
B <- 100L
set.seed(20260708)

poster_cluster_palette <- c(
  "1" = "#4F6FA8",
  "2" = "#7F9CCB",
  "3" = "#7A7F87",
  "4" = "#D1837D",
  "5" = "#A83F4D"
)

mirna_lcpm <- readRDS(file.path(mirna_data_dir, "mirna_lcpm_matrix.rds"))
subtypes <- read.csv(
  file.path(mirna_data_dir, "mirna_nmf_subtypes_k5.csv"),
  stringsAsFactors = FALSE
) |>
  filter(sample_id %in% colnames(mirna_lcpm)) |>
  mutate(cluster = as.integer(cluster)) |>
  arrange(sample_id)

sample_ids <- subtypes$sample_id
original_cluster <- setNames(subtypes$cluster, subtypes$sample_id)
original_sets <- split(names(original_cluster), original_cluster)

nmf_input <- mirna_lcpm[, sample_ids, drop = FALSE]
nmf_input <- nmf_input - min(nmf_input)
n_samples <- ncol(nmf_input)

jaccard <- function(a, b) {
  union_n <- length(union(a, b))
  if (union_n == 0L) return(NA_real_)
  length(intersect(a, b)) / union_n
}

bootstrap_one <- function(iteration) {
  boot_idx <- sample(seq_len(n_samples), n_samples, replace = TRUE)
  boot_samples <- sample_ids[boot_idx]
  inbag_samples <- unique(boot_samples)
  boot_input <- nmf_input[, boot_idx, drop = FALSE]
  colnames(boot_input) <- paste0(boot_samples, "__boot", seq_along(boot_idx))

  fit <- NMF::nmf(
    boot_input,
    rank = k,
    method = "brunet",
    nrun = 1,
    .opt = "v0"
  )

  boot_cluster_raw <- as.integer(NMF::predict(fit))
  boot_sample_ids <- names(NMF::predict(fit))
  if (is.null(boot_sample_ids)) {
    boot_sample_ids <- colnames(boot_input)
  }
  boot_df <- tibble(
    boot_sample_id = boot_sample_ids,
    sample_id = sub("__boot[0-9]+$", "", boot_sample_ids),
    boot_cluster = boot_cluster_raw
  ) |>
    distinct(sample_id, boot_cluster)

  boot_sets <- split(boot_df$sample_id, boot_df$boot_cluster)
  tibble(original_cluster = seq_len(k)) |>
    rowwise() |>
    mutate(
      n_original = length(original_sets[[as.character(original_cluster)]]),
      n_inbag_original = length(intersect(
        original_sets[[as.character(original_cluster)]],
        inbag_samples
      )),
      best_bootstrap_cluster = {
        original_inbag <- intersect(
          original_sets[[as.character(original_cluster)]],
          inbag_samples
        )
        scores <- vapply(boot_sets, function(x) jaccard(original_inbag, x), numeric(1))
        as.integer(names(scores)[which.max(scores)])
      },
      best_jaccard = {
        original_inbag <- intersect(
          original_sets[[as.character(original_cluster)]],
          inbag_samples
        )
        scores <- vapply(boot_sets, function(x) jaccard(original_inbag, x), numeric(1))
        unname(max(scores, na.rm = TRUE))
      }
    ) |>
    ungroup() |>
    mutate(iteration = iteration, .before = 1)
}

message("Running bootstrap-Jaccard stability: k=", k, ", B=", B, ", nrun=1")
all_results <- vector("list", B)
for (i in seq_len(B)) {
  elapsed <- system.time({
    all_results[[i]] <- bootstrap_one(i)
  })[["elapsed"]]
  message(
    sprintf(
      "bootstrap %02d/%02d complete (%.1f sec); mean Jaccard = %.3f",
      i,
      B,
      elapsed,
      mean(all_results[[i]]$best_jaccard)
    )
  )
}

jaccard_results <- bind_rows(all_results) |>
  mutate(original_cluster = factor(original_cluster, levels = seq_len(k)))

jaccard_summary <- jaccard_results |>
  group_by(original_cluster) |>
  summarize(
    n_original = first(n_original),
    mean_inbag_original = round(mean(n_inbag_original), 1),
    mean_jaccard = mean(best_jaccard),
    median_jaccard = median(best_jaccard),
    q25_jaccard = quantile(best_jaccard, 0.25),
    q75_jaccard = quantile(best_jaccard, 0.75),
    min_jaccard = min(best_jaccard),
    stable_ge_0.75 = mean(best_jaccard >= 0.75),
    stable_ge_0.85 = mean(best_jaccard >= 0.85),
    .groups = "drop"
  ) |>
  mutate(across(ends_with("jaccard") | starts_with("stable"), ~ round(.x, 3)))

write_csv(
  jaccard_results,
  file.path(mirna_data_dir, "mirna_nmf_k5_bootstrap_jaccard.csv")
)
write_csv(
  jaccard_summary,
  file.path(mirna_data_dir, "mirna_nmf_k5_bootstrap_jaccard_summary.csv")
)

jaccard_plot <- ggplot(jaccard_results, aes(original_cluster, best_jaccard, fill = original_cluster)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "#555555", linewidth = 0.4) +
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "#777777", linewidth = 0.4) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, height = 0, size = 1.4, alpha = 0.55, color = "#222222") +
  scale_fill_manual(values = poster_cluster_palette, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Bootstrap-Jaccard stability of miRNA NMF clusters (k5)",
    subtitle = paste0(B, " bootstrap resamples; rank-5 NMF refit with nrun=1 per resample"),
    x = "Original NMF cluster",
    y = "Best in-bag Jaccard similarity"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title = element_text(color = "#222222")
  )

ggsave(
  file.path(mirna_plot_dir, "poster_mirna_nmf_k5_bootstrap_jaccard.png"),
  jaccard_plot,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

print(jaccard_summary)
