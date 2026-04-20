#!/usr/bin/env Rscript
# plot_pvalue_distributions.R
# Create ECDF plots for p-value distributions
# - Noise p-values (uncorrected vs corrected) per cancer, comparison, and normalization
# - Overlay ECDFs across normalization methods
# - limma adj.P.Val distribution
# - edgeR FDR distribution

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(gridExtra)

source("config.R")
source("utils.R")

# ==================== CONFIGURATION ====================

CANCER_TYPES <- c(
  "breast cancer" = "TCGA-BRCA",
  "non small cell lung cancer" = "TCGA-LUAD",
  "small cell lung cancer" = "TCGA-LUSC",
  "bladder cancer" = "TCGA-BLCA",
  "colon cancer" = "TCGA-COAD",
  "kidney cancer" = "TCGA-KIRC",
  "stomach cancer" = "TCGA-STAD",
  "thyroid gland cancer" = "TCGA-THCA"
)

NORM_METHODS <- c("raw", "log", "std_log", "full")
NORM_DISPLAY <- c(
  "raw" = "raw",
  "log" = "log",
  "std_log" = "std_log",
  "full" = "full"
)

RAW_P_CUTOFF <- 0.01
ADJ_P_CUTOFF <- 0.01
ZOOM_XMAX <- 0.05
OVERLAY_ALPHA <- 0.75

# ==================== DATA LOADING ====================

load_tensoromics_data <- function(project_id) {
  sig_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")
  data_file <- file.path(sig_dir, "all_genes_pvalues.rds")

  if (!file.exists(data_file)) {
    message(sprintf("all_genes_pvalues file not found: %s", data_file))
    return(NULL)
  }

  all_data <- readRDS(data_file)

  required_cols <- c("cancer_id", "comparison", "norm_method", "noise_p_value", "noise_p_value_adj", "gene_id")
  missing_cols <- setdiff(required_cols, names(all_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing columns in all_genes_pvalues.rds: %s", paste(missing_cols, collapse = ", ")))
  }

  project_data <- all_data %>%
    filter(cancer_id == project_id) %>%
    mutate(
      raw_outlier = !is.na(noise_p_value) & noise_p_value < RAW_P_CUTOFF,
      corrected_outlier = !is.na(noise_p_value_adj) & noise_p_value_adj < ADJ_P_CUTOFF
    )

  if (nrow(project_data) == 0) {
    message(sprintf("No data found for project %s", project_id))
    return(NULL)
  }

  return(project_data)
}

load_de_data <- function(project_id, stage) {
  edger_file <- file.path(DE_OUTPUT_DIR, project_id, stage, "edgeR_results.csv")
  edger_data <- NULL
  if (file.exists(edger_file)) {
    edger_data <- read.csv(edger_file)
  }

  limma_file <- file.path(DE_OUTPUT_DIR, project_id, stage, "voom_results.csv")
  limma_data <- NULL
  if (file.exists(limma_file)) {
    limma_data <- read.csv(limma_file)
  }

  return(list(
    edger = edger_data,
    limma = limma_data
  ))
}

# ==================== HELPERS ====================

format_count_label <- function(data) {
  n_total <- nrow(data)
  n_genes <- dplyr::n_distinct(data$gene_id)
  n_raw <- sum(data$raw_outlier, na.rm = TRUE)
  n_adj <- sum(data$corrected_outlier, na.rm = TRUE)
  n_raw_genes <- dplyr::n_distinct(data$gene_id[data$raw_outlier])
  n_adj_genes <- dplyr::n_distinct(data$gene_id[data$corrected_outlier])

  paste0(
    "rows = ", format(n_total, big.mark = ","), "\n",
    "genes = ", format(n_genes, big.mark = ","), "\n",
    "raw p < ", RAW_P_CUTOFF, ": ", format(n_raw, big.mark = ","),
    " rows / ", format(n_raw_genes, big.mark = ","), " genes\n",
    "adj p < ", ADJ_P_CUTOFF, ": ", format(n_adj, big.mark = ","),
    " rows / ", format(n_adj_genes, big.mark = ","), " genes"
  )
}

make_ecdf_panel <- function(data, x_col, threshold, color, title, x_label, label_text) {
  plot_data <- data %>%
    filter(!is.na(.data[[x_col]]))

  if (nrow(plot_data) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No data", size = 5) +
        theme_void() +
        labs(title = title)
    )
  }

  ggplot(plot_data, aes(x = .data[[x_col]])) +
    stat_ecdf(geom = "step", linewidth = 1, color = color) +
    geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 0.8) +
    coord_cartesian(xlim = c(0, ZOOM_XMAX), ylim = c(0, 1)) +
    annotate(
      "label",
      x = ZOOM_XMAX * 0.98,
      y = 0.98,
      label = label_text,
      hjust = 1,
      vjust = 1,
      size = 3.2,
      label.size = 0.2,
      fill = "white",
      alpha = 0.9
    ) +
    labs(
      title = title,
      x = x_label,
      y = "Cumulative proportion"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )
}

# ==================== PLOT FUNCTIONS ====================

plot_noise_ecdf_individual <- function(data, norm_method_local, comp, project_id) {
  norm_data <- data %>%
    filter(norm_method == !!norm_method_local, comparison == !!comp)

  if (nrow(norm_data) == 0) {
    return(NULL)
  }

  label_text <- format_count_label(norm_data)
  subtitle_text <- sprintf("%s - %s - %s", project_id, NORM_DISPLAY[norm_method_local], comp)

  p1 <- make_ecdf_panel(
    data = norm_data,
    x_col = "noise_p_value",
    threshold = RAW_P_CUTOFF,
    color = "steelblue",
    title = "Noise p-values (uncorrected)",
    x_label = "p-value",
    label_text = label_text
  ) + labs(subtitle = subtitle_text)

  p2 <- make_ecdf_panel(
    data = norm_data,
    x_col = "noise_p_value_adj",
    threshold = ADJ_P_CUTOFF,
    color = "darkgreen",
    title = "Noise p-values (BH-corrected)",
    x_label = "adjusted p-value",
    label_text = label_text
  ) + labs(subtitle = subtitle_text)

  p1 + p2 + plot_layout(ncol = 2)
}

plot_noise_ecdf_overlay <- function(data, comp, project_id,
                                    method_order = NORM_METHODS,
                                    alpha_value = OVERLAY_ALPHA) {
  plot_data <- data %>%
    filter(comparison == !!comp, norm_method %in% method_order) %>%
    mutate(norm_method = factor(norm_method, levels = method_order))

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  counts_df <- plot_data %>%
    group_by(norm_method) %>%
    summarise(
      n_genes = n_distinct(gene_id),
      n_raw_genes = n_distinct(gene_id[raw_outlier]),
      n_adj_genes = n_distinct(gene_id[corrected_outlier]),
      .groups = "drop"
    ) %>%
    mutate(
      label = paste0(
        as.character(norm_method), ": ",
        "genes=", n_genes,
        ", uncorrected<", RAW_P_CUTOFF, "=", n_raw_genes,
        ", adjusted<", ADJ_P_CUTOFF, "=", n_adj_genes
      )
    )

  counts_text <- paste(counts_df$label, collapse = "\n")

  p_raw <- ggplot(
    plot_data %>% filter(!is.na(noise_p_value)),
    aes(x = noise_p_value, color = norm_method, group = norm_method)
  ) +
    stat_ecdf(geom = "step", linewidth = 1, alpha = alpha_value) +
    geom_vline(xintercept = RAW_P_CUTOFF, color = "red", linetype = "dashed", linewidth = 0.8) +
    coord_cartesian(xlim = c(0, ZOOM_XMAX), ylim = c(0, 1)) +
    labs(
      title = "Noise p-values (uncorrected)",
      subtitle = sprintf("%s - %s", project_id, comp),
      x = "p-value",
      y = "Cumulative proportion",
      color = "Normalization"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      legend.position = "bottom"
    )

  p_adj <- ggplot(
    plot_data %>% filter(!is.na(noise_p_value_adj)),
    aes(x = noise_p_value_adj, color = norm_method, group = norm_method)
  ) +
    stat_ecdf(geom = "step", linewidth = 1, alpha = alpha_value) +
    geom_vline(xintercept = ADJ_P_CUTOFF, color = "red", linetype = "dashed", linewidth = 0.8) +
    coord_cartesian(xlim = c(0, ZOOM_XMAX), ylim = c(0, 1)) +
    annotate(
      "label",
      x = ZOOM_XMAX * 0.98,
      y = 0.98,
      label = counts_text,
      hjust = 1,
      vjust = 1,
      size = 3.0,
      label.size = 0.2,
      fill = "white",
      alpha = 0.9
    ) +
    labs(
      title = "Noise p-values (BH-corrected)",
      subtitle = sprintf("%s - %s", project_id, comp),
      x = "adjusted p-value",
      y = "Cumulative proportion",
      color = "Normalization"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      plot.subtitle = element_text(size = 8),
      legend.position = "bottom"
    )

  p_raw + p_adj + plot_layout(ncol = 2)
}

plot_edger_histogram <- function(edger_data, project_id, stage) {
  if (is.null(edger_data)) {
    return(NULL)
  }

  p <- ggplot(edger_data, aes(x = FDR)) +
    stat_ecdf(geom = "step", linewidth = 1, color = "darkred") +
    geom_vline(xintercept = 0.05, color = "blue", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = 0.01, color = "purple", linetype = "dashed", linewidth = 1) +
    coord_cartesian(xlim = c(0, 0.05), ylim = c(0, 1)) +
    labs(
      title = sprintf("edgeR FDR ECDF - %s - %s", project_id, stage),
      x = "FDR",
      y = "Cumulative proportion"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}

plot_limma_histogram <- function(limma_data, project_id, stage) {
  if (is.null(limma_data)) {
    return(NULL)
  }

  p <- ggplot(limma_data, aes(x = adj.P.Val)) +
    stat_ecdf(geom = "step", linewidth = 1, color = "darkblue") +
    geom_vline(xintercept = 0.05, color = "blue", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = 0.01, color = "purple", linetype = "dashed", linewidth = 1) +
    coord_cartesian(xlim = c(0, 0.05), ylim = c(0, 1)) +
    labs(
      title = sprintf("limma adj.P.Val ECDF - %s - %s", project_id, stage),
      x = "adjusted p-value",
      y = "Cumulative proportion"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))

  return(p)
}

plot_de_histograms <- function(edger_data, limma_data, project_id, stage) {
  p_edger <- plot_edger_histogram(edger_data, project_id, stage)
  p_limma <- plot_limma_histogram(limma_data, project_id, stage)

  if (is.null(p_edger) && is.null(p_limma)) {
    return(NULL)
  } else if (!is.null(p_edger) && !is.null(p_limma)) {
    return(p_edger + p_limma + plot_layout(ncol = 2))
  } else if (!is.null(p_edger)) {
    return(p_edger)
  } else {
    return(p_limma)
  }
}

# ==================== MAIN LOOP ====================

plot_all_distributions <- function(output_dir) {

  plots_dir <- file.path(output_dir, "ecdf_plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  all_noise_plots <- list()
  all_overlay_plots <- list()
  all_de_plots <- list()

  for (cancer_name in names(CANCER_TYPES)) {
    project_id <- CANCER_TYPES[cancer_name]

    to_data <- load_tensoromics_data(project_id)
    if (!is.null(to_data)) {
      comparisons <- unique(to_data$comparison)

      for (comp in comparisons) {
        p_overlay <- plot_noise_ecdf_overlay(to_data, comp, project_id)

        if (!is.null(p_overlay)) {
          filename <- sprintf("noise_ecdf_overlay_%s_%s.png", project_id, comp)
          ggsave(file.path(plots_dir, filename), p_overlay, width = 13, height = 5.5, dpi = 300)
          all_overlay_plots[[paste(project_id, comp, sep = "_")]] <- p_overlay
        }

        for (norm in NORM_METHODS) {
          p <- plot_noise_ecdf_individual(to_data, norm, comp, project_id)

          if (!is.null(p)) {
            filename <- sprintf("noise_ecdf_%s_%s_%s.png", project_id, norm, comp)
            ggsave(file.path(plots_dir, filename), p, width = 12, height = 5, dpi = 300)
            all_noise_plots[[paste(project_id, norm, comp, sep = "_")]] <- p
          }
        }
      }
    }

    for (stage in STAGES) {
      de_data <- load_de_data(project_id, stage)

      if (!is.null(de_data$edger) || !is.null(de_data$limma)) {
        p <- plot_de_histograms(de_data$edger, de_data$limma, project_id, stage)

        if (!is.null(p)) {
          filename <- sprintf("de_ecdf_%s_%s.png", project_id, gsub(" ", "_", stage))
          ggsave(file.path(plots_dir, filename), p, width = 12, height = 5, dpi = 300)
          all_de_plots[[paste(project_id, stage, sep = "_")]] <- p
        }
      }
    }
  }


  return(list(
    noise_plots = all_noise_plots,
    overlay_plots = all_overlay_plots,
    de_plots = all_de_plots
  ))
}

# ==================== EXECUTION ====================

if (sys.nframe() == 0) {
  output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "pvalue_distributions")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  results <- plot_all_distributions(output_dir)

}