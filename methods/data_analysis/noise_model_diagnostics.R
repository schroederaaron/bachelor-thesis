#!/usr/bin/env Rscript
# noise_model_diagnostics.R
# Diagnostic plots for evaluating adaptive kNN noise model performance

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(data.table)
library(viridis)

# ==================== CONFIGURATION ====================

# Source config to get FAMILY_OUTPUT_DIR
source("config.R")

# Set base directory consistent with outlier significance analysis
BASE_DIR <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")
RESULTS_FILE <- file.path(BASE_DIR, "all_genes_pvalues.rds")
OUTPUT_DIR <- file.path(BASE_DIR, "diagnostics")

# Color schemes
NORM_COLORS <- c(
  "raw" = "#E41A1C",
  "log" = "#377EB8", 
  "std_log" = "#4DAF4A",
  "full" = "#984EA3"
)

NORM_DISPLAY <- c(
  "raw" = "mean(TPM)",
  "log" = "log(mean(TPM))",
  "std_log" = "log(mean(gene-wise scaled TPM))",
  "full" = "log(mean(quantile(gene-wise scaled TPM)))"
)

COMPARISON_SHAPES <- c(
  "own_healthy" = 16,
  "family_mean" = 17,
  "ortholog_mean" = 15
)

COMPARISON_DISPLAY <- c(
  "own_healthy" = "Gene vs own healthy",
  "family_mean" = "Gene vs family mean",
  "ortholog_mean" = "Gene vs ortholog mean"
)

# ==================== HELPER FUNCTIONS ====================

#' Load and prepare data
load_data <- function(file_path) {
  cat("Loading data from:", file_path, "\n")
  
  if (!file.exists(file_path)) {
    stop("Results file not found: ", file_path)
  }
  
  data <- readRDS(file_path)
  
  # Convert to data.table for faster operations if not already
  if (!is.data.table(data)) {
    data <- as.data.table(data)
  }
  
  # Ensure factors are properly set
  data$norm_method <- factor(data$norm_method, 
                              levels = c("raw", "log", "std_log", "full"))
  data$comparison <- factor(data$comparison, 
                            levels = c("own_healthy", "family_mean", "ortholog_mean"))
  
  # Add display names for plotting
  data$norm_display <- NORM_DISPLAY[as.character(data$norm_method)]
  data$comparison_display <- COMPARISON_DISPLAY[as.character(data$comparison)]
  
  cat(sprintf("Loaded %d rows\n", nrow(data)))
  cat("Normalizations:", paste(levels(data$norm_method), collapse = ", "), "\n")
  cat("Comparisons:", paste(levels(data$comparison), collapse = ", "), "\n")
  cat("Cancer types:", paste(unique(data$cancer_type), collapse = ", "), "\n\n")
  
  return(data)
}

#' Create histogram of neighborhood sizes
plot_neighborhood_histograms <- function(data, output_dir) {
  cat("Creating neighborhood size histograms...\n")
  
  # Cancer neighborhood sizes
  p_cancer <- ggplot(data, aes(x = neighborhood_size_cancer, fill = norm_method)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 50, color = "black", size = 0.1) +
    facet_grid(comparison_display ~ norm_display, scales = "free_y") +
    scale_fill_manual(values = NORM_COLORS) +
    labs(
      title = "Cancer Neighborhood Size Distribution",
      subtitle = "By normalization method and comparison type",
      x = "Neighborhood Size (Cancer)",
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = c(32, 1024), linetype = "dashed", color = "red", alpha = 0.5)
  
  ggsave(file.path(output_dir, "histogram_cancer_neighborhoods.png"), 
         p_cancer, width = 16, height = 10, dpi = 300)
  
  # Healthy neighborhood sizes
  p_healthy <- ggplot(data, aes(x = neighborhood_size_healthy, fill = norm_method)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 50, color = "black", size = 0.1) +
    facet_grid(comparison_display ~ norm_display, scales = "free_y") +
    scale_fill_manual(values = NORM_COLORS) +
    labs(
      title = "Healthy Neighborhood Size Distribution",
      subtitle = "By normalization method and comparison type",
      x = "Neighborhood Size (Healthy)",
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = c(32, 1024), linetype = "dashed", color = "red", alpha = 0.5)
  
  ggsave(file.path(output_dir, "histogram_healthy_neighborhoods.png"), 
         p_healthy, width = 16, height = 10, dpi = 300)
  
  # Combined density plot
  p_density <- ggplot(data, aes(x = neighborhood_size_cancer, color = norm_method, linetype = comparison_display)) +
    geom_density(size = 1.2, alpha = 0.7) +
    scale_color_manual(values = NORM_COLORS) +
    labs(
      title = "Cancer Neighborhood Size Density",
      subtitle = "By normalization and comparison type",
      x = "Neighborhood Size (Cancer)",
      y = "Density",
      color = "Normalization",
      linetype = "Comparison"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = c(32, 1024), linetype = "dashed", color = "gray50", alpha = 0.5)
  
  ggsave(file.path(output_dir, "density_cancer_neighborhoods.png"), 
         p_density, width = 14, height = 8, dpi = 300)
  
  cat("   Saved histogram plots\n")
}

#' Create correlation analysis between cancer and healthy neighborhoods
plot_neighborhood_correlations <- function(data, output_dir) {
  cat("Creating neighborhood correlation analysis...\n")
  
  # Sample data to max 1,000,000 points for plotting
  set.seed(42)
  if (nrow(data) > 1000000) {
    data_plot <- data[sample(1:nrow(data), 1000000), ]
    cat(sprintf("  Sampled %d points from %d total for plotting\n", nrow(data_plot), nrow(data)))
  } else {
    data_plot <- data
  }
  
  # Calculate correlations per group (use full data for accurate correlations)
  cor_summary <- data %>%
    group_by(norm_method, norm_display, comparison, comparison_display) %>%
    summarise(
      pearson_cor = cor(neighborhood_size_cancer, neighborhood_size_healthy, 
                        method = "pearson", use = "complete.obs"),
      spearman_cor = cor(neighborhood_size_cancer, neighborhood_size_healthy, 
                         method = "spearman", use = "complete.obs"),
      n = n(),
      .groups = "drop"
    )
  
  # Save correlation table
  write.csv(cor_summary, file.path(output_dir, "neighborhood_correlations.csv"), 
            row.names = FALSE)
  
  # Print summary
  cat("\n  Correlation Summary:\n")
  print(cor_summary %>% select(norm_method, comparison, spearman_cor, pearson_cor, n))
  
  # Heatmap of correlations
  p_heatmap <- ggplot(cor_summary, aes(x = norm_display, y = comparison_display, fill = spearman_cor)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", spearman_cor)), color = "white", size = 4, fontface = "bold") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    labs(
      title = "Spearman Correlation: Cancer vs Healthy Neighborhood Sizes",
      x = "Normalization Method",
      y = "Comparison Type",
      fill = "Correlation"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(output_dir, "neighborhood_correlation_heatmap.png"), 
         p_heatmap, width = 12, height = 7, dpi = 300)
  
  # Scatter plot with facets - using sampled data
  p_scatter <- ggplot(data_plot, 
                      aes(x = neighborhood_size_cancer, y = neighborhood_size_healthy)) +
    geom_point(alpha = 0.15, size = 0.4, aes(color = norm_method)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5, alpha = 0.2) +
    facet_grid(comparison_display ~ norm_display) +
    scale_color_manual(values = NORM_COLORS) +
    labs(
      title = "Cancer vs Healthy Neighborhood Sizes",
      subtitle = sprintf("Sampled %d points from %d total", nrow(data_plot), nrow(data)),
      x = "Cancer Neighborhood Size",
      y = "Healthy Neighborhood Size"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "neighborhood_scatter_facets.png"), 
         p_scatter, width = 16, height = 10, dpi = 300)
  
  cat("   Saved correlation plots\n")
}

#' Create p-value vs neighborhood size scatter plots
plot_pvalue_vs_neighborhood <- function(data, output_dir) {
  cat("Creating p-value vs neighborhood size plots...\n")
  
  # Sample data to max 1,000,000 points for plotting
  set.seed(42)
  if (nrow(data) > 1000000) {
    data_plot <- data[sample(1:nrow(data), 1000000), ]
    cat(sprintf("  Sampled %d points from %d total for plotting\n", nrow(data_plot), nrow(data)))
  } else {
    data_plot <- data
  }
  
  # Cancer neighborhood vs p-value 
  p_cancer <- ggplot(data_plot, 
                     aes(x = neighborhood_size_cancer, y = noise_p_value_adj)) +
    geom_point(alpha = 0.15, size = 0.4, aes(color = norm_method)) +
    geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8,
                method.args = list(span = 0.5), n = 100) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_grid(comparison_display ~ norm_display) +
    scale_color_manual(values = NORM_COLORS) +
    labs(
      title = "Adjusted Noise P-value vs Cancer Neighborhood Size",
      subtitle = sprintf("Sampled %d points from %d total | Red line: p = 0.01", 
                         nrow(data_plot), nrow(data)),
      x = "Cancer Neighborhood Size",
      y = "Adjusted Noise P-value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_vs_cancer_neighborhood.png"), 
         p_cancer, width = 16, height = 10, dpi = 300)
  
  # Healthy neighborhood vs p-value
  p_healthy <- ggplot(data_plot, 
                      aes(x = neighborhood_size_healthy, y = noise_p_value_adj)) +
    geom_point(alpha = 0.15, size = 0.4, aes(color = norm_method)) +
    geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 0.8,
                method.args = list(span = 0.5), n = 100) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_grid(comparison_display ~ norm_display) +
    scale_color_manual(values = NORM_COLORS) +
    labs(
      title = "Adjusted Noise P-value vs Healthy Neighborhood Size",
      subtitle = sprintf("Sampled %d points from %d total | Red line: p = 0.01", 
                         nrow(data_plot), nrow(data)),
      x = "Healthy Neighborhood Size",
      y = "Adjusted Noise P-value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_vs_healthy_neighborhood.png"), 
         p_healthy, width = 16, height = 10, dpi = 300)
  
  # Combined plot with both neighborhood types
  data_long <- data_plot %>%
    pivot_longer(
      cols = c(neighborhood_size_cancer, neighborhood_size_healthy),
      names_to = "neighborhood_type",
      values_to = "neighborhood_size"
    ) %>%
    mutate(
      neighborhood_type = case_when(
        neighborhood_type == "neighborhood_size_cancer" ~ "Cancer",
        neighborhood_type == "neighborhood_size_healthy" ~ "Healthy"
      )
    )
  
  p_combined <- ggplot(data_long, 
                       aes(x = neighborhood_size, y = noise_p_value_adj, 
                           color = norm_method, shape = neighborhood_type)) +
    geom_point(alpha = 0.15, size = 0.4) +
    geom_smooth(aes(linetype = neighborhood_type), method = "loess", 
                se = FALSE, linewidth = 1, color = "black",
                method.args = list(span = 0.5), n = 100) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red", alpha = 0.7) +
    facet_grid(comparison_display ~ norm_display) +
    scale_color_manual(values = NORM_COLORS) +
    scale_shape_manual(values = c("Cancer" = 16, "Healthy" = 17)) +
    labs(
      title = "Adjusted Noise P-value vs Neighborhood Size",
      subtitle = sprintf("Sampled %d points from %d total | Cancer (circles) vs Healthy (triangles)", 
                         nrow(data_plot), nrow(data)),
      x = "Neighborhood Size",
      y = "Adjusted Noise P-value"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_vs_neighborhood_combined.png"), 
         p_combined, width = 16, height = 12, dpi = 300)
  
  cat("   Saved p-value scatter plots\n")
}

#' Create p-value distribution plots
plot_pvalue_distributions <- function(data, output_dir) {
  cat("Creating p-value distribution plots...\n")
  
  # Histogram of p-values
  p_hist <- ggplot(data, aes(x = noise_p_value_adj, fill = norm_method)) +
    geom_histogram(position = "identity", alpha = 0.6, bins = 50, 
                   color = "black", size = 0.1) +
    facet_grid(comparison_display ~ norm_display) +
    scale_fill_manual(values = NORM_COLORS) +
    geom_vline(xintercept = 0.01, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(
      title = "Adjusted Noise P-value Distribution",
      subtitle = "Red line: p = 0.01 significance threshold",
      x = "Adjusted Noise P-value",
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_histograms.png"), 
         p_hist, width = 16, height = 10, dpi = 300)
  
  # Density plot
  p_density <- ggplot(data, aes(x = noise_p_value_adj, color = norm_method, linetype = comparison_display)) +
    geom_density(size = 1.2, alpha = 0.7) +
    scale_color_manual(values = NORM_COLORS) +
    geom_vline(xintercept = 0.01, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(
      title = "Adjusted Noise P-value Density",
      subtitle = "By normalization and comparison type",
      x = "Adjusted Noise P-value",
      y = "Density",
      color = "Normalization",
      linetype = "Comparison"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_density.png"), 
         p_density, width = 14, height = 8, dpi = 300)
  
  # QQ plot of p-values vs uniform
  p_qq <- data %>%
    group_by(norm_method, comparison, norm_display, comparison_display) %>%
    arrange(noise_p_value_adj) %>%
    mutate(
      theoretical = (1:n()) / (n() + 1)
    ) %>%
    ungroup() %>%
    ggplot(aes(x = theoretical, y = noise_p_value_adj, color = norm_method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(alpha = 0.3, size = 0.5) +
    facet_grid(comparison_display ~ norm_display) +
    scale_color_manual(values = NORM_COLORS) +
    labs(
      title = "QQ Plot: Observed vs Uniform P-values",
      subtitle = "Deviation from diagonal indicates non-uniform distribution",
      x = "Theoretical Uniform Quantiles",
      y = "Observed P-values"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_qq_plot.png"), 
         p_qq, width = 16, height = 10, dpi = 300)
  
  cat("   Saved p-value distribution plots\n")
}

#' Create summary statistics table
create_summary_table <- function(data, output_dir) {
  cat("Creating summary statistics table...\n")
  
  summary_stats <- data %>%
    group_by(norm_method, norm_display, comparison, comparison_display) %>%
    summarise(
      n_total = n(),
      n_valid = sum(!is.na(noise_p_value_adj)),
      n_significant_raw = sum(noise_p_value_adj < 0.01, na.rm = TRUE),
      pct_significant_raw = 100 * n_significant_raw / n_valid,
      median_p = median(noise_p_value_adj, na.rm = TRUE),
      q25_p = quantile(noise_p_value_adj, 0.25, na.rm = TRUE),
      q75_p = quantile(noise_p_value_adj, 0.75, na.rm = TRUE),
      min_p = min(noise_p_value_adj, na.rm = TRUE),
      max_p = max(noise_p_value_adj, na.rm = TRUE),
      median_cancer_pool = median(neighborhood_size_cancer, na.rm = TRUE),
      median_healthy_pool = median(neighborhood_size_healthy, na.rm = TRUE),
      mean_cancer_pool = mean(neighborhood_size_cancer, na.rm = TRUE),
      mean_healthy_pool = mean(neighborhood_size_healthy, na.rm = TRUE),
      cor_spearman = cor(noise_p_value_adj, neighborhood_size_cancer, 
                         method = "spearman", use = "complete.obs"),
      cor_pearson = cor(noise_p_value_adj, neighborhood_size_cancer, 
                        method = "pearson", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    arrange(norm_method, comparison)
  
  write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), 
            row.names = FALSE)
  
  cat("   Saved summary statistics\n")
  
  # Print key findings
  cat("\n  Key Findings:\n")
  cat("  =============\n")
  
  own_data <- summary_stats %>% filter(comparison == "own_healthy")
  fam_data <- summary_stats %>% filter(comparison == "family_mean")
  orth_data <- summary_stats %>% filter(comparison == "ortholog_mean")
  
  cat("\n  own_healthy significance rates:\n")
  for (i in 1:nrow(own_data)) {
    cat(sprintf("    %s: %.1f%% (%d / %d) | Q1 p = %.3f | Median p = %.3f\n", 
                own_data$norm_method[i],
                own_data$pct_significant_raw[i],
                own_data$n_significant_raw[i],
                own_data$n_valid[i],
                own_data$q25_p[i],
                own_data$median_p[i]))
  }
  
  cat("\n  family_mean significance rates:\n")
  for (i in 1:nrow(fam_data)) {
    cat(sprintf("    %s: %.1f%% (%d / %d) | Q1 p = %.3f | Median p = %.3f\n", 
                fam_data$norm_method[i],
                fam_data$pct_significant_raw[i],
                fam_data$n_significant_raw[i],
                fam_data$n_valid[i],
                fam_data$q25_p[i],
                fam_data$median_p[i]))
  }
  
  cat("\n  ortholog_mean significance rates:\n")
  for (i in 1:nrow(orth_data)) {
    cat(sprintf("    %s: %.1f%% (%d / %d) | Q1 p = %.3f | Median p = %.3f\n", 
                orth_data$norm_method[i],
                orth_data$pct_significant_raw[i],
                orth_data$n_significant_raw[i],
                orth_data$n_valid[i],
                orth_data$q25_p[i],
                orth_data$median_p[i]))
  }
}

#' Create boxplot comparisons
plot_boxplot_comparisons <- function(data, output_dir) {
  cat("Creating boxplot comparisons...\n")
  
  # P-value boxplots
  p_box <- ggplot(data, aes(x = norm_display, y = noise_p_value_adj, fill = comparison_display)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.1, outlier.size = 0.5) +
    geom_hline(yintercept = 0.01, linetype = "dashed", color = "red", alpha = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Noise P-value Distribution by Normalization and Comparison",
      x = "Normalization Method",
      y = "Adjusted Noise P-value",
      fill = "Comparison"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "pvalue_boxplots.png"), 
         p_box, width = 12, height = 7, dpi = 300)
  
  # Neighborhood size boxplots
  p_pool_box <- data %>%
    pivot_longer(
      cols = c(neighborhood_size_cancer, neighborhood_size_healthy),
      names_to = "pool_type",
      values_to = "pool_size"
    ) %>%
    mutate(pool_type = ifelse(pool_type == "neighborhood_size_cancer", "Cancer", "Healthy")) %>%
    ggplot(aes(x = norm_display, y = pool_size, fill = comparison_display)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.1, outlier.size = 0.5) +
    facet_wrap(~pool_type, scales = "free_y") +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Neighborhood Size Distribution",
      x = "Normalization Method",
      y = "Pool Size",
      fill = "Comparison"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(output_dir, "neighborhood_boxplots.png"), 
         p_pool_box, width = 14, height = 8, dpi = 300)
  
  cat("   Saved boxplot comparisons\n")
}

#' Create significance heatmap
plot_significance_heatmap <- function(data, output_dir) {
  cat("Creating significance heatmap...\n")
  
  # Calculate significance rates by cancer type and stage
  sig_rates <- data %>%
    group_by(cancer_type, stage, norm_method, norm_display, comparison, comparison_display) %>%
    summarise(
      sig_rate = mean(noise_p_value_adj < 0.01, na.rm = TRUE) * 100,
      n_genes = n(),
      .groups = "drop"
    )
  
  p_heatmap <- ggplot(sig_rates, aes(x = stage, y = cancer_type, fill = sig_rate)) +
    geom_tile(color = "white", size = 0.5) +
    facet_grid(comparison_display ~ norm_display) +
    scale_fill_viridis_c(option = "plasma", name = "% Significant") +
    labs(
      title = "Significance Rate by Cancer Type and Stage",
      subtitle = "Percentage of genes with adjusted noise p-value < 0.01",
      x = "Stage",
      y = "Cancer Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 8),
      panel.grid = element_blank(),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir, "significance_heatmap.png"), 
         p_heatmap, width = 18, height = 12, dpi = 300)
  
  cat("   Saved significance heatmap\n")
}

# ==================== MAIN EXECUTION ====================

main <- function() {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("NOISE MODEL DIAGNOSTICS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  cat("Base directory:", BASE_DIR, "\n")
  cat("Results file:", RESULTS_FILE, "\n")
  cat("Output directory:", OUTPUT_DIR, "\n\n")
  
  # Create output directory
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  data <- load_data(RESULTS_FILE)
  
  # Generate plots
  plot_neighborhood_histograms(data, OUTPUT_DIR)
  plot_neighborhood_correlations(data, OUTPUT_DIR)
  plot_pvalue_vs_neighborhood(data, OUTPUT_DIR)
  plot_pvalue_distributions(data, OUTPUT_DIR)
  plot_boxplot_comparisons(data, OUTPUT_DIR)
  plot_significance_heatmap(data, OUTPUT_DIR)
  
  # Create summary table
  create_summary_table(data, OUTPUT_DIR)
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("DIAGNOSTICS COMPLETE\n")
  cat("All outputs saved to:", OUTPUT_DIR, "\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
}

# Run if script is executed directly
if (sys.nframe() == 0) {
  main()
}