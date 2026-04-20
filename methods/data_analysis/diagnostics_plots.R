#!/usr/bin/env Rscript
# plot_diagnostic_plots.R
# Creates diagnostic plots for noise model calibration
# Uses SAVED data from outlier_significance_analysis.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(gridExtra)
library(RColorBrewer)

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
  "raw" = "mean(TPM)",
  "log" = "log(mean(TPM))",
  "std_log" = "log(gene-wise scaled TPM)",
  "full" = "log(mean(quantile(gene-wise scaled TPM)))"
)

# ==================== DATA LOADING ====================

load_healthy_replicates <- function(project_id) {
  
  # Load healthy data in RAW format (TPM)
  healthy_data <- load_stage_data(
    project_id, 
    STAGES[1], 
    "healthy", 
    use_constant_healthy = TRUE, 
    norm_method = "raw",
    apply_mean = FALSE
  )
  
  if (is.null(healthy_data) || !is.matrix(healthy_data$expression_vectors)) {
    cat(sprintf("  Healthy RAW Data not found for %s\n", project_id))
    return(NULL)
  }
  
  return(healthy_data$expression_vectors)  # samples × genes
}

# ==================== PLOT 5: Voom-like Plots (Healthy vs Cancerous) ====================

load_cancerous_replicates <- function(project_id) {
  
  # Load cancerous data in RAW format (TPM)
  # Note: Adjust based on how your cancer data is organized
  
  # Option 1: All tumor samples regardless of stage
  cancer_data <- load_stage_data(
    project_id, 
    STAGES[1],  # or NULL for all stages
    "cancer",   # Adjust this parameter in your load_stage_data function if needed
    use_constant_healthy = TRUE, 
    norm_method = "raw",
    apply_mean = FALSE
  )
  
  if (is.null(cancer_data) || !is.matrix(cancer_data$expression_vectors)) {
    cat(sprintf("  Cancerous RAW Data not found for %s\n", project_id))
    return(NULL)
  }
  
  return(cancer_data$expression_vectors)  # samples × genes
}

plot_voom_like_separate <- function(healthy_replicates, cancerous_replicates, project_id, output_dir) {
  
  cat(sprintf("\n  Creating Voom-like Plot: SD of log2(TPM+1) vs Mean of log2(TPM+1)\n"))
  cat(sprintf("    Healthy: %d samples, %d genes\n", nrow(healthy_replicates), ncol(healthy_replicates)))
  cat(sprintf("    Cancerous: %d samples, %d genes\n", nrow(cancerous_replicates), ncol(cancerous_replicates)))
  
  # ===== HELPER FUNCTION: Compute for one group =====
  compute_voom_data <- function(expression_matrix, group_name) {
    cat(sprintf("    Processing %s: %d samples, %d genes\n", 
                group_name, nrow(expression_matrix), ncol(expression_matrix)))
    
    # Apply log2(TPM + 1) transformation
    log_expr <- log2(expression_matrix + 1)
    
    # Compute gene means on log-transformed data
    gene_means <- colMeans(log_expr, na.rm = TRUE)
    
    # Compute gene standard deviations on log-transformed data
    gene_sd <- apply(log_expr, 2, sd, na.rm = TRUE)
    
    # Remove invalid values (NA, zero, negative)
    valid_idx <- which(!is.na(gene_means) & !is.na(gene_sd) & 
                        gene_means > 0 & gene_sd > 0)
    
    means_valid <- gene_means[valid_idx]
    sd_valid <- gene_sd[valid_idx]
    
    data.frame(
      group = group_name,
      log2_mean = means_valid,
      log2_sd = sd_valid,
      stringsAsFactors = FALSE
    )
  }
  
  # Compute data for both groups
  healthy_data <- compute_voom_data(healthy_replicates, "Healthy")
  cancerous_data <- compute_voom_data(cancerous_replicates, "Cancerous")
  
  # Check if data exists
  if (nrow(healthy_data) == 0 || nrow(cancerous_data) == 0) {
    cat("  Error: No valid data for Voom-like plot\n")
    return(NULL)
  }
  
  # Combine for consistent axis scaling
  all_data <- rbind(healthy_data, cancerous_data)
  
  # Determine common axis limits (1st-99th percentile)
  x_min <- quantile(all_data$log2_mean, 0.01, na.rm = TRUE)
  x_max <- quantile(all_data$log2_mean, 0.99, na.rm = TRUE)
  y_min <- quantile(all_data$log2_sd, 0.01, na.rm = TRUE)
  y_max <- quantile(all_data$log2_sd, 0.99, na.rm = TRUE)
  
  cat(sprintf("\n  Common axis scaling:\n"))
  cat(sprintf("    x-axis: [%.2f, %.2f] (log2)\n", x_min, x_max))
  cat(sprintf("    y-axis: [%.2f, %.2f] (log2)\n", y_min, y_max))
  
  # ===== LOESS function for one group =====
  fit_loess <- function(data) {
    if (nrow(data) < 100) return(NULL)
    
    # Remove Inf and NA
    data_clean <- data[is.finite(data$log2_mean) & is.finite(data$log2_sd), ]
    if (nrow(data_clean) < 100) return(NULL)
    
    loess_fit <- loess(log2_sd ~ log2_mean, data = data_clean, 
                       span = 0.3, control = loess.control(surface = "direct"))
    
    x_range <- seq(max(x_min, min(data_clean$log2_mean, na.rm = TRUE)), 
                   min(x_max, max(data_clean$log2_mean, na.rm = TRUE)), 
                   length.out = 200)
    
    loess_pred <- predict(loess_fit, newdata = data.frame(log2_mean = x_range), se = TRUE)
    
    data.frame(
      log2_mean = x_range,
      fitted = loess_pred$fit,
      lower = loess_pred$fit - 1.96 * loess_pred$se,
      upper = loess_pred$fit + 1.96 * loess_pred$se
    )
  }
  
  # Compute LOESS for both groups
  healthy_loess <- fit_loess(healthy_data)
  cancerous_loess <- fit_loess(cancerous_data)
  
  # ===== PLOT 1: HEALTHY =====
  # Downsampling for Healthy
  healthy_plot_data <- healthy_data
  if (nrow(healthy_plot_data) > 500000) {
    set.seed(42)
    healthy_plot_data <- healthy_plot_data[sample(1:nrow(healthy_plot_data), 500000), ]
  }
  
  p_healthy <- ggplot(healthy_plot_data, aes(x = log2_mean, y = log2_sd)) +
    geom_point(alpha = 0.15, size = 0.5, color = "darkgreen")
  
  # Add LOESS only if available
  if (!is.null(healthy_loess)) {
    p_healthy <- p_healthy +
      geom_ribbon(data = healthy_loess, 
                  aes(x = log2_mean, ymin = lower, ymax = upper),
                  fill = "red", alpha = 0.25, inherit.aes = FALSE) +
      geom_line(data = healthy_loess, 
                aes(x = log2_mean, y = fitted),
                color = "red", linewidth = 1.2, inherit.aes = FALSE)
  }
  
  p_healthy <- p_healthy +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    labs(
      title = sprintf("Healthy Tissue - %s", project_id),
      x = expression("Mean of log"[2] * "(TPM + 1)"),
      y = expression("SD of log"[2] * "(TPM + 1)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
      panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5)
    )
  
  # ===== PLOT 2: CANCEROUS =====
  # Downsampling for Cancerous
  cancerous_plot_data <- cancerous_data
  if (nrow(cancerous_plot_data) > 500000) {
    set.seed(42)
    cancerous_plot_data <- cancerous_plot_data[sample(1:nrow(cancerous_plot_data), 500000), ]
  }
  
  p_cancerous <- ggplot(cancerous_plot_data, aes(x = log2_mean, y = log2_sd)) +
    geom_point(alpha = 0.15, size = 0.5, color = "darkred")
  
  # Add LOESS only if available
  if (!is.null(cancerous_loess)) {
    p_cancerous <- p_cancerous +
      geom_ribbon(data = cancerous_loess, 
                  aes(x = log2_mean, ymin = lower, ymax = upper),
                  fill = "blue", alpha = 0.25, inherit.aes = FALSE) +
      geom_line(data = cancerous_loess, 
                aes(x = log2_mean, y = fitted),
                color = "blue", linewidth = 1.2, inherit.aes = FALSE)
  }
  
  p_cancerous <- p_cancerous +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    labs(
      title = sprintf("Cancerous Tissue - %s", project_id),
      subtitle = "Cancerous Expression in Stage I",
      x = expression("Mean of log"[2] * "(TPM + 1)"),
      y = expression("SD of log"[2] * "(TPM + 1)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40"),
      panel.border = element_rect(fill = NA, color = "gray70", linewidth = 0.5)
    )
  
  # ===== DIRECT COMPARISON PLOT (both groups overlaid) =====
  # For comparison, plot both groups together
  comparison_data <- rbind(
    cbind(healthy_plot_data, type = "Healthy"),
    cbind(cancerous_plot_data, type = "Cancerous")
  )
  
  p_comparison <- ggplot(comparison_data, aes(x = log2_mean, y = log2_sd, color = type)) +
    geom_point(alpha = 0.1, size = 0.3) +
    geom_smooth(method = "loess", span = 0.3, se = TRUE, alpha = 0.2, linewidth = 1.2) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    scale_color_manual(values = c("Healthy" = "darkgreen", "Cancerous" = "darkred"),
                       name = "Tissue Type") +
    labs(
      title = sprintf("Direct Comparison: Healthy vs Cancerous - %s", project_id),
      subtitle = "LOESS curves with 95% confidence interval",
      x = expression("Mean of log"[2] * "(TPM + 1)"),
      y = expression("SD of log"[2] * "(TPM + 1)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "bottom"
    )
  
  # ===== COMBINED PLOT (side by side) =====
  combined <- (p_healthy | p_cancerous) / p_comparison + 
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = sprintf("Expression Trend Analysis: %s", project_id),
      subtitle = expression("Standard Deviation vs. Mean of log"[2] * "(TPM + 1) transformed expression"),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    )
  
  # ===== SAVE =====
  # Ensure output_dir exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Individual plots
  ggsave(file.path(output_dir, sprintf("voom_like_healthy_%s.png", project_id)), 
         p_healthy, width = 6, height = 5, dpi = 300)
  
  ggsave(file.path(output_dir, sprintf("voom_like_cancerous_%s.png", project_id)), 
         p_cancerous, width = 6, height = 5, dpi = 300)
  
  ggsave(file.path(output_dir, sprintf("voom_like_comparison_%s.png", project_id)), 
         p_comparison, width = 7, height = 6, dpi = 300)
  
  ggsave(file.path(output_dir, sprintf("voom_like_combined_%s.png", project_id)), 
         combined, width = 10, height = 10, dpi = 300)
  
  cat(sprintf("\n  Voom-like plots saved to: %s\n", output_dir))
  cat(sprintf("    - voom_like_healthy_%s.png\n", project_id))
  cat(sprintf("    - voom_like_cancerous_%s.png\n", project_id))
  cat(sprintf("    - voom_like_comparison_%s.png\n", project_id))
  cat(sprintf("    - voom_like_combined_%s.png\n", project_id))
  
  # Return all plots
  return(list(
    plot_healthy = p_healthy,
    plot_cancerous = p_cancerous,
    plot_comparison = p_comparison,
    plot_combined = combined,
    data_healthy = healthy_data,
    data_cancerous = cancerous_data
  ))
}

# ==================== PLOT 1: Mean Expression vs Replicate Values ====================

plot_mean_vs_replicates <- function(healthy_replicates, project_id, output_dir) {
  
  cat(sprintf("\n  Creating Plot 1: Mean Expression vs Replicate Values\n"))
  
  n_genes <- ncol(healthy_replicates)
  n_samples <- nrow(healthy_replicates)
  
  # Compute gene means
  gene_means <- colMeans(healthy_replicates, na.rm = TRUE)
  
  # Use all genes
  set.seed(42)
  sample_idx <- 1:n_genes
  
  # Extract data for selected genes (vectorized!)
  sampled_means <- gene_means[sample_idx]
  sampled_replicates <- healthy_replicates[, sample_idx, drop = FALSE]
  
  # Create data frame by reshaping the matrix
  plot_data <- data.frame(
    mean_expr = rep(sampled_means, each = n_samples),
    replicate = as.vector(sampled_replicates)
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$replicate), ]
  
  # Downsampling to maximum 1,000,000 points for better performance
  if (nrow(plot_data) > 1000000) {
    set.seed(42)
    plot_data <- plot_data[sample(1:nrow(plot_data), 1000000), ]
  }
  
  # Remove extreme outliers
  q99 <- quantile(plot_data$replicate, 0.99, na.rm = TRUE)
  plot_data <- plot_data %>% filter(replicate < q99)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = mean_expr, y = replicate)) +
    geom_point(alpha = 0.1, size = 0.5, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(
      title = sprintf("Mean Expression vs Replicate Values - %s", project_id),
      subtitle = sprintf("%d Genes (all), 1M Points (downsampled)", n_genes),
      x = "Mean Expression (TPM)",
      y = "Replicate Value (TPM)"
    ) +
    theme_minimal()
  
  # Log-Log version
  p_log <- ggplot(plot_data, aes(x = log1p(mean_expr), y = log1p(replicate))) +
    geom_point(alpha = 0.1, size = 0.5, color = "steelblue") +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
    labs(
      title = sprintf("(log-log): Mean Expression vs Replicate Values - %s", project_id),
      x = "log(Mean Expression + 1)",
      y = "log(Replicate Value + 1)"
    ) +
    theme_minimal()
  
  combined <- p + p_log + plot_layout(ncol = 2)
  
  ggsave(file.path(output_dir, sprintf("plot1_mean_vs_replicates_%s.png", project_id)), 
         combined, width = 14, height = 6, dpi = 300)
  
  return(combined)
}

# ==================== PLOT 2: Mean Expression vs Signed Residuals ====================

plot_mean_vs_residuals <- function(healthy_replicates, project_id, output_dir) {
  
  cat(sprintf("\n  Creating Plot 2: Mean Expression vs Signed Residuals\n"))
  
  n_genes <- ncol(healthy_replicates)
  n_samples <- nrow(healthy_replicates)
  
  # Compute gene means
  gene_means <- colMeans(healthy_replicates, na.rm = TRUE)
  
  # Compute residuals matrix (vectorized)
  residuals_matrix <- healthy_replicates - matrix(gene_means, nrow = n_samples, ncol = n_genes, byrow = TRUE)
  
  # Select all genes
  set.seed(42)
  sample_idx <- 1:n_genes
  
  # Extract data for selected genes
  sampled_means <- gene_means[sample_idx]
  sampled_residuals <- residuals_matrix[, sample_idx, drop = FALSE]
  
  # Create data frame
  plot_data <- data.frame(
    mean_expr = rep(sampled_means, each = n_samples),
    residual = as.vector(sampled_residuals)
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$residual), ]
  
  # Downsampling
  if (nrow(plot_data) > 1000000) {
    set.seed(42)
    plot_data <- plot_data[sample(1:nrow(plot_data), 1000000), ]
  }
  
  # ===== NEW: Compute percentiles for better axis scaling =====
  q1 <- quantile(plot_data$mean_expr, 0.25, na.rm = TRUE)
  q99 <- quantile(plot_data$mean_expr, 0.99, na.rm = TRUE)
  median_expr <- median(plot_data$mean_expr, na.rm = TRUE)
  
  # Remove extreme outliers for better visualization (99th percentile)
  q99_resid <- quantile(abs(plot_data$residual), 0.99, na.rm = TRUE)
  plot_data <- plot_data %>% filter(abs(residual) < q99_resid)
  
  # ===== MAIN PLOT WITH LOG SCALING =====
  p <- ggplot(plot_data, aes(x = mean_expr, y = residual)) +
    geom_point(alpha = 0.1, size = 0.5, color = "darkgreen") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    scale_x_log10(  # Log10 scaling for better distribution!
      breaks = c(0.1, 1, 10, 100, 1000, 10000),
      labels = c("0.1", "1", "10", "100", "1000", "10000")
    ) +
    annotation_logticks(sides = "b") +
    labs(
      title = sprintf("Mean Expression vs Signed Residuals - %s", project_id),
      subtitle = sprintf("%d Genes, max 1M Points | Log10 x-Axis | 99th Percentile cap", 
                         n_genes),
      x = "Mean Expression (TPM, log10 Scale)",
      y = "Signed Residual"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  
  # ===== ALTERNATIVE: Boxplot version for better overview =====
  # Create bins for x-axis
  plot_data_binned <- plot_data %>%
    mutate(
      expr_bin = cut(log10(mean_expr), 
                     breaks = seq(-2, 5, by = 0.5),
                     labels = paste0("10^", seq(-2, 4.5, by = 0.5)))
    ) %>%
    group_by(expr_bin) %>%
    summarise(
      median_resid = median(residual, na.rm = TRUE),
      q25_resid = quantile(residual, 0.25, na.rm = TRUE),
      q75_resid = quantile(residual, 0.75, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(expr_bin))
  
  p_binned <- ggplot(plot_data_binned, aes(x = expr_bin, y = median_resid)) +
    geom_errorbar(aes(ymin = q25_resid, ymax = q75_resid), width = 0.2, alpha = 0.5) +
    geom_point(size = 2, color = "darkgreen") +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.8) +
    labs(
      title = sprintf("Binned Mean Expression vs Signed Residuals - %s", project_id),
      subtitle = "Medians and Quartiles per Expression Bin",
      x = "Mean Expression (log10 Bins)",
      y = "Signed Residual (Median with IQR)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    )
  
  # Combined version
  combined <- p + p_binned + plot_layout(ncol = 2)
  
  ggsave(file.path(output_dir, sprintf("plot2_mean_vs_residuals_%s.png", project_id)), 
         combined, width = 14, height = 6, dpi = 300)
  
  return(combined)
}

plot_residual_distribution <- function(healthy_replicates, project_id, output_dir) {
  
  gene_means <- colMeans(healthy_replicates, na.rm = TRUE)
  residuals <- as.vector(healthy_replicates - matrix(gene_means, 
                                                     nrow = nrow(healthy_replicates), 
                                                     ncol = ncol(healthy_replicates), 
                                                     byrow = TRUE))
  residuals <- residuals[!is.na(residuals)]
  
  # Remove extreme values for better display
  q99 <- quantile(abs(residuals), 0.99, na.rm = TRUE)
  residuals <- residuals[abs(residuals) < q99]
  
  p <- ggplot(data.frame(residual = residuals), aes(x = residual)) +
    geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = sprintf("Residual Distribution - %s", project_id),
      x = "Residual",
      y = "Frequency"
    ) +
    theme_minimal()
  
  ggsave(file.path(output_dir, sprintf("plot2_residual_distribution_%s.png", project_id)), 
         p, width = 10, height = 6, dpi = 300)
  
  return(p)
}

# ==================== MAIN FUNCTION ====================

create_diagnostic_plots <- function(output_dir) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("DIAGNOSTIC PLOTS FOR EMPIRICAL NOISE CALIBRATION\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  plots_dir <- file.path(output_dir, "diagnostic_plots_adaptive_knn")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  all_plots <- list()
  
  # ===== PROJECT LOOP =====
  for (cancer_name in names(CANCER_TYPES)) {
    project_id <- CANCER_TYPES[cancer_name]
    cat(sprintf("\n>>> Processing: %s (%s)\n", cancer_name, project_id))
    
    # ===== PLOTS 1 & 2: Require Healthy Replicates =====
    healthy_replicates <- load_healthy_replicates(project_id)
    
    if (!is.null(healthy_replicates)) {
      # Plot 1
      p1 <- plot_mean_vs_replicates(healthy_replicates, project_id, plots_dir)
      all_plots[[paste(project_id, "plot1", sep = "_")]] <- p1
      
      # Plot 2
      p2 <- plot_mean_vs_residuals(healthy_replicates, project_id, plots_dir)
      all_plots[[paste(project_id, "plot2", sep = "_")]] <- p2
      
      # ===== NEW: PLOT 5 - Voom-like Plots =====
      # Load Cancerous data
      cancerous_replicates <- load_cancerous_replicates(project_id)
      
      if (!is.null(cancerous_replicates)) {
        # Ensure both matrices have the same genes
        common_genes <- intersect(colnames(healthy_replicates), colnames(cancerous_replicates))
        
        if (length(common_genes) > 1000) {
          healthy_common <- healthy_replicates[, common_genes]
          cancerous_common <- cancerous_replicates[, common_genes]
          
          p5 <- plot_voom_like_separate(healthy_common, cancerous_common, project_id, plots_dir)
          all_plots[[paste(project_id, "plot5", sep = "_")]] <- p5
        } else {
          cat(sprintf("  Too few common genes for Voom-like plot: %d\n", length(common_genes)))
        }
      } else {
        cat(sprintf("  No Cancerous data available for %s\n", project_id))
      }
    }
  }
  
  # ===== Summary =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("DIAGNOSTIC PLOTS COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("Plots saved to: %s\n", plots_dir))
  cat(sprintf("Number of plots created: %d\n", length(all_plots)))
  
  return(invisible(all_plots))
}

# ==================== EXECUTION ====================

if (sys.nframe() == 0) {
  
  # Output directory
  output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR))
  # dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create plots for DEFAULT parameters
  results <- create_diagnostic_plots(output_dir)
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("DIAGNOSTIC PLOTS COMPLETED SUCCESSFULLY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
}
