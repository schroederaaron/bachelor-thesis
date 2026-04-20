# plot_family_trajectories.R
# Complete updated plotting script for Family Trajectory Analysis
# - top families
# - Stage-specific outlier heatmaps
# - Velocity multi-line and heatmap plots
# - Meaningful stage/velocity percentages (vary by stage/transition)

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(gridExtra)
library(RColorBrewer)

# ==================== CONFIGURATION ====================

source("config.R")
source("utils.R")

# ==================== LOAD RESULTS ====================

#' Load family results for a specific normalization method
#' @param project_id Project ID
#' @param norm_method Normalization method: "raw", "full", "log", or "std_log"
#' @return List with all family results
load_family_results_norm <- function(project_id, norm_method) {
  
  norm_display <- get_norm_display(norm_method)
  norm_suffix <- get_norm_suffix(norm_method)
  
  # Get normalization-specific directory
  norm_dir <- get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method)
  proj_dir <- file.path(norm_dir, project_id)
  
  if (!dir.exists(proj_dir)) {
    cat(sprintf("  Family directory not found for [%s]: %s\n", norm_display, proj_dir))
    return(NULL)
  }
  
  # Check for required files
  required_files <- c(
    paste0(project_id, "_family_results", norm_suffix, ".rds"),
    paste0(project_id, "_family_shifts", norm_suffix, ".rds"),
    paste0(project_id, "_family_summary", norm_suffix, ".rds")
  )
  
  missing_files <- required_files[!file.exists(file.path(proj_dir, required_files))]
  if (length(missing_files) > 0) {
    cat(sprintf("  Missing family files for %s [%s]: %s\n", 
                project_id, norm_display, paste(missing_files, collapse = ", ")))
    return(NULL)
  }
  
  results <- list(
    project_id = project_id,
    norm_method = norm_method,
    norm_display = norm_display,
    norm_suffix = norm_suffix,
    trajectory = readRDS(file.path(proj_dir, paste0(project_id, "_family_results", norm_suffix, ".rds"))),
    stage_shifts = readRDS(file.path(proj_dir, paste0(project_id, "_family_shifts", norm_suffix, ".rds"))),
    summary = readRDS(file.path(proj_dir, paste0(project_id, "_family_summary", norm_suffix, ".rds")))
  )
  
  # Load optional files
  outliers_file <- file.path(proj_dir, paste0(project_id, "_family_outliers", norm_suffix, ".rds"))
  velocities_file <- file.path(proj_dir, paste0(project_id, "_family_velocities", norm_suffix, ".rds"))
  loess_file <- file.path(proj_dir, paste0(project_id, "_family_loess_stats", norm_suffix, ".rds"))
  
  if (file.exists(outliers_file)) {
    results$outliers <- readRDS(outliers_file)
  }
  
  if (file.exists(velocities_file)) {
    results$velocities <- readRDS(velocities_file)
  }
  
  if (file.exists(loess_file)) {
    results$loess_stats <- readRDS(loess_file)
  }
  
  return(results)
}

# ==================== DISTRIBUTION PLOTS ====================

#' Plot distribution of scaled total distances
plot_family_scaled_total_distance <- function(results, method = "all") {
  
  if (method == "all") {
    distance_col <- "total_distance_scaled_all"
    outlier_col <- "is_total_outlier_all"
    title_suffix <- "All Genes"
  } else if (method == "ortholog") {
    distance_col <- "total_distance_scaled_ortholog"
    outlier_col <- "is_total_outlier_ortholog"
    title_suffix <- "Orthologs Only"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  if (!distance_col %in% names(results$trajectory)) {
    cat(sprintf("  Warning: Scaled distance column not found for %s [%s]\n", 
                method, results$norm_display))
    return(NULL)
  }
  
  plot_data <- results$trajectory
  plot_data <- plot_data[!is.na(plot_data[[distance_col]]), ]
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  percentile_95 <- quantile(plot_data[[distance_col]], 0.95, na.rm = TRUE)
  n_outliers <- sum(plot_data[[outlier_col]], na.rm = TRUE)
  total_families <- nrow(plot_data)
  families_with_orthologs <- if (method == "ortholog") sum(plot_data$has_orthologs, na.rm = TRUE) else total_families
  
  p <- ggplot(plot_data, aes(x = .data[[distance_col]])) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = percentile_95, 
               color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = percentile_95 * 1.1, y = Inf, 
             label = paste("95th percentile =", round(percentile_95, 2)),
             vjust = 2, hjust = 0, size = 3, color = "red") +
    labs(
      title = sprintf("Family Total Trajectory Distances (Scaled) - %s", results$project_id),
      subtitle = sprintf("%s [%s]: %d outlier families (top 5%%)\nTotal families: %d", 
                        title_suffix, results$norm_display, n_outliers, families_with_orthologs),
      x = "Scaled Total Distance (Distance / Family SD)",
      y = "Number of Families"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  return(p)
}

# ==================== STAGE-SPECIFIC HEATMAP ====================

#' Heatmap showing per-stage outlier status for top families
plot_family_stage_outlier_heatmap <- function(results, method = "all", n_top = 20) {
  
  if (method == "all") {
    shift_matrix <- results$stage_shifts$all_genes
    if (is.null(results$outliers) || is.null(results$outliers$stage_outliers$all_genes)) {
      cat(sprintf("  No stage outlier matrix found for all genes [%s]\n", results$norm_display))
      return(NULL)
    }
    outlier_matrix <- results$outliers$stage_outliers$all_genes
    title_suffix <- "All Genes"
    distance_col <- "total_distance_scaled_all"
  } else if (method == "ortholog") {
    shift_matrix <- results$stage_shifts$orthologs_only
    if (is.null(results$outliers) || is.null(results$outliers$stage_outliers$orthologs_only)) {
      cat(sprintf("  No stage outlier matrix found for orthologs [%s]\n", results$norm_display))
      return(NULL)
    }
    outlier_matrix <- results$outliers$stage_outliers$orthologs_only
    title_suffix <- "Orthologs Only"
    distance_col <- "total_distance_scaled_ortholog"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  # Find families that are outliers in at least one stage
  families_with_outliers <- which(rowSums(outlier_matrix, na.rm = TRUE) > 0)
  
  if (length(families_with_outliers) == 0) {
    cat(sprintf("  No families with stage-specific outliers found for %s [%s]\n", 
                method, results$norm_display))
    return(NULL)
  }
  
  # Get top families by number of stages where they're outliers
  outlier_counts <- rowSums(outlier_matrix[families_with_outliers, , drop = FALSE], na.rm = TRUE)
  top_families <- families_with_outliers[order(outlier_counts, decreasing = TRUE)][1:min(n_top, length(families_with_outliers))]
  
  # Prepare data for heatmap
  heatmap_data <- data.frame()
  
  for (fam_id in top_families) {
    for (stage in colnames(outlier_matrix)) {
      if (stage %in% colnames(shift_matrix)) {
        shift_val <- shift_matrix[fam_id, stage]
        is_outlier <- outlier_matrix[fam_id, stage]
        
        if (!is.na(shift_val)) {
          total_dist <- results$trajectory[[distance_col]][results$trajectory$family_id == fam_id]
          heatmap_data <- rbind(heatmap_data, data.frame(
            family_id = fam_id,
            family_label = sprintf("F%d (%.1f)", fam_id, total_dist),
            stage = stage,
            shift = shift_val,
            is_outlier = is_outlier
          ))
        }
      }
    }
  }
  
  if (nrow(heatmap_data) == 0) {
    return(NULL)
  }
  
  # Order families by total distance
  family_order <- heatmap_data %>%
    group_by(family_id) %>%
    summarise(mean_abs_shift = mean(abs(shift), na.rm = TRUE)) %>%
    arrange(desc(mean_abs_shift)) %>%
    pull(family_id)
  
  heatmap_data <- heatmap_data %>%
    mutate(
      stage = factor(stage, levels = STAGES),
      family_id = factor(family_id, levels = family_order),
      family_label = factor(family_label, levels = unique(family_label[order(family_id)]))
    )
  
  p <- ggplot(heatmap_data, aes(x = stage, y = family_label)) +
    geom_tile(aes(fill = shift)) +
    geom_point(data = subset(heatmap_data, is_outlier), 
               aes(shape = "outlier"), size = 3, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Normalized\nShift") +
    scale_shape_manual(values = c("outlier" = 8), name = "") +
    labs(
      title = sprintf("Stage-Specific Outlier Families - %s", results$project_id),
      subtitle = sprintf("%s [%s]: Families with per-stage outliers (★ = outlier)\nTop %d families by outlier frequency", 
                        title_suffix, results$norm_display, length(top_families)),
      x = "Stage",
      y = "Family (Total Scaled Distance)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "right",
          legend.box = "vertical")
  
  return(p)
}

# ==================== TOP FAMILIES ====================

#' Show top outlier families across multiple consistency thresholds
#' Creates a 2x2 grid of plots with different consistency thresholds
plot_family_top_outliers <- function(results, method = "all", n_per_tail = 5,
                                     shift_threshold = 0.5) {
  
  # Define the four consistency thresholds to test
  consistency_thresholds <- c(0.2, 0.4, 0.6, 0.8)
  threshold_labels <- c("≥20% (any direction)", "≥40% (≥2 stages)", 
                        "≥60% (≥3 stages)", "≥80% (all stages)")
  
  if (method == "all") {
    shift_matrix <- results$stage_shifts$all_genes
    title_suffix <- "All Genes"
    distance_col <- "total_distance_scaled_all"
    outlier_col <- "is_total_outlier_all"
  } else if (method == "ortholog") {
    shift_matrix <- results$stage_shifts$orthologs_only
    title_suffix <- "Orthologs Only"
    distance_col <- "total_distance_scaled_ortholog"
    outlier_col <- "is_total_outlier_ortholog"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  # Get the top outlier families based on scaled total distance
  outlier_families <- results$trajectory %>%
    filter(.data[[outlier_col]]) %>%
    arrange(desc(.data[[distance_col]])) %>%
    pull(family_id)
  
  if (length(outlier_families) == 0) {
    cat(sprintf("  No outlier families found for %s [%s]\n", title_suffix, results$norm_display))
    return(NULL)
  }
  
  cat(sprintf("  Found %d outlier families total, will show top %d per threshold\n", 
              length(outlier_families), n_per_tail))
  
  # Create a list to store the four plots
  all_plots <- list()
  
  for (t_idx in seq_along(consistency_thresholds)) {
    consistency_threshold <- consistency_thresholds[t_idx]
    threshold_label <- threshold_labels[t_idx]
    
    # Find outlier families that meet the consistency threshold
    selected_families <- c()
    selected_means <- c()
    selected_consistency <- c()
    
    for (fam_id in outlier_families) {
      # Get shifts for this family across stages
      fam_shifts <- c()
      for (stage in STAGES) {
        if (stage %in% colnames(shift_matrix)) {
          shift_val <- shift_matrix[fam_id, stage]
          if (!is.na(shift_val)) {
            fam_shifts <- c(fam_shifts, shift_val)
          }
        }
      }
      
      if (length(fam_shifts) >= 2) {
        mean_abs_shift <- mean(abs(fam_shifts), na.rm = TRUE)
        n_stages <- length(fam_shifts)
        
        # Check consistency of direction
        n_positive <- sum(fam_shifts > 0, na.rm = TRUE)
        n_negative <- sum(fam_shifts < 0, na.rm = TRUE)
        consistency <- max(n_positive, n_negative) / n_stages
        
        if (consistency >= consistency_threshold) {
          selected_families <- c(selected_families, fam_id)
          selected_means <- c(selected_means, mean_abs_shift)
          selected_consistency <- c(selected_consistency, consistency)
        }
      }
    }
    
    # Take top N by mean absolute shift
    if (length(selected_families) > 0) {
      names(selected_means) <- selected_families
      top_families <- as.numeric(names(sort(selected_means, decreasing = TRUE)[1:min(n_per_tail, length(selected_families))]))
    } else {
      top_families <- c()
    }
    
    if (length(top_families) == 0) {
      # Create a placeholder plot when no families found
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, 
                 label = sprintf("No families at %.0f%% consistency", consistency_threshold*100),
                 size = 5, hjust = 0.5) +
        theme_void() +
        labs(title = sprintf("Consistency ≥%.0f%%", consistency_threshold*100),
             subtitle = threshold_label)
      all_plots[[t_idx]] <- p
      next
    }
    
    # Prepare data for this threshold
    plot_data <- data.frame()
    
    for (fam_id in top_families) {
      # Calculate mean absolute shift for this family
      fam_shifts <- c()
      for (stage in STAGES) {
        if (stage %in% colnames(shift_matrix)) {
          shift_val <- shift_matrix[fam_id, stage]
          if (!is.na(shift_val)) {
            fam_shifts <- c(fam_shifts, shift_val)
          }
        }
      }
      mean_abs_shift <- mean(abs(fam_shifts), na.rm = TRUE)
      
      # Calculate consistency for this family
      n_positive <- sum(fam_shifts > 0, na.rm = TRUE)
      n_negative <- sum(fam_shifts < 0, na.rm = TRUE)
      consistency_val <- max(n_positive, n_negative) / length(fam_shifts)
      direction <- ifelse(n_positive > n_negative, "Up", "Down")
      
      for (stage in STAGES) {
        if (stage %in% colnames(shift_matrix)) {
          shift_val <- shift_matrix[fam_id, stage]
          if (!is.na(shift_val)) {
            plot_data <- rbind(plot_data, data.frame(
              family_id = fam_id,
              family_label = sprintf("F%d (%.2f)", 
                                    fam_id, mean_abs_shift),
              stage = stage,
              shift = shift_val,
              stage_num = which(STAGES == stage),
              direction = direction
            ))
          }
        }
      }
    }
    
    if (nrow(plot_data) == 0) {
      next
    }
    
    # Create the plot for this threshold
    p <- ggplot(plot_data, aes(x = stage_num, y = shift, color = family_label, group = family_label)) +
      geom_line(size = 1.2) +
      geom_point(aes(shape = direction), size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
      scale_x_continuous(breaks = 1:length(STAGES), labels = STAGES) +
      scale_shape_manual(values = c("Up" = 16, "Down" = 17)) +
      scale_color_brewer(palette = "Set1") +
      labs(
        title = sprintf("Consistency ≥%.0f%% (n=%d)", consistency_threshold*100, length(top_families)),
        subtitle = threshold_label,
        x = "Stage",
        y = "Normalized Shift",
        color = "Family (mean)",
        shape = "Direction"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = 7),
            legend.box = "vertical",
            plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9))
    
    all_plots[[t_idx]] <- p
  }
  
  # Combine all four thresholds into a 2x2 grid
  if (length(all_plots) == 4) {
    final_plot <- wrap_plots(all_plots, ncol = 2) +
      plot_annotation(
        title = sprintf("Top Outlier Families by Consistency Threshold - %s", results$project_id),
        subtitle = sprintf("%s [%s] - Top %d outlier families (▲ up, ▼ down)", 
                          title_suffix, results$norm_display, n_per_tail),
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                     plot.subtitle = element_text(hjust = 0.5, size = 11))
      )
    return(final_plot)
  } else {
    return(NULL)
  }
}

# ==================== MEANINGFUL STAGE OUTLIER PERCENTAGES ====================

#' Show percentage of families that are outliers at each stage
plot_family_stage_outlier_percentages <- function(results, method = "all") {
  
  if (method == "all") {
    if (is.null(results$outliers) || is.null(results$outliers$stage_outliers$all_genes)) {
      return(NULL)
    }
    outlier_matrix <- results$outliers$stage_outliers$all_genes
    title_suffix <- "All Genes"
  } else if (method == "ortholog") {
    if (is.null(results$outliers) || is.null(results$outliers$stage_outliers$orthologs_only)) {
      return(NULL)
    }
    outlier_matrix <- results$outliers$stage_outliers$orthologs_only
    title_suffix <- "Orthologs Only"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  # Calculate percentage of outliers at each stage
  plot_data <- data.frame()
  
  for (stage in STAGES) {
    if (stage %in% colnames(outlier_matrix)) {
      n_outliers <- sum(outlier_matrix[, stage], na.rm = TRUE)
      total_families <- sum(!is.na(outlier_matrix[, stage]))
      
      if (total_families > 0) {
        plot_data <- rbind(plot_data, data.frame(
          stage = stage,
          n_outliers = n_outliers,
          total_families = total_families,
          percentage = n_outliers / total_families * 100
        ))
      }
    }
  }
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  plot_data$stage <- factor(plot_data$stage, levels = STAGES)
  
  p <- ggplot(plot_data, aes(x = stage, y = percentage)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percentage, n_outliers)), 
              vjust = -0.5, size = 3) +
    labs(
      title = sprintf("Stage-Specific Outlier Percentages - %s", results$project_id),
      subtitle = sprintf("%s [%s]: Percentage of families flagged as outliers at each stage\n(based on per-stage shift distributions)", 
                        title_suffix, results$norm_display),
      x = "Stage",
      y = "Families Flagged as Outliers (%)"
    ) +
    ylim(0, min(100, max(plot_data$percentage, na.rm = TRUE) * 1.2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  
  return(p)
}

# ==================== VELOCITY VISUALIZATIONS ====================

#' Multi-line plot showing velocity trajectories for top outlier families
plot_family_velocity_multiline <- function(results, method = "all", n_per_tail = 5) {
  
  if (method == "all") {
    if (is.null(results$velocities)) {
      return(NULL)
    }
    velocities <- results$velocities$all_genes
    outliers <- results$velocities$outliers_all
    title_suffix <- "All Genes"
    distance_col <- "total_distance_scaled_all"
  } else if (method == "ortholog") {
    if (is.null(results$velocities)) {
      return(NULL)
    }
    velocities <- results$velocities$orthologs_only
    outliers <- results$velocities$outliers_ortholog
    title_suffix <- "Orthologs Only"
    distance_col <- "total_distance_scaled_ortholog"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  if (is.null(velocities) || length(velocities) == 0) {
    return(NULL)
  }
  
  # Find families that are velocity outliers in any transition
  all_outlier_fams <- c()
  for (trans in names(outliers)) {
    all_outlier_fams <- c(all_outlier_fams, which(outliers[[trans]]))
  }
  all_outlier_fams <- unique(all_outlier_fams)
  
  if (length(all_outlier_fams) == 0) {
    cat(sprintf("  No velocity outlier families found for %s [%s]\n", 
                method, results$norm_display))
    return(NULL)
  }
  
  # Get top families by mean absolute velocity
  fam_velocities <- data.frame()
  for (fam_id in all_outlier_fams) {
    vel_values <- c()
    for (trans in names(velocities)) {
      if (!is.null(velocities[[trans]][fam_id]) && !is.na(velocities[[trans]][fam_id])) {
        vel_values <- c(vel_values, abs(velocities[[trans]][fam_id]))
      }
    }
    if (length(vel_values) > 0) {
      fam_velocities <- rbind(fam_velocities, data.frame(
        family_id = fam_id,
        mean_abs_vel = mean(vel_values, na.rm = TRUE)
      ))
    }
  }
  
  top_families <- fam_velocities %>%
    arrange(desc(mean_abs_vel)) %>%
    head(n_per_tail) %>%
    pull(family_id)
  
  # Prepare data for plotting
  plot_data <- data.frame()
  
  for (fam_id in top_families) {
    # Get velocity for each transition
    for (i in 1:(length(STAGES)-1)) {
      trans_name <- paste(STAGES[i], "->", STAGES[i+1])
      if (!is.null(velocities[[trans_name]])) {
        vel_val <- velocities[[trans_name]][fam_id]
        if (!is.na(vel_val)) {
          is_outlier_val <- outliers[[trans_name]][fam_id]
          total_dist <- results$trajectory[[distance_col]][results$trajectory$family_id == fam_id]
          plot_data <- rbind(plot_data, data.frame(
            family_id = fam_id,
            family_label = sprintf("F%d (%.1f)", fam_id, total_dist),
            transition = trans_name,
            transition_num = i,
            velocity = vel_val,
            is_outlier = is_outlier_val
          ))
        }
      }
    }
  }
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # Order families by total distance
  family_order <- plot_data %>%
    group_by(family_id) %>%
    summarise(mean_abs_vel = mean(abs(velocity), na.rm = TRUE)) %>%
    arrange(desc(mean_abs_vel)) %>%
    pull(family_id)
  
  plot_data <- plot_data %>%
    mutate(
      family_label = factor(family_label, levels = unique(family_label[order(family_id)])),
      transition = factor(transition, levels = names(velocities))
    )
  
  p <- ggplot(plot_data, aes(x = transition, y = velocity, color = family_label, group = family_label)) +
    geom_line(size = 1) +
    geom_point(aes(shape = is_outlier), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 16), 
                       labels = c("Outlier", "Normal")) +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = sprintf("Family Velocity Trajectories - %s", results$project_id),
      subtitle = sprintf("%s [%s]: Top %d families by mean |velocity|\n★ = velocity outlier in that transition", 
                        title_suffix, results$norm_display, length(top_families)),
      x = "Stage Transition",
      y = "Velocity (Δ Normalized Shift)",
      color = "Family",
      shape = "Status"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom",
          legend.box = "vertical")
  
  return(p)
}

#' Heatmap showing velocity values for top outlier families
plot_family_velocity_heatmap <- function(results, method = "all", n_top = 20) {
  
  if (method == "all") {
    if (is.null(results$velocities)) {
      return(NULL)
    }
    velocities <- results$velocities$all_genes
    outliers <- results$velocities$outliers_all
    title_suffix <- "All Genes"
    distance_col <- "total_distance_scaled_all"
  } else if (method == "ortholog") {
    if (is.null(results$velocities)) {
      return(NULL)
    }
    velocities <- results$velocities$orthologs_only
    outliers <- results$velocities$outliers_ortholog
    title_suffix <- "Orthologs Only"
    distance_col <- "total_distance_scaled_ortholog"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  if (is.null(velocities) || length(velocities) == 0) {
    return(NULL)
  }
  
  # Find families that are velocity outliers in any transition
  all_outlier_fams <- c()
  for (trans in names(outliers)) {
    all_outlier_fams <- c(all_outlier_fams, which(outliers[[trans]]))
  }
  all_outlier_fams <- unique(all_outlier_fams)
  
  if (length(all_outlier_fams) == 0) {
    return(NULL)
  }
  
  # Get top families by total distance
  top_families <- results$trajectory %>%
    filter(family_id %in% all_outlier_fams) %>%
    arrange(desc(.data[[distance_col]])) %>%
    head(n_top) %>%
    pull(family_id)
  
  # Prepare heatmap data
  heatmap_data <- data.frame()
  
  for (fam_id in top_families) {
    total_dist <- results$trajectory[[distance_col]][results$trajectory$family_id == fam_id]
    for (i in 1:(length(STAGES)-1)) {
      trans_name <- paste(STAGES[i], "->", STAGES[i+1])
      if (!is.null(velocities[[trans_name]])) {
        vel_val <- velocities[[trans_name]][fam_id]
        if (!is.na(vel_val)) {
          is_outlier_val <- outliers[[trans_name]][fam_id]
          heatmap_data <- rbind(heatmap_data, data.frame(
            family_label = sprintf("F%d (%.1f)", fam_id, total_dist),
            transition = trans_name,
            velocity = vel_val,
            is_outlier = is_outlier_val
          ))
        }
      }
    }
  }
  
  if (nrow(heatmap_data) == 0) {
    return(NULL)
  }
  
  # Order families by mean absolute velocity
  family_order <- heatmap_data %>%
    group_by(family_label) %>%
    summarise(mean_abs_vel = mean(abs(velocity), na.rm = TRUE)) %>%
    arrange(desc(mean_abs_vel)) %>%
    pull(family_label)
  
  heatmap_data <- heatmap_data %>%
    mutate(
      transition = factor(transition, levels = names(velocities)),
      family_label = factor(family_label, levels = family_order)
    )
  
  p <- ggplot(heatmap_data, aes(x = transition, y = family_label)) +
    geom_tile(aes(fill = velocity)) +
    geom_point(data = subset(heatmap_data, is_outlier), 
               aes(shape = "outlier"), size = 2, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Velocity") +
    scale_shape_manual(values = c("outlier" = 8), name = "") +
    labs(
      title = sprintf("Family Velocity Heatmap - %s", results$project_id),
      subtitle = sprintf("%s [%s]: Top %d velocity outlier families\n★ = velocity outlier in that transition", 
                        title_suffix, results$norm_display, length(top_families)),
      x = "Stage Transition",
      y = "Family (Total Scaled Distance)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "right")
  
  return(p)
}

#' Show percentage of families that are velocity outliers at each transition
plot_family_velocity_outlier_percentages <- function(results, method = "all") {
  
  if (method == "all") {
    if (is.null(results$velocities) || is.null(results$velocities$outliers_all)) {
      return(NULL)
    }
    outliers <- results$velocities$outliers_all
    title_suffix <- "All Genes"
  } else if (method == "ortholog") {
    if (is.null(results$velocities) || is.null(results$velocities$outliers_ortholog)) {
      return(NULL)
    }
    outliers <- results$velocities$outliers_ortholog
    title_suffix <- "Orthologs Only"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  plot_data <- data.frame()
  
  for (trans in names(outliers)) {
    n_outliers <- sum(outliers[[trans]], na.rm = TRUE)
    total_families <- length(outliers[[trans]])
    
    if (total_families > 0) {
      plot_data <- rbind(plot_data, data.frame(
        transition = trans,
        n_outliers = n_outliers,
        total_families = total_families,
        percentage = n_outliers / total_families * 100
      ))
    }
  }
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  plot_data$transition <- factor(plot_data$transition, levels = names(outliers))
  
  p <- ggplot(plot_data, aes(x = transition, y = percentage)) +
    geom_bar(stat = "identity", fill = "darkorange", alpha = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", percentage, n_outliers)), 
              vjust = -0.5, size = 3) +
    labs(
      title = sprintf("Velocity Outlier Percentages - %s", results$project_id),
      subtitle = sprintf("%s [%s]: Percentage of families with |velocity| > 95th percentile\n(thresholds determined per transition)", 
                        title_suffix, results$norm_display),
      x = "Stage Transition",
      y = "Families with Outlier Velocity (%)"
    ) +
    ylim(0, min(100, max(plot_data$percentage, na.rm = TRUE) * 1.2)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  
  return(p)
}

# ==================== LOESS SMOOTHING PLOTS ====================

plot_family_loess_smoothing <- function(results, method = "all") {
  
  if (is.null(results$loess_stats)) {
    cat(sprintf("  No LOESS stats available for [%s]\n", results$norm_display))
    return(NULL)
  }
  
  if (method == "all") {
    x_vals <- results$loess_stats$loess_x_all
    y_vals <- results$loess_stats$loess_y_all
    y_smoothed <- results$loess_stats$loess_y_smoothed_all
    family_means <- results$loess_stats$gene_family_means_all
    family_sds <- results$loess_stats$loess_y_all
    valid_families <- results$loess_stats$valid_families_all
    title_suffix <- "All Genes"
  } else if (method == "ortholog") {
    x_vals <- results$loess_stats$loess_x_ortholog
    y_vals <- results$loess_stats$loess_y_ortholog
    y_smoothed <- results$loess_stats$loess_y_smoothed_ortholog
    family_means <- results$loess_stats$gene_family_means_ortholog
    family_sds <- results$loess_stats$loess_y_ortholog
    valid_families <- results$loess_stats$valid_families_ortholog
    title_suffix <- "Orthologs Only"
  } else {
    stop("Method must be 'all' or 'ortholog'")
  }
  
  if (is.null(x_vals) || length(x_vals) == 0) {
    return(NULL)
  }
  
  valid_idx <- !is.na(x_vals) & !is.na(y_vals) & !is.na(y_smoothed) & 
              x_vals > 0 & y_smoothed > 0
  
  loess_data <- data.frame(
    x = x_vals[valid_idx],
    y_raw = y_vals[valid_idx],
    y_smoothed = y_smoothed[valid_idx]
  )
  
  # Sortiere nach x
  loess_data <- loess_data[order(loess_data$x), ]
  
  # Prepare plot_data (alle Familien)
  plot_data <- data.frame(
    family_means = family_means,
    family_sds = family_sds
  )
  
  plot_data <- plot_data[!is.na(plot_data$family_means) & !is.na(plot_data$family_sds) & 
                         plot_data$family_means > 0 & plot_data$family_sds > 0, ]
  
  p <- ggplot() +
    geom_point(data = plot_data, aes(x = family_means, y = family_sds), 
               alpha = 0.3, color = "grey50", size = 1) +
    geom_point(data = loess_data, aes(x = x, y = y_raw), 
               color = "steelblue", size = 2, alpha = 0.7) +
    geom_line(data = loess_data, aes(x = x, y = y_smoothed), 
              color = "red", size = 1.5) +
    labs(
      title = sprintf("LOESS Smoothing of Family Gene SDs - %s", results$project_id),
      subtitle = sprintf("%s [%s]: %d families used in LOESS\nSpan = %.2f, Degree = %d", 
                        title_suffix, results$norm_display,
                        length(valid_families),
                        LOESS_SPAN, LOESS_DEGREE),
      x = "Family Mean Expression",
      y = "Within-Family Gene SD"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  
  return(p)
}

# ==================== CREATE ALL PLOTS ====================

#' Create all family plots for a project across normalization methods
create_family_plots <- function(project_id, norm_methods = c("raw", "full", "std_log"), 
                                save_dir = NULL) {
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("CREATING FAMILY PLOTS FOR: %s\n", project_id))
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  # Load results for all normalization methods
  results_list <- list()
  for (norm_method in norm_methods) {
    res <- load_family_results_norm(project_id, norm_method)
    if (!is.null(res)) {
      results_list[[norm_method]] <- res
    }
  }
  
  if (length(results_list) == 0) {
    cat(sprintf("  No family results found for %s\n", project_id))
    return(NULL)
  }
  
  all_plots <- list()
  
  # Create plots for each normalization method
  for (norm_method in names(results_list)) {
    results <- results_list[[norm_method]]
    
    cat(sprintf("\n  Processing [%s]...\n", results$norm_display))
    
    # Create method-specific save directory
    if (is.null(save_dir)) {
      norm_dir <- get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method)
      plot_dir <- file.path(norm_dir, project_id, "plots")
    } else {
      plot_dir <- file.path(save_dir, norm_method)
    }
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    method_plots <- list()
    
    # 1. Distribution plots
    cat("    Creating distribution plots...\n")
    p1 <- plot_family_scaled_total_distance(results, "all")
    p2 <- plot_family_scaled_total_distance(results, "ortholog")
    
    if (!is.null(p1)) {
      ggsave(file.path(plot_dir, "01_family_total_distance_scaled_all.png"), p1, 
             width = 10, height = 6, dpi = 300)
      method_plots$total_distance_scaled_all <- p1
    }
    if (!is.null(p2)) {
      ggsave(file.path(plot_dir, "02_family_total_distance_scaled_ortholog.png"), p2, 
             width = 10, height = 6, dpi = 300)
      method_plots$total_distance_scaled_ortholog <- p2
    }
    
    # 2. Stage outlier heatmap (based on per-stage outliers)
    cat("    Creating stage outlier heatmaps...\n")
    p3 <- plot_family_stage_outlier_heatmap(results, "all")
    p4 <- plot_family_stage_outlier_heatmap(results, "ortholog")
    
    if (!is.null(p3)) {
      ggsave(file.path(plot_dir, "03_family_stage_outlier_heatmap_all.png"), p3, 
             width = 12, height = 10, dpi = 300)
      method_plots$stage_heatmap_all <- p3
    }
    if (!is.null(p4)) {
      ggsave(file.path(plot_dir, "04_family_stage_outlier_heatmap_ortholog.png"), p4, 
             width = 12, height = 10, dpi = 300)
      method_plots$stage_heatmap_ortholog <- p4
    }
    
    # 3. Top example trajectories
    cat("    Creating top family example trajectories...\n")
    p5 <- plot_family_top_outliers(results, "all")
    p6 <- plot_family_top_outliers(results, "ortholog")
    
    if (!is.null(p5)) {
      ggsave(file.path(plot_dir, "05_top_family_examples_all.png"), p5, 
             width = 14, height = 8, dpi = 300)
      method_plots$top_families_all <- p5
    }
    if (!is.null(p6)) {
      ggsave(file.path(plot_dir, "06_top_family_examples_ortholog.png"), p6, 
             width = 14, height = 8, dpi = 300)
      method_plots$top_families_ortholog <- p6
    }
    
    # 4. Stage outlier percentages (meaningful variation)
    cat("    Creating stage outlier percentage plots...\n")
    p7 <- plot_family_stage_outlier_percentages(results, "all")
    p8 <- plot_family_stage_outlier_percentages(results, "ortholog")
    
    if (!is.null(p7)) {
      ggsave(file.path(plot_dir, "07_family_stage_outlier_percentages_all.png"), p7, 
             width = 10, height = 7, dpi = 300)
      method_plots$stage_percentages_all <- p7
    }
    if (!is.null(p8)) {
      ggsave(file.path(plot_dir, "08_family_stage_outlier_percentages_ortholog.png"), p8, 
             width = 10, height = 7, dpi = 300)
      method_plots$stage_percentages_ortholog <- p8
    }
    
    # 5. Velocity multi-line plots
    cat("    Creating velocity multi-line plots...\n")
    p9 <- plot_family_velocity_multiline(results, "all")
    p10 <- plot_family_velocity_multiline(results, "ortholog")
    
    if (!is.null(p9)) {
      ggsave(file.path(plot_dir, "09_family_velocity_multiline_all.png"), p9, 
             width = 12, height = 8, dpi = 300)
      method_plots$velocity_multiline_all <- p9
    }
    if (!is.null(p10)) {
      ggsave(file.path(plot_dir, "10_family_velocity_multiline_ortholog.png"), p10, 
             width = 12, height = 8, dpi = 300)
      method_plots$velocity_multiline_ortholog <- p10
    }
    
    # 6. Velocity heatmaps
    cat("    Creating velocity heatmaps...\n")
    p11 <- plot_family_velocity_heatmap(results, "all")
    p12 <- plot_family_velocity_heatmap(results, "ortholog")
    
    if (!is.null(p11)) {
      ggsave(file.path(plot_dir, "11_family_velocity_heatmap_all.png"), p11, 
             width = 10, height = 10, dpi = 300)
      method_plots$velocity_heatmap_all <- p11
    }
    if (!is.null(p12)) {
      ggsave(file.path(plot_dir, "12_family_velocity_heatmap_ortholog.png"), p12, 
             width = 10, height = 10, dpi = 300)
      method_plots$velocity_heatmap_ortholog <- p12
    }
    
    # 7. Velocity outlier percentages (meaningful variation)
    cat("    Creating velocity percentage plots...\n")
    p13 <- plot_family_velocity_outlier_percentages(results, "all")
    p14 <- plot_family_velocity_outlier_percentages(results, "ortholog")
    
    if (!is.null(p13)) {
      ggsave(file.path(plot_dir, "13_family_velocity_percentages_all.png"), p13, 
             width = 10, height = 7, dpi = 300)
      method_plots$velocity_percentages_all <- p13
    }
    if (!is.null(p14)) {
      ggsave(file.path(plot_dir, "14_family_velocity_percentages_ortholog.png"), p14, 
             width = 10, height = 7, dpi = 300)
      method_plots$velocity_percentages_ortholog <- p14
    }
    
    # 8. LOESS smoothing plots
    cat("    Creating LOESS smoothing plots...\n")
    p15 <- plot_family_loess_smoothing(results, "all")
    p16 <- plot_family_loess_smoothing(results, "ortholog")
    
    if (!is.null(p15)) {
      ggsave(file.path(plot_dir, "15_family_loess_smoothing_all.png"), p15, 
             width = 10, height = 7, dpi = 300)
      method_plots$loess_all <- p15
    }
    if (!is.null(p16)) {
      ggsave(file.path(plot_dir, "16_family_loess_smoothing_ortholog.png"), p16, 
             width = 10, height = 7, dpi = 300)
      method_plots$loess_ortholog <- p16
    }
    
    all_plots[[norm_method]] <- method_plots
    cat(sprintf("    ✓ Created %d plots for [%s]\n", length(method_plots), results$norm_display))
  }
  
  cat(sprintf("\n  ✓ All family plots complete for %s\n", project_id))
  
  return(all_plots)
}

# ==================== MAIN EXECUTION ====================

if (sys.nframe() == 0) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("FAMILY TRAJECTORY ANALYSIS - IMPROVED PLOTTING\n")
  cat("Top families | Velocity visualizations | Stage-specific heatmaps\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  start_time <- Sys.time()
  
  # Create plots for all projects
  all_project_plots <- list()
  
  for (project_id in PROJECTS) {
    tryCatch({
      plots <- create_family_plots(project_id, norm_methods = c("raw", "full", "std_log", "log"))
      if (!is.null(plots)) {
        all_project_plots[[project_id]] <- plots
      }
    }, error = function(e) {
      cat(sprintf("\n  ✗ Error plotting %s: %s\n", project_id, e$message))
      print(e)
    })
  }
  
  # Summary
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("FAMILY PLOTTING COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("Total time: %.1f minutes\n", duration))
  cat(sprintf("Projects processed: %d/%d\n", length(all_project_plots), length(PROJECTS)))
  cat("\nOutput directories:\n")
  cat("  - Raw data:       ", get_norm_output_dir(FAMILY_OUTPUT_DIR, "raw"), "/[PROJECT]/plots/\n")
  cat("  - Full norm:      ", get_norm_output_dir(FAMILY_OUTPUT_DIR, "full"), "/[PROJECT]/plots/\n")
  cat("  - Std+log norm:   ", get_norm_output_dir(FAMILY_OUTPUT_DIR, "std_log"), "/[PROJECT]/plots/\n")
  cat("  - Log norm:       ", get_norm_output_dir(FAMILY_OUTPUT_DIR, "log"), "/[PROJECT]/plots/\n")
}
