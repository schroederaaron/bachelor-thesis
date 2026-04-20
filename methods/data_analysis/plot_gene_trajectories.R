# plot_gene_trajectories.R
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(gridExtra)
library(RColorBrewer)
library(colorspace)
library(ggrepel)

# ==================== CONFIGURATION ====================

source("config.R")
source("utils.R")

# Helper function to get short normalization display name
get_short_norm_name <- function(norm_method) {
  switch(norm_method,
    "raw" = "raw",
    "log" = "log",
    "std_log" = "std_log",
    "full" = "full",
    "gene_wise_scaling_quantile_log2" = "std_log",
    "gene_wise_scaling_quantile" = "full",
    norm_method
  )
}

SHIFT_LABEL <- "Normalized signed difference"
SHIFT_LABEL_ABS <- "|Normalized signed difference|"
VELOCITY_LABEL <- "Δ normalized signed difference"

# ==================== LOAD RESULTS ====================

compute_fixed_label_positions <- function(n_labels, y_limits, span_fraction = 0.50) {
  if (!all(is.finite(y_limits))) y_limits <- c(-1, 1)
  y_range <- diff(y_limits)
  if (!is.finite(y_range) || y_range == 0) y_range <- max(abs(y_limits), na.rm = TRUE) * 2
  if (!is.finite(y_range) || y_range == 0) y_range <- 2

  y_center <- mean(y_limits)
  half_span <- (y_range * span_fraction) / 2

  if (n_labels <= 1) {
    return(y_center)
  }

  seq(y_center + half_span, y_center - half_span, length.out = n_labels)
}

load_gene_results_all <- function(project_id, norm_methods = c("raw", "full", "std_log", "log")) {
  results_list <- list()
  for (norm_method in norm_methods) {
    res <- load_gene_results_norm(project_id, norm_method)
    if (!is.null(res)) {
      results_list[[norm_method]] <- res
    }
  }
  return(results_list)
}

# define themes
theme_top_multiline <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9, color = "grey20"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.y = grid::unit(0.03, "cm"),
      plot.margin = margin(8, 12, 8, 8)
    )
}

# ==================== DISTRIBUTION PLOTS ====================

plot_gene_scaled_total_distance <- function(results, comp_type = "own") {
  if (comp_type == "own") {
    distance_col <- "total_distance_scaled_own"
    outlier_col <- "is_outlier_own"
    title_suffix <- "Own"
  } else if (comp_type == "fam") {
    distance_col <- "total_distance_scaled_fam"
    outlier_col <- "is_outlier_fam"
    title_suffix <- "Family"
  } else if (comp_type == "orth") {
    distance_col <- "total_distance_scaled_orth"
    outlier_col <- "is_outlier_orth"
    title_suffix <- "Ortholog"
  } else {
    stop("comp_type must be 'own', 'fam', or 'orth'")
  }
  
  if (!distance_col %in% names(results$trajectory)) return(NULL)
  
  plot_data <- results$trajectory[!is.na(results$trajectory[[distance_col]]), ]
  if (nrow(plot_data) == 0) return(NULL)
  
  percentile_95 <- quantile(plot_data[[distance_col]], 0.95, na.rm = TRUE)
  n_outliers <- sum(plot_data[[outlier_col]], na.rm = TRUE)
  total_genes <- nrow(plot_data)
  
  short_norm <- get_short_norm_name(results$norm_method)
  
  p <- ggplot(plot_data, aes(x = .data[[distance_col]])) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = percentile_95, color = "red", linetype = "dashed", size = 1) +
    annotate("text", x = percentile_95 * 1.1, y = Inf, 
             label = paste("95th %ile =", round(percentile_95, 2)),
             vjust = 2, hjust = 0, size = 3, color = "red") +
    labs(
      title = sprintf("%s [%s] %s", results$project_id, short_norm, title_suffix),
      subtitle = sprintf("%d outliers / %d genes", n_outliers, total_genes),
      x = "Scaled Total Distance",
      y = "Count"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  
  return(p)
}

# ==================== TOP OUTLIER GENES ====================

plot_gene_top_outliers <- function(results, comp_type = "own", n_genes = 5, color_manager = NULL) {
  
  if (comp_type == "own") {
    matrix_name   <- "gene_vs_own_healthy"
    title_suffix  <- "Own"
    distance_col  <- "total_distance_scaled_own"
    outlier_col   <- "is_outlier_own"
  } else if (comp_type == "fam") {
    matrix_name   <- "gene_vs_family_mean"
    title_suffix  <- "Family"
    distance_col  <- "total_distance_scaled_fam"
    outlier_col   <- "is_outlier_fam"
  } else if (comp_type == "orth") {
    matrix_name   <- "gene_vs_family_ortholog_mean"
    title_suffix  <- "Ortholog"
    distance_col  <- "total_distance_scaled_orth"
    outlier_col   <- "is_outlier_orth"
  } else {
    stop("comp_type must be 'own', 'fam', or 'orth'")
  }
  
  if (is.null(results$stage_shifts$shifts[[matrix_name]])) return(NULL)
  
  shift_matrix <- results$stage_shifts$shifts[[matrix_name]]
  
  top_tbl <- results$trajectory %>%
    filter(.data[[outlier_col]]) %>%
    arrange(desc(.data[[distance_col]])) %>%
    slice_head(n = n_genes) %>%
    select(gene_id, ranking_value = all_of(distance_col))
  
  if (nrow(top_tbl) == 0) return(NULL)
  
  plot_data <- lapply(top_tbl$gene_id, function(gene_id) {
    gene_idx <- which(results$trajectory$gene_id == gene_id)
    if (length(gene_idx) == 0) return(NULL)
    
    gene_shifts <- shift_matrix[gene_idx, , drop = TRUE]
    gene_shifts <- gene_shifts[!is.na(gene_shifts)]
    if (length(gene_shifts) == 0) return(NULL)
    
    mean_abs_shift <- mean(abs(gene_shifts), na.rm = TRUE)
    
    do.call(rbind, lapply(seq_along(STAGES), function(i) {
      stage <- STAGES[i]
      if (!stage %in% colnames(shift_matrix)) return(NULL)
      shift_val <- shift_matrix[gene_idx, stage]
      if (is.na(shift_val)) return(NULL)
      
      data.frame(
        gene_id = gene_id,
        stage = stage,
        stage_num = i,
        shift = shift_val,
        mean_abs_shift = mean_abs_shift,
        stringsAsFactors = FALSE
      )
    }))
  }) %>% bind_rows()
  
  if (nrow(plot_data) == 0) return(NULL)
  
  legend_tbl <- plot_data %>%
    group_by(gene_id) %>%
    summarise(mean_abs_shift = first(mean_abs_shift), .groups = "drop") %>%
    arrange(desc(mean_abs_shift), gene_id) %>%
    mutate(legend_label = sprintf("%s (%.2f)", gene_id, mean_abs_shift))
  
  plot_data <- plot_data %>%
    left_join(legend_tbl, by = "gene_id")
  
  plot_data$legend_label <- factor(plot_data$legend_label, levels = legend_tbl$legend_label)
  
  short_norm <- get_short_norm_name(results$norm_method)
  
  if (!is.null(color_manager)) {
    gene_ids_in_plot <- legend_tbl$gene_id
    gene_colors <- color_manager$get_colors(gene_ids_in_plot)
    color_values <- setNames(unname(gene_colors[legend_tbl$gene_id]), legend_tbl$legend_label)
  } else {
    pal <- scales::hue_pal()(nrow(legend_tbl))
    color_values <- setNames(pal, legend_tbl$legend_label)
  }
  
  p <- ggplot(plot_data, aes(x = stage_num, y = shift, color = legend_label, group = gene_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_line(linewidth = 1.1, alpha = 0.95) +
    geom_point(size = 2.3) +
    scale_x_continuous(
      breaks = seq_along(STAGES),
      labels = STAGES,
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_color_manual(values = color_values, drop = FALSE) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 1.1, size = 2.3))) +
    labs(
      title = sprintf("%s [%s] %s", results$project_id, short_norm, title_suffix),
      subtitle = sprintf("Top %d outliers by mean absolute %s across stages
Legend: Gene ID (mean %s)",
                         n_genes, tolower(SHIFT_LABEL), SHIFT_LABEL_ABS),
      x = "Stage",
      y = SHIFT_LABEL,
      color = "Gene ID"
    ) +
    theme_top_multiline()
  
  return(p)
}

# ==================== STAGE-SPECIFIC GENE HEATMAP ====================

plot_gene_stage_outlier_heatmap <- function(results, comp_type = "own", n_top = 20) {
  if (comp_type == "own") {
    matrix_name <- "gene_vs_own_healthy"
    title_suffix <- "Own"
    distance_col <- "total_distance_scaled_own"
  } else if (comp_type == "fam") {
    matrix_name <- "gene_vs_family_mean"
    title_suffix <- "Family"
    distance_col <- "total_distance_scaled_fam"
  } else if (comp_type == "orth") {
    matrix_name <- "gene_vs_family_ortholog_mean"
    title_suffix <- "Ortholog"
    distance_col <- "total_distance_scaled_orth"
  } else {
    stop("comp_type must be 'own', 'fam', or 'orth'")
  }
  
  if (is.null(results$stage_shifts$shifts[[matrix_name]])) return(NULL)
  
  shift_matrix <- results$stage_shifts$shifts[[matrix_name]]
  
  if (is.null(results$outliers) || 
      is.null(results$outliers$stage_specific_outliers[[matrix_name]])) return(NULL)
  
  stage_outliers <- results$outliers$stage_specific_outliers[[matrix_name]]
  
  outlier_matrix <- matrix(FALSE, nrow = nrow(shift_matrix), ncol = length(STAGES))
  colnames(outlier_matrix) <- STAGES
  rownames(outlier_matrix) <- rownames(shift_matrix)
  
  for (stage in STAGES) {
    if (!is.null(stage_outliers[[stage]])) {
      outlier_matrix[, stage] <- stage_outliers[[stage]]$is_outlier
    }
  }
  
  genes_with_outliers <- which(rowSums(outlier_matrix, na.rm = TRUE) > 0)
  if (length(genes_with_outliers) == 0) return(NULL)
  
  # Calculate mean absolute shift for each gene with outliers
  gene_mean_shifts <- data.frame(
    gene_id = results$trajectory$gene_id[genes_with_outliers],
    mean_abs_shift = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(genes_with_outliers)) {
    gene_idx <- genes_with_outliers[i]
    gene_shifts <- shift_matrix[gene_idx, ]
    gene_shifts <- gene_shifts[!is.na(gene_shifts)]
    gene_mean_shifts$mean_abs_shift[i] <- mean(abs(gene_shifts), na.rm = TRUE)
  }
  
  # Select top genes based on mean absolute shift (same as trajectory plot)
  top_genes <- gene_mean_shifts %>%
    arrange(desc(mean_abs_shift)) %>%
    head(n_top) %>%
    pull(gene_id)
  
  # Build heatmap data
  heatmap_data <- data.frame()
  for (gene_id in top_genes) {
    gene_idx <- which(results$trajectory$gene_id == gene_id)
    gene_display <- ifelse(nchar(gene_id) > 15, 
                          paste0(substr(gene_id, 1, 12), "..."), 
                          gene_id)
    
    # Get mean absolute shift for this gene
    mean_shift <- gene_mean_shifts$mean_abs_shift[gene_mean_shifts$gene_id == gene_id]
    
    for (stage in STAGES) {
      if (stage %in% colnames(shift_matrix)) {
        shift_val <- shift_matrix[gene_idx, stage]
        if (!is.na(shift_val)) {
          is_outlier <- outlier_matrix[gene_idx, stage]
          heatmap_data <- rbind(heatmap_data, data.frame(
            gene_label = sprintf("%s (%.2f)", gene_display, mean_shift),
            stage = stage,
            shift = shift_val,
            is_outlier = is_outlier,
            mean_abs_shift = mean_shift,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  if (nrow(heatmap_data) == 0) return(NULL)
  
  # Order by mean absolute shift (same as trajectory plot)
  gene_order <- heatmap_data %>%
    group_by(gene_label) %>%
    summarise(mean_abs_shift = unique(mean_abs_shift)) %>%
    arrange(desc(mean_abs_shift)) %>%
    pull(gene_label)
  
  heatmap_data <- heatmap_data %>%
    mutate(stage = factor(stage, levels = STAGES),
           gene_label = factor(gene_label, levels = gene_order))
  
  short_norm <- get_short_norm_name(results$norm_method)
  
  p <- ggplot(heatmap_data, aes(x = stage, y = gene_label)) +
    geom_tile(aes(fill = shift)) +
    geom_point(data = subset(heatmap_data, is_outlier), 
               aes(shape = "outlier"), size = 2, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Norm. signed diff.") +
    scale_shape_manual(values = c("outlier" = 8), name = "") +
    labs(
      title = sprintf("%s [%s] %s", results$project_id, short_norm, title_suffix),
      subtitle = sprintf("Top %d genes by mean absolute %s", length(top_genes), tolower(SHIFT_LABEL)),
      x = "Stage",
      y = "Gene"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "right")
  
  return(p)
}

# ==================== VELOCITY VISUALIZATIONS ====================

plot_gene_velocity_multiline <- function(results, comp_type = "own", n_genes = 5,
                                         color_manager = NULL) {
  
  if (comp_type == "own") {
    matrix_name <- "gene_vs_own_healthy"
    title_suffix <- "Own"
  } else if (comp_type == "fam") {
    matrix_name <- "gene_vs_family_mean"
    title_suffix <- "Family"
  } else if (comp_type == "orth") {
    matrix_name <- "gene_vs_family_ortholog_mean"
    title_suffix <- "Ortholog"
  } else {
    stop("comp_type must be 'own', 'fam', or 'orth'")
  }
  
  if (is.null(results$stage_shifts$shifts[[matrix_name]])) return(NULL)
  
  shift_matrix <- results$stage_shifts$shifts[[matrix_name]]
  
  outlier_genes <- rownames(shift_matrix)[apply(!is.na(shift_matrix), 1, any)]
  if (length(outlier_genes) == 0) return(NULL)
  
  gene_scores <- lapply(outlier_genes, function(gene_id) {
    gene_idx <- which(rownames(shift_matrix) == gene_id)
    if (length(gene_idx) == 0) return(NULL)
    gene_shifts <- shift_matrix[gene_idx, , drop = TRUE]
    vel_vals <- c(gene_shifts[1], diff(gene_shifts))
    data.frame(gene_id = gene_id, mean_abs_vel = mean(abs(vel_vals), na.rm = TRUE), stringsAsFactors = FALSE)
  }) %>% bind_rows()
  
  if (nrow(gene_scores) == 0) return(NULL)
  
  top_genes <- gene_scores %>%
    arrange(desc(mean_abs_vel), gene_id) %>%
    slice_head(n = n_genes) %>%
    pull(gene_id)
  
  plot_data <- lapply(top_genes, function(gene_id) {
    gene_idx <- which(rownames(shift_matrix) == gene_id)
    if (length(gene_idx) == 0) return(NULL)
    
    gene_shifts <- shift_matrix[gene_idx, , drop = TRUE]
    mean_abs_vel <- mean(abs(c(gene_shifts[1], diff(gene_shifts))), na.rm = TRUE)
    gene_rows <- list()
    
    for (i in seq_along(gene_shifts)) {
      vel_val <- if (i == 1) gene_shifts[1] else gene_shifts[i] - gene_shifts[i - 1]
      if (is.na(vel_val)) next
      
      transition_label <- if (i == 1) "H→S1" else paste0("S", i - 1, "→S", i)
      gene_rows[[length(gene_rows) + 1]] <- data.frame(
        gene_id = gene_id,
        transition = transition_label,
        transition_num = i,
        velocity = vel_val,
        mean_abs_vel = mean_abs_vel,
        is_outlier = TRUE,
        stringsAsFactors = FALSE
      )
    }
    
    bind_rows(gene_rows)
  }) %>% bind_rows()
  
  if (nrow(plot_data) == 0) return(NULL)
  
  all_transitions <- c("H→S1", paste0("S", 1:(length(STAGES) - 1), "→S", 2:length(STAGES)))
  plot_data$transition <- factor(plot_data$transition, levels = all_transitions)
  
  legend_tbl <- plot_data %>%
    group_by(gene_id) %>%
    summarise(mean_abs_vel = mean(abs(velocity), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_abs_vel), gene_id) %>%
    mutate(legend_label = sprintf("%s (%.2f)", gene_id, mean_abs_vel))
  
  plot_data <- plot_data %>%
    left_join(legend_tbl, by = "gene_id")
  
  plot_data$legend_label <- factor(plot_data$legend_label, levels = legend_tbl$legend_label)
  
  short_norm <- get_short_norm_name(results$norm_method)
  
  if (!is.null(color_manager)) {
    gene_ids_in_plot <- legend_tbl$gene_id
    gene_colors <- color_manager$get_colors(gene_ids_in_plot)
    color_values <- setNames(unname(gene_colors[legend_tbl$gene_id]), legend_tbl$legend_label)
  } else {
    pal <- scales::hue_pal()(nrow(legend_tbl))
    color_values <- setNames(pal, legend_tbl$legend_label)
  }
  
  p <- ggplot(plot_data, aes(x = transition_num, y = velocity, color = legend_label, group = gene_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.5) +
    geom_line(linewidth = 1.1, alpha = 0.95) +
    geom_point(size = 2.4, stroke = 0.6) +
    scale_x_continuous(
      breaks = seq_along(all_transitions),
      labels = all_transitions,
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_color_manual(values = color_values, drop = FALSE) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 1.1, size = 2.4))) +
    labs(
      title = sprintf("%s [%s] %s", results$project_id, short_norm, title_suffix),
      subtitle = sprintf("Top %d genes by mean absolute %s
Legend: Gene ID (mean |%s|)",
                         n_genes, VELOCITY_LABEL, VELOCITY_LABEL),
      x = "Transition",
      y = VELOCITY_LABEL,
      color = "Gene ID"
    ) +
    theme_top_multiline() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5)
    )
  
  return(p)
}

# ==================== VELOCITY HEATMAP ====================

plot_gene_velocity_heatmap <- function(results, comp_type = "own", n_top = 30) {
  
  # Select the appropriate velocity data
  if (comp_type == "own") {
    matrix_name <- "gene_vs_own_healthy"
    title_suffix <- "Own"
    distance_col <- "total_distance_scaled_own"
  } else if (comp_type == "fam") {
    matrix_name <- "gene_vs_family_mean"
    title_suffix <- "Family"
    distance_col <- "total_distance_scaled_fam"
  } else if (comp_type == "orth") {
    matrix_name <- "gene_vs_family_ortholog_mean"
    title_suffix <- "Ortholog"
    distance_col <- "total_distance_scaled_orth"
  } else {
    stop("comp_type must be 'own', 'fam', or 'orth'")
  }
  
  if (is.null(results$velocities) || 
      is.null(results$velocities$velocities[[matrix_name]])) {
    return(NULL)
  }
  
  velocities <- results$velocities$velocities[[matrix_name]]
  outliers <- results$velocities$outliers[[matrix_name]]
  
  if (is.null(velocities) || length(velocities) == 0) {
    return(NULL)
  }
  
  # Get shift matrix for Healthy -> Stage1 transition
  shift_matrix <- results$stage_shifts$shifts[[matrix_name]]
  
  # Get stage-specific outlier information to determine if Healthy->Stage1 is outlier
  stage_outliers <- NULL
  if (!is.null(results$outliers) && 
      !is.null(results$outliers$stage_specific_outliers[[matrix_name]])) {
    stage_outliers <- results$outliers$stage_specific_outliers[[matrix_name]]
  }
  
  # Find genes that are velocity outliers in any transition
  all_outlier_genes <- c()
  for (trans in names(outliers)) {
    all_outlier_genes <- c(all_outlier_genes, which(outliers[[trans]]))
  }
  all_outlier_genes <- unique(all_outlier_genes)
  
  if (length(all_outlier_genes) == 0) {
    return(NULL)
  }
  
  # Get top genes by scaled total distance
  top_genes <- results$trajectory %>%
    filter(row_number() %in% all_outlier_genes) %>%
    arrange(desc(.data[[distance_col]])) %>%
    head(n_top) %>%
    pull(gene_id)
  
  # Map stage names to S1, S2, etc.
  stage_map <- setNames(1:length(STAGES), STAGES)
  
  # Prepare heatmap data
  heatmap_data <- data.frame()
  
  for (gene_id in top_genes) {
    gene_idx <- which(results$trajectory$gene_id == gene_id)
    gene_display <- ifelse(nchar(gene_id) > 15, 
                          paste0(substr(gene_id, 1, 12), "..."), 
                          gene_id)
    total_dist <- results$trajectory[[distance_col]][gene_idx]
    
    # Add Healthy -> Stage1 transition
    stage1_shift <- shift_matrix[gene_idx, STAGES[1]]
    if (!is.na(stage1_shift)) {
      is_outlier_stage1 <- FALSE
      if (!is.null(stage_outliers) && !is.null(stage_outliers[[STAGES[1]]])) {
        is_outlier_stage1 <- stage_outliers[[STAGES[1]]]$is_outlier[gene_idx]
      }
      
      heatmap_data <- rbind(heatmap_data, data.frame(
        gene_label = sprintf("%s (%.1f)", gene_display, total_dist),
        transition = "H -> S1",
        velocity = stage1_shift,
        is_outlier = is_outlier_stage1
      ))
    }
    
    # Add existing transitions with simplified names
    for (i in 1:(length(STAGES)-1)) {
      trans_name_original <- paste(STAGES[i], "->", STAGES[i+1])
      trans_name_simplified <- paste0("S", stage_map[STAGES[i]], " -> S", stage_map[STAGES[i+1]])
      
      if (!is.null(velocities[[trans_name_original]])) {
        vel_val <- velocities[[trans_name_original]][gene_idx]
        if (!is.na(vel_val)) {
          is_outlier_val <- outliers[[trans_name_original]][gene_idx]
          heatmap_data <- rbind(heatmap_data, data.frame(
            gene_label = sprintf("%s (%.1f)", gene_display, total_dist),
            transition = trans_name_simplified,
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
  
  # Define factor levels including new transition
  all_transitions <- c("H -> S1", 
                       paste0("S", 1:(length(STAGES)-1), " -> S", 2:length(STAGES)))
  
  # Order genes by total distance
  gene_order <- heatmap_data %>%
    group_by(gene_label) %>%
    summarise(mean_abs_vel = mean(abs(velocity), na.rm = TRUE)) %>%
    arrange(desc(mean_abs_vel)) %>%
    pull(gene_label)
  
  heatmap_data <- heatmap_data %>%
    mutate(
      transition = factor(transition, levels = all_transitions),
      gene_label = factor(gene_label, levels = gene_order)
    )
  
  short_norm <- get_short_norm_name(results$norm_method)
  
  p <- ggplot(heatmap_data, aes(x = transition, y = gene_label)) +
    geom_tile(aes(fill = velocity)) +
    geom_point(data = subset(heatmap_data, is_outlier), 
               aes(shape = "outlier"), size = 2, color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = "Δ norm. signed diff.") +
    scale_shape_manual(values = c("outlier" = 8), name = "Outlier") +
    labs(
      title = sprintf("%s [%s] %s", results$project_id, short_norm, title_suffix),
      subtitle = sprintf("Top %d genes by mean absolute %s", length(top_genes), VELOCITY_LABEL),
      x = "Transition",
      y = "Gene"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 11),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "right")
  
  return(p)
}

# ==================== JACCARD SIMILARITY FUNCTIONS ====================

#' Calculate Jaccard similarity between two outlier sets
#' @param outliers1 Logical vector for method 1
#' @param outliers2 Logical vector for method 2
#' @return Jaccard similarity coefficient
calculate_jaccard <- function(outliers1, outliers2) {
  if (length(outliers1) != length(outliers2)) {
    return(NA)
  }
  
  intersection <- sum(outliers1 & outliers2, na.rm = TRUE)
  union <- sum(outliers1 | outliers2, na.rm = TRUE)
  
  if (union == 0) {
    return(0)
  }
  
  return(intersection / union)
}

#' Calculate Jaccard similarities between all pairs of normalization methods
#' for a given project and comparison type
#' 
#' @param results_list List of results for different normalization methods
#' @param comp_type Comparison type ("own", "fam", or "orth")
#' @return Matrix of Jaccard similarities
calculate_jaccard_matrix <- function(results_list, comp_type = "own") {
  
  # Map comparison type to outlier column
  outlier_col <- paste0("is_outlier_", comp_type)
  
  # Extract outlier vectors for each normalization method
  outlier_vectors <- list()
  norm_names <- c()
  
  for (norm_method in names(results_list)) {
    res <- results_list[[norm_method]]
    
    if (!is.null(res$trajectory) && outlier_col %in% names(res$trajectory)) {
      outlier_vectors[[norm_method]] <- res$trajectory[[outlier_col]]
      norm_names <- c(norm_names, norm_method)
    }
  }
  
  if (length(outlier_vectors) < 2) {
    return(NULL)
  }
  
  # Calculate pairwise Jaccard similarities
  n_methods <- length(outlier_vectors)
  jaccard_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods)
  rownames(jaccard_matrix) <- norm_names
  colnames(jaccard_matrix) <- norm_names
  
  for (i in 1:n_methods) {
    for (j in 1:n_methods) {
      if (i == j) {
        jaccard_matrix[i, j] <- 1
      } else {
        jaccard_matrix[i, j] <- calculate_jaccard(outlier_vectors[[i]], outlier_vectors[[j]])
      }
    }
  }
  
  return(jaccard_matrix)
}

#' Create combined Jaccard heatmap for a cancer project
#' Shows similarity between normalization methods across all comparison types
#' 
#' @param results_list List of results for different normalization methods
#' @param project_id Project ID
#' @param output_dir Output directory for saving plots
plot_combined_jaccard_heatmap <- function(results_list, project_id, output_dir) {
  
  comp_types <- c("own", "fam", "orth")
  comp_labels <- c("Gene vs Own", "Gene vs Family", "Gene vs Ortholog")
  
  # Calculate Jaccard matrices for each comparison type
  jaccard_matrices <- list()
  valid_comps <- c()
  
  for (i in seq_along(comp_types)) {
    comp <- comp_types[i]
    jmat <- calculate_jaccard_matrix(results_list, comp)
    
    if (!is.null(jmat)) {
      jaccard_matrices[[comp]] <- jmat
      valid_comps <- c(valid_comps, comp)
    }
  }
  
  if (length(jaccard_matrices) == 0) {
    cat("    No Jaccard matrices calculated\n")
    return(NULL)
  }
  
  # Prepare data for combined heatmap
  all_data <- data.frame()
  
  for (comp in names(jaccard_matrices)) {
    jmat <- jaccard_matrices[[comp]]
    
    # Convert to long format
    for (i in 1:nrow(jmat)) {
      for (j in 1:ncol(jmat)) {
        all_data <- rbind(all_data, data.frame(
          Method1 = rownames(jmat)[i],
          Method2 = colnames(jmat)[j],
          Similarity = jmat[i, j],
          Comparison = comp
        ))
      }
    }
  }
  
  # Create label mappings
  norm_labels <- c(
    "raw" = "raw",
    "log" = "log",
    "std_log" = "std_log",
    "full" = "full"
  )
  
  all_data$Method1_label <- norm_labels[all_data$Method1]
  all_data$Method2_label <- norm_labels[all_data$Method2]
  
  # Create factor levels for methods
  method_order <- c("raw", "log", "std_log", "full")
  all_data$Method1_label <- factor(all_data$Method1_label, levels = norm_labels[method_order])
  all_data$Method2_label <- factor(all_data$Method2_label, levels = norm_labels[method_order])
  
  # Comparison labels
  comp_labels_full <- c(
    "own" = "Own Healthy",
    "fam" = "Family Mean",
    "orth" = "Ortholog Mean"
  )
  all_data$Comparison_label <- comp_labels_full[all_data$Comparison]
  
  # Create individual heatmaps for each comparison type and combine
  plots <- list()
  
  for (comp in names(jaccard_matrices)) {
    plot_data <- all_data[all_data$Comparison == comp, ]
    
    p <- ggplot(plot_data, aes(x = Method1_label, y = Method2_label, fill = Similarity)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%.2f", Similarity)), size = 3) +
      scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkred",
                          midpoint = 0.5, limits = c(0, 1),
                          name = "Jaccard\nSimilarity") +
      labs(
        title = comp_labels_full[comp],
        x = "",
        y = ""
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        legend.position = "right",
        panel.grid = element_blank()
      ) +
      coord_fixed()
    
    plots[[comp]] <- p
  }
  
  # Combine plots based on number of valid comparisons
  if (length(plots) == 3) {
    combined <- (plots[[1]] | plots[[2]] | plots[[3]]) +
      plot_annotation(
        title = sprintf("Jaccard Similarity: Normalization Methods - %s", project_id),
        subtitle = "Outlier overlap between normalization methods",
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      )
    width <- 15
    height <- 6
  } else if (length(plots) == 2) {
    combined <- (plots[[1]] | plots[[2]]) +
      plot_annotation(
        title = sprintf("Jaccard Similarity: Normalization Methods - %s", project_id),
        subtitle = "Outlier overlap between normalization methods",
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      )
    width <- 10
    height <- 6
  } else {
    combined <- plots[[1]] +
      plot_annotation(
        title = sprintf("Jaccard Similarity: Normalization Methods - %s", project_id),
        subtitle = "Outlier overlap between normalization methods",
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      )
    width <- 6
    height <- 5
  }
  
  # Save combined plot
  plot_dir <- file.path(output_dir, "jaccard_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  filename <- sprintf("jaccard_combined_%s.png", project_id)
  ggsave(file.path(plot_dir, filename), combined, width = width, height = height, dpi = 300)
  cat(sprintf("    Saved: %s\n", filename))
  
  # Also save individual matrices as CSV
  for (comp in names(jaccard_matrices)) {
    csv_filename <- sprintf("jaccard_matrix_%s_%s.csv", project_id, comp)
    write.csv(jaccard_matrices[[comp]], file.path(plot_dir, csv_filename))
  }
  
  return(combined)
}

# ==================== MULTI-CANCER OUTLIER HEATMAP ====================

#' Create a multi-cancer outlier heatmap using distance-based outliers
plot_multicancer_outlier_heatmap <- function(all_results_list, norm_method, comp_type, 
                                              output_dir, n_top = 30, sort_by = "n_stages") {
  
  library(data.table)
  
  # Map comp_type to matrix name
  if (comp_type == "own") {
    matrix_name <- "gene_vs_own_healthy"
    title_suffix <- "Own"
  } else if (comp_type == "fam") {
    matrix_name <- "gene_vs_family_mean"
    title_suffix <- "Family"
  } else if (comp_type == "orth") {
    matrix_name <- "gene_vs_family_ortholog_mean"
    title_suffix <- "Ortholog"
  } else {
    stop("comp_type must be 'own', 'fam', or 'orth'")
  }
  
  cat(sprintf("\n  Creating outlier heatmap: %s [%s] (sorted by %s)\n", 
              norm_method, title_suffix, 
              ifelse(sort_by == "n_stages", "n stages", paste("mean", SHIFT_LABEL_ABS))))
  
  # Collect data for all projects
  all_dt_list <- list()
  project_order <- c()
  
  for (project_id in names(all_results_list)) {
    if (!norm_method %in% names(all_results_list[[project_id]])) next
    
    results <- all_results_list[[project_id]][[norm_method]]
    
    if (is.null(results$outliers) || 
        is.null(results$outliers$stage_specific_outliers[[matrix_name]])) {
      next
    }
    
    stage_outliers <- results$outliers$stage_specific_outliers[[matrix_name]]
    shift_matrix <- results$stage_shifts$shifts[[matrix_name]]
    gene_ids <- results$trajectory$gene_id
    
    dt_list <- list()
    for (stage_idx in seq_along(STAGES)) {
      stage <- STAGES[stage_idx]
      if (!is.null(stage_outliers[[stage]])) {
        is_outlier <- stage_outliers[[stage]]$is_outlier
        
        if (!is.null(shift_matrix) && stage %in% colnames(shift_matrix)) {
          shift_vals <- shift_matrix[, stage]
          if (is.matrix(shift_vals)) {
            shift_vals <- as.vector(shift_vals)
          }
        } else {
          shift_vals <- rep(0, length(is_outlier))
        }
        
        direction <- rep(0, length(is_outlier))
        direction[is_outlier & shift_vals > 0] <- 1
        direction[is_outlier & shift_vals < 0] <- -1
        shift_abs <- abs(shift_vals)
        
        dt <- data.table(
          gene_id = gene_ids,
          stage = stage,
          direction = direction,
          is_outlier = is_outlier,
          shift_abs = shift_abs
        )
        dt_list[[stage_idx]] <- dt
      }
    }
    
    if (length(dt_list) > 0) {
      project_dt <- rbindlist(dt_list)
      project_dt[, project := project_id]
      all_dt_list[[length(all_dt_list) + 1]] <- project_dt
      project_order <- c(project_order, project_id)
      cat(sprintf("    Processed %s: %d genes\n", project_id, length(gene_ids)))
    }
  }
  
  if (length(all_dt_list) == 0) {
    cat("    No outlier data found\n")
    return(NULL)
  }
  
  all_data <- rbindlist(all_dt_list)
  
  # Calculate statistics per gene
  mean_shifts <- all_data[, .(mean_abs_shift = mean(shift_abs, na.rm = TRUE)), by = gene_id]
  gene_stats <- all_data[is_outlier == TRUE, .(total_outlier_stages = .N), by = gene_id]
  
  all_genes <- unique(all_data$gene_id)
  gene_stats <- data.table(gene_id = all_genes)
  gene_stats <- merge(gene_stats, mean_shifts, by = "gene_id", all.x = TRUE)
  gene_stats <- merge(gene_stats, 
                      all_data[is_outlier == TRUE, .(total_outlier_stages = .N), by = gene_id], 
                      by = "gene_id", all.x = TRUE)
  gene_stats[is.na(total_outlier_stages), total_outlier_stages := 0]
  gene_stats[is.na(mean_abs_shift), mean_abs_shift := 0]
  
  # Sort and select top genes
  gene_stats_sorted <- copy(gene_stats)
  if (sort_by == "n_stages") {
    setorder(gene_stats_sorted, -total_outlier_stages, -mean_abs_shift)
  } else {
    setorder(gene_stats_sorted, -mean_abs_shift, -total_outlier_stages)
  }
  
  top_genes <- gene_stats_sorted[1:min(n_top, nrow(gene_stats_sorted)), gene_id]
  
  # Prepare heatmap data
  heatmap_data <- all_data[gene_id %in% top_genes]
  heatmap_data <- merge(heatmap_data, gene_stats, by = "gene_id")
  
  heatmap_data[, gene_label := sprintf("%s (%.2f)", 
                                       ifelse(nchar(gene_id) > 15, 
                                              paste0(substr(gene_id, 1, 12), "..."), 
                                              gene_id),
                                       mean_abs_shift)]
  heatmap_data[, stage := factor(stage, levels = STAGES)]
  heatmap_data[, col_label := paste(project, stage, sep = "\n")]
  
  # Order genes
  gene_order_dt <- unique(heatmap_data[, .(gene_id, gene_label, total_outlier_stages, mean_abs_shift)])
  if (sort_by == "n_stages") {
    setorder(gene_order_dt, -total_outlier_stages, -mean_abs_shift)
  } else {
    setorder(gene_order_dt, -mean_abs_shift, -total_outlier_stages)
  }
  heatmap_data[, gene_label := factor(gene_label, levels = gene_order_dt$gene_label)]
  
  # Order columns
  col_order <- CJ(project = project_order, stage = STAGES)[, col_label := paste(project, stage, sep = "\n")]$col_label
  heatmap_data[, col_label := factor(col_label, levels = col_order)]
  
  # Create plot
  p <- ggplot(heatmap_data, aes(x = col_label, y = gene_label, fill = factor(direction))) +
    geom_tile(color = "white", size = 0.3) +
    scale_fill_manual(
      values = c("1" = "firebrick", "-1" = "steelblue", "0" = "grey90"),
      name = "Outlier direction",
      labels = c("1" = "Up", "-1" = "Down", "0" = "Not Outlier")
    ) +
    labs(
      title = sprintf("Multi-Cancer Outliers [%s] %s", norm_method, title_suffix),
      subtitle = sprintf("Top %d genes by %s", n_top, 
                        ifelse(sort_by == "n_stages", "# stages", SHIFT_LABEL_ABS)),
      x = "Cancer Type and Stage",
      y = "Gene"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 6),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  # Save
  output_subdir <- file.path(output_dir, "multicancer_heatmaps_uncorrected")
  dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
  
  sort_suffix <- ifelse(sort_by == "n_stages", "bynstages", "bymeanshift")
  filename <- sprintf("multicancer_outliers_%s_%s_%s_top%d.png", 
                      norm_method, comp_type, sort_suffix, n_top)
  height <- max(6, min(15, 6 + n_top * 0.12))
  width <- max(8, length(col_order) * 0.4)
  ggsave(file.path(output_subdir, filename), p, width = width, height = height, dpi = 300)
  cat(sprintf("    Saved: %s\n", filename))
  
  # Save CSV
  csv_data <- dcast(heatmap_data[, .(gene_label, col_label, direction, total_outlier_stages, mean_abs_shift)], 
                    gene_label + total_outlier_stages + mean_abs_shift ~ col_label, 
                    value.var = "direction", fill = 0)
  
  if (sort_by == "n_stages") {
    setorder(csv_data, -total_outlier_stages, -mean_abs_shift)
  } else {
    setorder(csv_data, -mean_abs_shift, -total_outlier_stages)
  }
  
  csv_filename <- sprintf("multicancer_outliers_%s_%s_%s_top%d.csv", 
                          norm_method, comp_type, sort_suffix, n_top)
  fwrite(csv_data, file.path(output_subdir, csv_filename))
  
  return(p)
}

#' Create all multi-cancer outlier heatmaps
create_all_multicancer_outlier_heatmaps <- function(all_results_list, output_dir, n_top = 30) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CREATING MULTI-CANCER OUTLIER HEATMAPS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  norm_methods <- c()
  for (project_id in names(all_results_list)) {
    norm_methods <- unique(c(norm_methods, names(all_results_list[[project_id]])))
  }
  
  comp_types <- c("own", "fam", "orth")
  sort_options <- c("n_stages", "mean_shift")
  total_plots <- 0
  
  for (norm_method in norm_methods) {
    for (comp_type in comp_types) {
      for (sort_by in sort_options) {
        p <- plot_multicancer_outlier_heatmap(
          all_results_list = all_results_list,
          norm_method = norm_method,
          comp_type = comp_type,
          output_dir = output_dir,
          n_top = n_top,
          sort_by = sort_by
        )
        if (!is.null(p)) total_plots <- total_plots + 1
      }
    }
  }
  
  cat(sprintf("\n  Total plots created: %d\n", total_plots))
}

# ==================== CREATE ALL GENE PLOTS ====================

create_gene_plots <- function(project_id, norm_methods = c("raw", "full", "std_log", "log"), 
                              save_dir = NULL) {
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("CREATING GENE PLOTS FOR: %s\n", project_id))
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  results_list <- load_gene_results_all(project_id, norm_methods)
  
  if (length(results_list) == 0) {
    cat(sprintf("  No gene results found for %s\n", project_id))
    return(NULL)
  }
  color_map_file <- file.path(dirname(GENE_OUTPUT_DIR), project_id, "gene_color_mapping.rds")
  existing_color_map <- NULL
  if (file.exists(color_map_file)) {
    existing_color_map <- readRDS(color_map_file)
    cat(sprintf("  Loaded existing color mapping: %d genes\n", length(existing_color_map)))
  }

  color_manager <- create_gene_color_manager(initial_map = existing_color_map)
  cat("  Initialized color manager\n")
  color_manager$diagnose_colors()
  
  # Process each normalization method
  for (norm_method in names(results_list)) {
    results <- results_list[[norm_method]]
    cat(sprintf("\n  Processing [%s]...\n", norm_method))
    
    if (is.null(save_dir)) {
      norm_dir <- get_norm_output_dir(GENE_OUTPUT_DIR, norm_method)
      plot_dir <- file.path(norm_dir, project_id, "plots")
    } else {
      plot_dir <- file.path(save_dir, norm_method)
    }
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    # 1. Distance distributions
    for (comp in c("own", "fam", "orth")) {
      p <- plot_gene_scaled_total_distance(results, comp)
      if (!is.null(p)) {
        ggsave(file.path(plot_dir, paste0("01_distance_", comp, ".png")), 
               p, width = 5, height = 4, dpi = 300)
      }
    }
    
    # 2. Top outliers
    for (comp in c("own", "fam", "orth")) {
      p <- plot_gene_top_outliers(results, comp, n_genes = 5, color_manager = color_manager)
      if (!is.null(p)) {
        ggsave(file.path(plot_dir, paste0("02_top_outliers_", comp, ".png")), 
               p, width = 6, height = 5, dpi = 300)
      }
    }
    
    # 3. Stage heatmaps
    for (comp in c("own", "fam", "orth")) {
      p <- plot_gene_stage_outlier_heatmap(results, comp, n_top = 20)
      if (!is.null(p)) {
        height <- max(6, min(15, 6 + 20 * 0.12))
        ggsave(file.path(plot_dir, paste0("03_stage_heatmap_", comp, ".png")), 
               p, width = 6, height = height, dpi = 300)
      }
    }
    
    # 4. Velocity multi-line
    for (comp in c("own", "fam", "orth")) {
      p <- plot_gene_velocity_multiline(results, comp, n_genes = 5, color_manager = color_manager)
      if (!is.null(p)) {
        ggsave(file.path(plot_dir, paste0("04_velocity_lines_", comp, ".png")), 
               p, width = 6, height = 5, dpi = 300)
      }
    }
    
    # 5. Velocity heatmaps
    for (comp in c("own", "fam", "orth")) {
      p <- plot_gene_velocity_heatmap(results, comp, n_top = 30)
      if (!is.null(p)) {
        height <- max(6, min(15, 6 + 20 * 0.12))
        ggsave(file.path(plot_dir, paste0("05_velocity_heatmap_", comp, ".png")), 
               p, width = 6, height = height, dpi = 300)
      }
    }
  }
  
  # Create Jaccard plot ONCE (after processing all normalizations)
  cat("\n    Creating Jaccard similarity heatmaps...\n")
  jaccard_dir <- file.path(dirname(GENE_OUTPUT_DIR), "jaccard_plots")
  p_jaccard <- plot_combined_jaccard_heatmap(results_list, project_id, jaccard_dir)
  # Save color mapping
  color_manager$save_map(color_map_file)
  cat(sprintf("\n  Saved color mapping: %d genes\n", length(color_manager$get_full_map())))
  
  cat(sprintf("\n  Plots complete for %s\n", project_id))
  
  return(invisible(TRUE))
}

# ==================== MAIN EXECUTION ====================

if (sys.nframe() == 0) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("GENE TRAJECTORY ANALYSIS - COMPLETE PLOTTING SUITE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  start_time <- Sys.time()
  
  if (!requireNamespace("colorspace", quietly = TRUE)) {
    cat("Installing colorspace package...\n")
    install.packages("colorspace")
  }
  
  # Collect results for multi-cancer heatmaps
  all_results_for_heatmap <- list()
  all_project_plots <- list()
  
  # Process each project
  for (project_id in PROJECTS) {
    cat(sprintf("\n>>> Processing project: %s\n", project_id))
    
    # Load results for multi-cancer heatmap
    project_results <- list()
    for (norm_method in c("raw", "full", "std_log", "log")) {
      res <- load_gene_results_norm(project_id, norm_method)
      if (!is.null(res)) {
        project_results[[norm_method]] <- res
      }
    }
    if (length(project_results) > 0) {
      all_results_for_heatmap[[project_id]] <- project_results
    }
    
    # Create individual project plots
    tryCatch({
      plots <- create_gene_plots(project_id, norm_methods = c("raw", "full", "std_log", "log"))
      if (!is.null(plots)) {
        all_project_plots[[project_id]] <- plots
      }
    }, error = function(e) {
      cat(sprintf("\n  Error plotting %s: %s\n", project_id, e$message))
    })
  }
  
  # Create multi-cancer outlier heatmaps
  if (length(all_results_for_heatmap) > 0) {
    output_base_dir <- dirname(GENE_OUTPUT_DIR)
    create_all_multicancer_outlier_heatmaps(all_results_for_heatmap, output_base_dir, n_top = 30)
  }
  
  # Summary
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("GENE PLOTTING COMPLETE\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("Total time: %.1f minutes\n", duration))
  cat(sprintf("Projects processed: %d/%d\n", length(all_project_plots), length(PROJECTS)))
}