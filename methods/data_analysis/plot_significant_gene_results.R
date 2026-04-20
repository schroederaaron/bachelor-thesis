#!/usr/bin/env Rscript
# plot_significant_gene_results.R
# Gene Trajectory Analysis for SIGNIFICANT Genes (noise_p_adj < 0.01)
# Replicates plot structure from plot_gene_trajectories.R
# + Comparison plots of signed difference distributions
# + Consistent colors from main analysis
# + Capped (Top 20) and uncapped (all) heatmaps for significant genes
# + Multi-cancer significance heatmaps

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(data.table)
library(ggrepel)

source("config.R")
source("utils.R")

# ==================== CONFIGURATION ====================

Q_CUTOFF <- 0.01

COMP_TYPES <- c("own_healthy", "family_mean", "ortholog_mean")
COMP_DISPLAY <- c(
  "own_healthy" = "Gene vs Own Healthy",
  "family_mean" = "Gene vs Family Mean",
  "ortholog_mean" = "Gene vs Ortholog Mean"
)

# Mapping from significance output comparison names to gene trajectory output
COMP_TYPE_MAPPING <- c(
  "own_healthy" = "own",
  "family_mean" = "fam",
  "ortholog_mean" = "orth"
)

SHIFT_LABEL <- "Normalized signed difference"
SHIFT_LABEL_ABS <- "|Normalized signed difference|"
SHIFT_LEGEND_SHORT <- "Norm. signed diff."

#' Load saved color mapping from main analysis (if available)
#' @param project_id Project ID
#' @return Named character vector with gene colors, or NULL if not found
load_gene_color_map <- function(project_id) {
  # Try multiple possible locations for color mapping file
  possible_paths <- c(
    file.path(dirname(FAMILY_OUTPUT_DIR), project_id, "gene_color_mapping.rds"),
    file.path(dirname(GENE_OUTPUT_DIR), project_id, "gene_color_mapping.rds"),
    file.path("output", "gene_trajectory_analysis", project_id, "gene_color_mapping.rds")
  )
  
  for (color_file in possible_paths) {
    if (file.exists(color_file)) {
      color_map <- readRDS(color_file)
      cat(sprintf("  Loaded color mapping for %d genes from: %s\n", 
                  length(color_map), basename(dirname(dirname(color_file)))))
      return(color_map)
    }
  }
  
  cat("  No saved color mapping found, will create new colors\n")
  return(NULL)
}

#' Create or load color manager for a project
#' @param project_id Project ID
#' @param use_saved_colors Whether to try loading saved colors from main analysis
#' @return Color manager list
create_or_load_color_manager <- function(project_id, use_saved_colors = TRUE) {
  base_color_manager <- create_gene_color_manager()
  saved_colors <- NULL

  if (use_saved_colors) {
    saved_colors <- load_gene_color_map(project_id)
    if (!is.null(saved_colors) && length(saved_colors) > 0) {
      saved_colors <- saved_colors[!is.na(names(saved_colors)) & names(saved_colors) != ""]
      cat(sprintf("  Using %d saved colors when available\n", length(saved_colors)))
    }
  }

  color_state <- new.env(parent = emptyenv())
  color_state$color_map <- if (!is.null(saved_colors)) saved_colors else character(0)

  MIN_COLOR_DISTANCE <- 23
  MIN_LUMINANCE_DISTANCE <- 8

  get_lab_coords <- function(colors) {
    colors <- as.character(colors)
    if (length(colors) == 0) {
      return(matrix(numeric(0), ncol = 3,
                    dimnames = list(NULL, c("L", "A", "B"))))
    }

    rgb <- t(grDevices::col2rgb(colors)) / 255
    lab <- grDevices::convertColor(rgb, from = "sRGB", to = "Lab", scale.in = 1)
    colnames(lab) <- c("L", "A", "B")
    lab
  }

  perceptual_distance <- function(color1, color2) {
    if (requireNamespace("farver", quietly = TRUE)) {
      lab1 <- farver::decode_colour(color1, to = "lab")
      lab2 <- farver::decode_colour(color2, to = "lab")
      return(sqrt(sum((lab1[1, ] - lab2[1, ])^2)))
    }

    lab <- get_lab_coords(c(color1, color2))
    sqrt(sum((lab[1, ] - lab[2, ])^2))
  }

  select_distinct_color <- function(existing_colors) {
    candidate_grid <- expand.grid(
      h = seq(0, 355, by = 5),
      c = c(50, 65, 80, 95),
      l = c(35, 45, 55, 65, 75),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )

    candidates <- unique(grDevices::hcl(
      h = candidate_grid$h,
      c = candidate_grid$c,
      l = candidate_grid$l,
      fixup = TRUE
    ))
    candidates <- candidates[!(candidates %in% existing_colors)]

    if (length(candidates) == 0) {
      stop("No candidate colors available for assignment")
    }

    if (length(existing_colors) == 0) {
      return(candidates[1])
    }

    existing_lab <- get_lab_coords(existing_colors)

    score_tbl <- lapply(candidates, function(candidate) {
      dists <- vapply(existing_colors, function(existing) {
        perceptual_distance(candidate, existing)
      }, numeric(1))
      candidate_lab <- get_lab_coords(candidate)
      l_diffs <- abs(existing_lab[, "L"] - candidate_lab[1, "L"])
      data.frame(
        candidate = candidate,
        min_dist = min(dists),
        mean_dist = mean(dists),
        min_l_diff = min(l_diffs),
        stringsAsFactors = FALSE
      )
    })
    score_tbl <- do.call(rbind, score_tbl)

    strong_candidates <- score_tbl[
      score_tbl$min_dist >= MIN_COLOR_DISTANCE &
        score_tbl$min_l_diff >= MIN_LUMINANCE_DISTANCE,
      , drop = FALSE
    ]

    if (nrow(strong_candidates) > 0) {
      strong_candidates <- strong_candidates[
        order(-strong_candidates$min_dist,
              -strong_candidates$mean_dist,
              -strong_candidates$min_l_diff,
              strong_candidates$candidate),
        , drop = FALSE
      ]
      return(strong_candidates$candidate[[1]])
    }

    fallback_tbl <- score_tbl[
      order(-score_tbl$min_dist,
            -score_tbl$mean_dist,
            -score_tbl$min_l_diff,
            score_tbl$candidate),
      , drop = FALSE
    ]

    warning(sprintf(
      "Could not satisfy ΔE >= %.1f and ΔL* >= %.1f; using best available color with ΔE = %.2f",
      MIN_COLOR_DISTANCE,
      MIN_LUMINANCE_DISTANCE,
      fallback_tbl$min_dist[[1]]
    ))

    fallback_tbl$candidate[[1]]
  }

  get_color <- function(gene_id) {
    gene_id <- as.character(gene_id)[1]
    if (is.na(gene_id) || gene_id == "") return("#CCCCCC")

    if (gene_id %in% names(color_state$color_map)) {
      return(unname(color_state$color_map[[gene_id]]))
    }

    existing_colors <- unname(color_state$color_map)
    new_color <- select_distinct_color(existing_colors)
    color_state$color_map[[gene_id]] <- new_color
    return(unname(new_color))
  }

  get_colors <- function(gene_ids) {
    gene_ids <- as.character(gene_ids)
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
    cols <- vapply(gene_ids, get_color, FUN.VALUE = character(1))
    names(cols) <- gene_ids
    cols
  }

  get_full_map <- function() {
    unname_map <- unlist(color_state$color_map)
    names(unname_map) <- names(color_state$color_map)
    unname_map
  }

  save_map <- function(output_dir = NULL) {
    if (is.null(output_dir) || is.na(output_dir) || output_dir == "") return(invisible(FALSE))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(get_full_map(), file.path(output_dir, "gene_color_mapping.rds"))
    invisible(TRUE)
  }

  reset_colors <- function() {
    color_state$color_map <- character(0)
  }

  diagnose_colors <- function(threshold = 23) {
    cols <- get_full_map()
    if (length(cols) < 2) {
      cat(sprintf("Color manager: %d colors assigned
", length(cols)))
      return(invisible(NULL))
    }
    cat(sprintf("Color manager: %d colors assigned
", length(cols)))
    for (i in seq_len(length(cols) - 1)) {
      for (j in (i + 1):length(cols)) {
        dist <- perceptual_distance(cols[[i]], cols[[j]])
        if (dist < threshold) {
          cat(sprintf("  Similar colors: %s vs %s (ΔE=%.2f)
", names(cols)[i], names(cols)[j], dist))
        }
      }
    }
    invisible(cols)
  }

  return(list(
    get_color = get_color,
    get_colors = get_colors,
    get_full_map = get_full_map,
    save_map = save_map,
    reset_colors = reset_colors,
    diagnose_colors = diagnose_colors,
    fix_similar_colors = function(...) invisible(FALSE)
  ))
}

# ==================== DATA LOADING ====================

#' Load significant genes from outlier significance analysis

load_significant_genes <- function() {
  sig_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")

  data_file <- file.path(sig_dir, "all_genes_pvalues.rds")

  if (!file.exists(data_file)) {
    stop(sprintf("Data file not found: %s\\n", data_file))
  }

  all_data <- readRDS(data_file)

  # Harmonize adjusted p-value column name across older/newer outputs
  if (!"noise_p_adj" %in% names(all_data) && "noise_p_value_adj" %in% names(all_data)) {
    all_data <- all_data %>% rename(noise_p_adj = noise_p_value_adj)
  }

  if (!"noise_p_adj" %in% names(all_data)) {
    stop("Expected adjusted p-value column 'noise_p_adj' (or legacy 'noise_p_value_adj') not found")
  }

  all_data <- all_data %>%
    mutate(
      is_outlier_stage = !is.na(noise_p_adj) &
        noise_p_adj < Q_CUTOFF &
        !is.na(signed_difference),
      outlier_marker = ifelse(is_outlier_stage, "*", "")
    )

  cat(sprintf("Loaded %d rows from all_genes_pvalues.rds\\n", nrow(all_data)))
  cat(sprintf("Stage-specific outlier rows: %d\\n", sum(all_data$is_outlier_stage, na.rm = TRUE)))

  return(all_data)
}

get_significant_rows <- function(data) {
  data %>% filter(is_outlier_stage)
}


#' Load ALL genes (uncorrected p-values) for comparison plots
load_all_genes <- function() {
  sig_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")
  
  data_file <- file.path(sig_dir, "all_genes_pvalues.rds")
  
  if (!file.exists(data_file)) {
    stop(sprintf("Data file not found: %s\n", data_file))
  }
  
  all_data <- readRDS(data_file)
  all_data <- all_data %>%
    mutate(
      is_outlier_stage = !is.na(noise_p_value_adj) &
        noise_p_value_adj < Q_CUTOFF &
        !is.na(signed_difference),
      outlier_marker = ifelse(is_outlier_stage, "*", "")
    )
  cat(sprintf("Loaded all genes: %d rows\n", nrow(all_data)))
  
  return(all_data)
}

# ==================== COMPARISON PLOTS: SIGNIFICANT vs ALL ====================

#' Plot comparing signed difference distributions
#' Left: Density plot, Right: Boxplot/Violin plot
#' Uses log2 scale for x-axis
plot_signed_diff_comparison <- function(all_data, sig_data = NULL, project_id, norm_method, 
                                         comp_type, output_dir) {
  
  all_filtered <- all_data %>%
    filter(norm_method == !!norm_method, comparison == comp_type, is_outlier_stage)
  
  if (is.null(sig_data)) {
    sig_filtered <- all_filtered %>%
      filter(is_outlier_stage)
  } else {
    sig_filtered <- sig_data %>%
      filter(norm_method == !!norm_method, comparison == comp_type, is_outlier_stage)
  }
  
  if (project_id != "all") {
    all_filtered <- all_filtered %>% filter(cancer_id == project_id)
    sig_filtered <- sig_filtered %>% filter(cancer_id == project_id)
  }
  
  if (nrow(all_filtered) == 0) {
    return(NULL)
  }
  
  # Calculate statistics for annotation
  all_stats <- all_filtered %>%
    summarise(
      n = n(),
      mean = mean(abs(signed_difference), na.rm = TRUE),
      median = median(abs(signed_difference), na.rm = TRUE),
      q95 = quantile(abs(signed_difference), 0.95, na.rm = TRUE)
    )
  
  sig_stats <- if (nrow(sig_filtered) > 0) {
    sig_filtered %>%
      summarise(
        n = n_distinct(gene_id),
        mean = mean(abs(signed_difference), na.rm = TRUE),
        median = median(abs(signed_difference), na.rm = TRUE),
        q95 = quantile(abs(signed_difference), 0.95, na.rm = TRUE)
      )
  } else {
    data.frame(n = 0, mean = NA, median = NA, q95 = NA)
  }
  
  # Prepare data frames for plots
  plot_data_all <- all_filtered %>%
    select(signed_difference) %>%
    mutate(type = "All Genes", abs_diff = abs(signed_difference))
  
  plot_data_sig <- sig_filtered %>%
    select(signed_difference) %>%
    mutate(type = "Significant Genes", abs_diff = abs(signed_difference))
  
  plot_data <- bind_rows(plot_data_all, plot_data_sig) %>%
    mutate(abs_diff_log2 = log2(abs_diff + 1e-10))
  
  # Title and subtitle
  if (project_id == "all") {
    title <- "Normalized Signed Difference Distribution - All Projects Combined"
    subtitle <- sprintf("%s | %s | log2 scale", COMP_DISPLAY[comp_type], norm_method)
  } else {
    title <- sprintf("Normalized Signed Difference Distribution - %s", project_id)
    subtitle <- sprintf("%s | %s | log2 scale", COMP_DISPLAY[comp_type], norm_method)
  }
  
  # Left plot: Density with log2 x-axis
  p_abs <- ggplot(plot_data, aes(x = abs_diff_log2, fill = type, color = type)) +
    geom_density(alpha = 0.5, adjust = 1.5) +
    scale_fill_manual(values = c("All Genes" = "grey70", "Significant Genes" = "steelblue")) +
    scale_color_manual(values = c("All Genes" = "grey50", "Significant Genes" = "steelblue")) +
    labs(
      title = title,
      subtitle = subtitle,
      x = sprintf("%s (log2 scale)", SHIFT_LABEL_ABS),
      y = "Density",
      fill = "Gene Set",
      color = "Gene Set"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "bottom"
    )
  
  # Right plot: Boxplot/Violin with log2 y-axis
  p_box <- ggplot(plot_data, aes(x = type, y = abs_diff_log2, fill = type)) +
    geom_violin(alpha = 0.5, trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.2, fill = "white", alpha = 0.7, outlier.size = 0.5) +
    scale_fill_manual(values = c("All Genes" = "grey70", "Significant Genes" = "steelblue")) +
    labs(
      title = "Distribution Comparison",
      subtitle = sprintf("Significant: n=%d genes | All: n=%s genes",
                         sig_stats$n, format(all_stats$n, big.mark = ",")),
      x = "Gene Set",
      y = sprintf("%s (log2 scale)", SHIFT_LABEL_ABS)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "none"
    )
  
  combined <- p_abs + p_box + plot_layout(widths = c(2, 1))
  
  # Save
  if (project_id == "all") {
    filename <- sprintf("07_signed_diff_comparison_%s_%s_all.png", 
                        gsub(" ", "_", norm_method), comp_type)
  } else {
    filename <- sprintf("07_signed_diff_comparison_%s_%s_%s.png", 
                        project_id, norm_method, comp_type)
  }
  
  ggsave(file.path(output_dir, filename), combined, width = 8, height = 5, dpi = 300)
  
  return(combined)
}

#' Combined signed difference distributions across all projects
#' One plot per normalization method with all 3 comparison types side by side
plot_signed_diff_by_norm <- function(all_data, sig_data, norm_method, output_dir) {
  
  plots <- list()
  
  for (comp_type in COMP_TYPES) {
    all_filtered <- all_data %>%
      filter(norm_method == !!norm_method, comparison == comp_type)
    
    sig_filtered <- sig_data %>%
      filter(norm_method == !!norm_method, comparison == comp_type)
    
    if (nrow(all_filtered) == 0) next
    
    plot_data_all <- all_filtered %>%
      select(signed_difference, cancer_id) %>%
      mutate(type = "All Genes", abs_diff = abs(signed_difference))
    
    plot_data_sig <- sig_filtered %>%
      select(signed_difference, cancer_id) %>%
      mutate(type = "Significant Genes", abs_diff = abs(signed_difference))
    
    plot_data <- bind_rows(plot_data_all, plot_data_sig) %>%
      mutate(abs_diff_log2 = log2(abs_diff + 1e-10))
    
    n_all <- n_distinct(plot_data_all$cancer_id)
    
    p <- ggplot(plot_data, aes(x = abs_diff_log2, fill = type, color = type)) +
      geom_density(alpha = 0.5, adjust = 1.5) +
      scale_fill_manual(values = c("All Genes" = "grey70", "Significant Genes" = "steelblue")) +
      scale_color_manual(values = c("All Genes" = "grey50", "Significant Genes" = "steelblue")) +
      labs(
        title = sprintf("Normalized Signed Difference Distribution - %s", COMP_DISPLAY[comp_type]),
        subtitle = sprintf("%s | %d projects | noise_p_adj < %.2f | log2 scale", 
                          norm_method, n_all, Q_CUTOFF),
        x = sprintf("%s (log2 scale)", SHIFT_LABEL_ABS),
        y = "Density",
        fill = "Gene Set",
        color = "Gene Set"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom"
      )
    
    plots[[comp_type]] <- p
  }
  
  if (length(plots) == 0) return(NULL)
  
  combined <- wrap_plots(plots, ncol = 3) +
    plot_annotation(
      title = sprintf("Normalized Signed Difference Distributions by Comparison Type - %s (log2 scale)", norm_method),
      subtitle = sprintf("All projects combined | Significant genes: noise_p_adj < %.2f", Q_CUTOFF),
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 11)
      )
    )
  
  filename <- sprintf("08_signed_diff_by_norm_%s_all_log2.png", gsub(" ", "_", norm_method))
  ggsave(file.path(output_dir, filename), combined, width = 15, height = 6, dpi = 300)
  
  return(combined)
}

# ==================== SIGNIFICANT GENE PLOTS ====================

#' Plot 1: Distribution of mean absolute shifts for significant genes
plot_significant_distance_distribution <- function(sig_data, project_id, norm_method, comp_type, output_file = NULL) {
  
  comp_data <- sig_data %>% filter(comparison == comp_type, is_outlier_stage)
  if (nrow(comp_data) == 0) return(NULL)
  
  gene_stats <- comp_data %>%
    group_by(gene_id) %>%
    summarise(
      mean_abs_shift = mean(abs(signed_difference), na.rm = TRUE),
      n_stages = n(),
      direction = first(direction)
    ) %>%
    arrange(desc(mean_abs_shift))
  
  n_genes <- n_distinct(comp_data$gene_id)
  
  p <- ggplot(gene_stats, aes(x = mean_abs_shift)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    labs(
      title = sprintf("Significant Genes - %s", project_id),
      subtitle = sprintf("%s [%s]: %d significant genes", 
                        COMP_DISPLAY[comp_type], norm_method, n_genes),
      x = sprintf("Mean %s across stages", SHIFT_LABEL_ABS),
      y = "Number of Genes"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 5, height = 5, dpi = 300)
  }
  
  return(p)
}

#' Plot 2: Top significant gene trajectories (n=5) with consistent colors

plot_significant_top_trajectories <- function(sig_data, project_id, norm_method, comp_type, 
                                              n_top = 5, output_file = NULL, color_manager = NULL) {

  comp_data_all <- sig_data %>% filter(comparison == comp_type)
  comp_data_sig <- comp_data_all %>% filter(is_outlier_stage)

  if (nrow(comp_data_sig) == 0) return(NULL)

  gene_stats <- comp_data_sig %>%
    group_by(gene_id) %>%
    summarise(
      mean_abs_shift = mean(abs(signed_difference), na.rm = TRUE),
      n_outlier_stages = n_distinct(stage),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_shift), desc(n_outlier_stages), gene_id)

  top_genes <- gene_stats %>%
    slice_head(n = n_top) %>%
    pull(gene_id)

  legend_tbl <- gene_stats %>%
    filter(gene_id %in% top_genes) %>%
    arrange(desc(mean_abs_shift), desc(n_outlier_stages), gene_id) %>%
    mutate(legend_label = sprintf("%s (%.2f)", gene_id, mean_abs_shift))

  plot_data <- comp_data_all %>%
    filter(gene_id %in% top_genes, !is.na(signed_difference)) %>%
    left_join(select(gene_stats, gene_id, mean_abs_shift, n_outlier_stages), by = "gene_id") %>%
    left_join(select(legend_tbl, gene_id, legend_label), by = "gene_id") %>%
    mutate(
      stage_num = match(stage, STAGES),
      stage = factor(stage, levels = STAGES),
      outlier_shape = ifelse(is_outlier_stage, "Outlier", "Not outlier"),
      legend_label = factor(legend_label, levels = legend_tbl$legend_label)
    ) %>%
    arrange(legend_label, stage_num)

  if (nrow(plot_data) == 0) return(NULL)

  if (!is.null(color_manager)) {
    gene_colors <- color_manager$get_colors(legend_tbl$gene_id)
    color_values <- setNames(unname(gene_colors[legend_tbl$gene_id]), legend_tbl$legend_label)
  } else {
    fallback_palette <- scales::hue_pal()(nrow(legend_tbl))
    color_values <- setNames(fallback_palette, legend_tbl$legend_label)
  }

  p <- ggplot(plot_data, aes(x = stage_num, y = signed_difference, 
                             color = legend_label, group = gene_id)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_line(linewidth = 1.2, alpha = 0.95, show.legend = TRUE) +
    geom_point(aes(shape = outlier_shape), size = 2.6, stroke = 0.6, show.legend = TRUE) +
    scale_x_continuous(
      breaks = seq_along(STAGES),
      labels = STAGES,
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_shape_manual(
      values = c("Not outlier" = 16, "Outlier" = 17),
      breaks = c("Not outlier", "Outlier")
    ) +
    scale_color_manual(values = color_values, drop = FALSE) +
    guides(
      color = guide_legend(ncol = 1, order = 1,
                           override.aes = list(shape = 16, linewidth = 1.2, alpha = 1)),
      shape = guide_legend(order = 2, nrow = 1, byrow = TRUE)
    ) +
    labs(
      title = sprintf("Top %d Significant Genes - %s", n_top, project_id),
      subtitle = sprintf(
        "%s [%s] - noise_p_adj < %.2f
Legend: Gene ID (mean %s)",
        COMP_DISPLAY[comp_type], norm_method, Q_CUTOFF, SHIFT_LABEL_ABS
      ),
      x = "Stage",
      y = SHIFT_LABEL,
      color = "Gene ID",
      shape = "Stage status"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 7),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.spacing.y = grid::unit(0.03, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 9, color = "grey20"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4)
    )

  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 5.6, height = 4.0, dpi = 300)
  }

  return(p)
}

#' Plot 2b: Capped heatmap (Top 20 significant genes)
plot_significant_heatmap_capped <- function(sig_data, project_id, norm_method, comp_type, 
                                            n_top = 20, output_file = NULL) {
  
  comp_data <- sig_data %>% filter(comparison == comp_type, is_outlier_stage)
  if (nrow(comp_data) == 0) return(NULL)
  
  gene_stats <- comp_data %>%
    group_by(gene_id) %>%
    summarise(
      mean_abs_shift = mean(abs(signed_difference), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_shift)) %>%
    head(n_top)
  
  top_genes <- gene_stats$gene_id
  
  heatmap_data <- comp_data %>%
    filter(gene_id %in% top_genes) %>%
    left_join(gene_stats, by = "gene_id") %>%
    mutate(
      gene_label = sprintf("%s (%.2f)", substr(gene_id, 1, 12), mean_abs_shift),
      stage = factor(stage, levels = STAGES)
    )
  
  gene_order <- gene_stats %>%
    mutate(gene_label = sprintf("%s (%.2f)", substr(gene_id, 1, 12), mean_abs_shift)) %>%
    arrange(desc(mean_abs_shift)) %>%
    pull(gene_label)
  
  heatmap_data <- heatmap_data %>%
    mutate(gene_label = factor(gene_label, levels = gene_order))
  
  p <- ggplot(heatmap_data, aes(x = stage, y = gene_label, fill = signed_difference)) +
    geom_tile() +
    geom_point(aes(shape = direction), size = 1.5, color = "black", stroke = 0.3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = SHIFT_LEGEND_SHORT) +
    scale_shape_manual(values = c("up" = 16, "down" = 17), name = "Direction") +
    labs(
      title = sprintf("Top %d Significant Genes Heatmap (Capped) - %s", n_top, project_id),
      subtitle = sprintf("%s [%s] - noise_p_adj < %.2f\nCell values: %s", 
                        COMP_DISPLAY[comp_type], norm_method, Q_CUTOFF, SHIFT_LABEL),
      x = "Stage",
      y = sprintf("Gene ID (mean %s)", SHIFT_LABEL_ABS)
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom")
  
  if (!is.null(output_file)) {
    height <- max(8, min(20, 8 + n_top * 0.15))
    ggsave(output_file, p, width = 6, height = height, dpi = 300)
  }
  
  return(p)
}

#' Plot 2c: UNCAPPED heatmap - ALL significant genes (no limit)
plot_significant_heatmap_uncapped <- function(sig_data, project_id, norm_method, comp_type, 
                                              output_file = NULL) {
  
  comp_data <- sig_data %>% filter(comparison == comp_type, is_outlier_stage)
  if (nrow(comp_data) == 0) return(NULL)
  
  gene_stats <- comp_data %>%
    group_by(gene_id) %>%
    summarise(
      mean_abs_shift = mean(abs(signed_difference), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_abs_shift))
  
  n_total_genes <- nrow(gene_stats)
  
  heatmap_data <- comp_data %>%
    left_join(gene_stats, by = "gene_id") %>%
    mutate(
      gene_label = sprintf("%s (%.2f)", substr(gene_id, 1, 12), mean_abs_shift),
      stage = factor(stage, levels = STAGES)
    )
  
  gene_order <- gene_stats %>%
    mutate(gene_label = sprintf("%s (%.2f)", substr(gene_id, 1, 12), mean_abs_shift)) %>%
    pull(gene_label)
  
  heatmap_data <- heatmap_data %>%
    mutate(gene_label = factor(gene_label, levels = gene_order))
  
  # Dynamic height based on number of genes
  plot_height <- max(8, min(50, 8 + n_total_genes * 0.08))
  
  p <- ggplot(heatmap_data, aes(x = stage, y = gene_label, fill = signed_difference)) +
    geom_tile() +
    geom_point(aes(shape = direction), size = 1, color = "black", stroke = 0.2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, name = SHIFT_LEGEND_SHORT) +
    scale_shape_manual(values = c("up" = 16, "down" = 17), name = "Direction") +
    labs(
      title = sprintf("ALL Significant Genes Heatmap (Uncapped) - %s", project_id),
      subtitle = sprintf("%s [%s] - %d significant genes (noise_p_adj < %.2f)\nCell values: %s", 
                        COMP_DISPLAY[comp_type], norm_method, n_total_genes, Q_CUTOFF, SHIFT_LABEL),
      x = "Stage",
      y = sprintf("Gene ID (mean %s)", SHIFT_LABEL_ABS)
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom")
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 8, height = plot_height, dpi = 300, limitsize = FALSE)
  }
  
  return(p)
}

#' Plot 3: Stage distribution heatmap for significant genes
plot_significant_stage_heatmap <- function(sig_data, project_id, norm_method, comp_type, output_file = NULL) {
  
  comp_data <- sig_data %>% filter(comparison == comp_type, is_outlier_stage)
  if (nrow(comp_data) == 0) return(NULL)
  
  plot_data <- comp_data %>%
    group_by(stage, direction) %>%
    summarise(
      n_genes = n_distinct(gene_id),
      .groups = "drop"
    ) %>%
    mutate(
      stage = factor(stage, levels = STAGES),
      direction = factor(direction, levels = c("up", "down"))
    )
  
  p <- ggplot(plot_data, aes(x = stage, y = direction, fill = n_genes)) +
    geom_tile() +
    geom_text(aes(label = n_genes), size = 4) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(
      title = sprintf("Significant Genes per Stage - %s", project_id),
      subtitle = sprintf("%s [%s] - noise_p_adj < %.2f\nCell values: %s", 
                        COMP_DISPLAY[comp_type], norm_method, Q_CUTOFF, SHIFT_LABEL),
      x = "Stage",
      y = "Direction",
      fill = "Count"
    ) +
    theme_minimal()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 6, height = 4, dpi = 300)
  }
  
  return(p)
}

#' Plot 4: Direction pie chart for significant genes
plot_significant_direction_pie <- function(sig_data, project_id, norm_method, comp_type, output_file = NULL) {
  
  comp_data <- sig_data %>% filter(comparison == comp_type, is_outlier_stage)
  if (nrow(comp_data) == 0) return(NULL)
  
  plot_data <- comp_data %>%
    distinct(gene_id, direction) %>%
    group_by(direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(percentage = n / sum(n) * 100)
  
  p <- ggplot(plot_data, aes(x = "", y = n, fill = direction)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = paste0(round(percentage, 1), "%\n(n=", n, ")")),
              position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    labs(
      title = sprintf("Direction of Significant Genes - %s", project_id),
      subtitle = sprintf("%s [%s] - noise_p_adj < %.2f\nCell values: %s", 
                        COMP_DISPLAY[comp_type], norm_method, Q_CUTOFF, SHIFT_LABEL),
      fill = "Direction"
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 6, height = 6, dpi = 300)
  }
  
  return(p)
}

#' Plot 5: Signed difference distribution histogram
plot_significant_distribution <- function(sig_data, project_id, norm_method, comp_type, output_file = NULL) {
  
  comp_data <- sig_data %>% filter(comparison == comp_type, is_outlier_stage)
  if (nrow(comp_data) == 0) return(NULL)
  
  p <- ggplot(comp_data, aes(x = signed_difference, fill = direction)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    labs(
      title = sprintf("Signed Difference Distribution - %s", project_id),
      subtitle = sprintf("%s [%s] - noise_p_adj < %.2f\nCell values: %s", 
                        COMP_DISPLAY[comp_type], norm_method, Q_CUTOFF, SHIFT_LABEL),
      x = "Signed Difference",
      y = "Count",
      fill = "Direction"
    ) +
    theme_minimal()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 5, height = 5, dpi = 300)
  }
  
  return(p)
}

#' Plot 6: Comparison barplot across comparison types
plot_significant_comparison_barplot <- function(sig_data, project_id, norm_method, output_file = NULL) {
  
  if (nrow(sig_data) == 0) return(NULL)
  
  plot_data <- sig_data %>%
    filter(is_outlier_stage) %>%
    group_by(comparison, direction) %>%
    summarise(
      n_genes = n_distinct(gene_id),
      .groups = "drop"
    )
  
  p <- ggplot(plot_data, aes(x = comparison, y = n_genes, fill = direction)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    labs(
      title = sprintf("Significant Genes by Comparison - %s", project_id),
      subtitle = sprintf("[%s] - noise_p_adj < %.2f", norm_method, Q_CUTOFF),
      x = "Comparison Type",
      y = "Number of Significant Genes",
      fill = "Direction"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = 5, height = 4, dpi = 300)
  }
  
  return(p)
}

# ==================== MULTI-CANCER SIGNIFICANCE HEATMAP (SIMPLIFIED) ====================

#' Create a multi-cancer significance heatmap using noise-filtered significant genes
#' Shows which genes are significant (over/under expressed) across all projects and stages
#' 
#' @param sig_data Significant genes data frame from load_significant_genes()
#' @param norm_method Normalization method to use
#' @param comp_type Comparison type ("own_healthy", "family_mean", "ortholog_mean")
#' @param output_dir Output directory for saving the plot
#' @param n_top Number of top genes to show (default: 30)
#' @param sort_by Options: "n_stages" (default) or "mean_shift"
#' @return ggplot object
plot_multicancer_significance_heatmap <- function(sig_data, norm_method_local, comp_type, 
                                                   output_dir, n_top = 30, sort_by = "n_stages") {
  
  library(data.table)
  
  # Map comp_type for display
  comp_display <- c(
    "own_healthy" = "Gene vs Own Healthy",
    "family_mean" = "Gene vs Family Mean",
    "ortholog_mean" = "Gene vs Ortholog Mean"
  )
  
  cat(sprintf("\n  Creating heatmap: %s [%s] (sorted by %s)\n", 
              norm_method_local, comp_display[comp_type], 
              ifelse(sort_by == "n_stages", "number of stages", "mean shift")))
  
  # Filter for specific normalization and comparison
  filtered_data <- sig_data %>%
    filter(norm_method == !!norm_method_local, comparison == comp_type, is_outlier_stage)
  
  if (nrow(filtered_data) == 0) {
    cat(sprintf("    No data for %s - %s\n", norm_method_local, comp_type))
    return(NULL)
  }
  
  # Convert to data.table for efficient processing
  dt <- as.data.table(filtered_data)
  
  # Get unique projects
  project_order <- unique(dt$cancer_id)
  
  # Calculate mean absolute shift per gene (across all stages and projects)
  mean_shifts <- dt[, .(mean_abs_shift = mean(abs(signed_difference), na.rm = TRUE)), by = gene_id]
  
  # Calculate total significant stages per gene (count of significant observations)
  gene_stats <- dt[, .(total_sig_stages = .N), by = gene_id]
  
  # Join with mean shifts
  gene_stats <- merge(gene_stats, mean_shifts, by = "gene_id", all.x = TRUE)
  gene_stats[is.na(mean_abs_shift), mean_abs_shift := 0]
  
  # Sort based on chosen method
  gene_stats_sorted <- copy(gene_stats)
  if (sort_by == "n_stages") {
    setorder(gene_stats_sorted, -total_sig_stages, -mean_abs_shift)
  } else {
    setorder(gene_stats_sorted, -mean_abs_shift, -total_sig_stages)
  }
  
  top_genes <- gene_stats_sorted[1:min(n_top, nrow(gene_stats_sorted)), gene_id]
  
  cat(sprintf("    Top genes selected: %d\n", length(top_genes)))
  
  # Prepare heatmap data
  heatmap_data <- dt[gene_id %in% top_genes]
  
  # Add direction (1 = up, -1 = down)
  heatmap_data[, direction := ifelse(signed_difference > 0, 1, -1)]
  
  # Merge with gene stats for labels
  heatmap_data <- merge(heatmap_data, gene_stats, by = "gene_id")
  
  # Create labels
  heatmap_data[, gene_label := sprintf("%s (%.2f)", 
                                       ifelse(nchar(gene_id) > 15, 
                                              paste0(substr(gene_id, 1, 12), "..."), 
                                              gene_id),
                                       mean_abs_shift)]
  heatmap_data[, stage := factor(stage, levels = STAGES)]
  heatmap_data[, col_label := paste(cancer_id, stage, sep = "\n")]
  
  # Order genes based on chosen method
  gene_order_dt <- unique(heatmap_data[, .(gene_id, gene_label, total_sig_stages, mean_abs_shift)])
  if (sort_by == "n_stages") {
    setorder(gene_order_dt, -total_sig_stages, -mean_abs_shift)
  } else {
    setorder(gene_order_dt, -mean_abs_shift, -total_sig_stages)
  }
  heatmap_data[, gene_label := factor(gene_label, levels = gene_order_dt$gene_label)]
  
  # Order columns (projects × stages)
  col_order <- CJ(project = project_order, stage = STAGES)[, col_label := paste(project, stage, sep = "\n")]$col_label
  heatmap_data[, col_label := factor(col_label, levels = col_order)]
  
  # Create subtitle based on sorting method
  if (sort_by == "n_stages") {
    subtitle_text <- sprintf("%s | Top %d genes by number of significant stages (noise_p_adj < %.2f)", 
                            norm_method_local, n_top, Q_CUTOFF)
  } else {
    subtitle_text <- sprintf("%s | Top %d genes by mean %s (noise_p_adj < %.2f)", 
                            norm_method_local, n_top, SHIFT_LABEL_ABS, Q_CUTOFF)
  }
  
  # Create plot
  p <- ggplot(heatmap_data, aes(x = col_label, y = gene_label, fill = factor(direction))) +
    geom_tile(color = "white", size = 0.3) +
    scale_fill_manual(
      values = c("1" = "firebrick", "-1" = "steelblue"),
      name = "Outlier direction",
      labels = c("1" = "Up-regulated", "-1" = "Down-regulated")
    ) +
    labs(
      title = sprintf("Multi-Cancer Significance Heatmap - %s", comp_display[comp_type]),
      subtitle = subtitle_text,
      x = "Cancer Type and Stage",
      y = sprintf("Gene ID (mean %s)", SHIFT_LABEL_ABS)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom",
      panel.grid = element_blank()
    )
  
  # Save
  output_subdir <- file.path(output_dir, "multicancer_heatmaps")
  dir.create(output_subdir, recursive = TRUE, showWarnings = FALSE)
  
  sort_suffix <- ifelse(sort_by == "n_stages", "bynstages", "bymeanshift")
  comp_suffix <- gsub("_", "", comp_type)
  filename <- sprintf("multicancer_heatmap_%s_%s_%s_top%d.png", 
                      norm_method_local, comp_suffix, sort_suffix, n_top)
  height <- max(8, min(20, 8 + n_top * 0.15))
  width <- max(10, length(col_order) * 0.6)
  ggsave(file.path(output_subdir, filename), p, width = width, height = height, dpi = 300)
  cat(sprintf("    Saved: %s\n", filename))
  
  # Save CSV
  csv_data <- dcast(heatmap_data[, .(gene_label, col_label, direction, total_sig_stages, mean_abs_shift)], 
                    gene_label + total_sig_stages + mean_abs_shift ~ col_label, 
                    value.var = "direction", fill = 0)
  
  if (sort_by == "n_stages") {
    setorder(csv_data, -total_sig_stages, -mean_abs_shift)
  } else {
    setorder(csv_data, -mean_abs_shift, -total_sig_stages)
  }
  
  csv_filename <- sprintf("multicancer_heatmap_%s_%s_%s_top%d.csv", 
                          norm_method_local, comp_suffix, sort_suffix, n_top)
  fwrite(csv_data, file.path(output_subdir, csv_filename))
  
  return(p)
}

#' Create multi-cancer significance heatmaps for all normalization methods and comparison types
create_all_multicancer_heatmaps <- function(sig_data, output_dir, n_top = 30) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CREATING MULTI-CANCER SIGNIFICANCE HEATMAPS\n")
  cat("Using noise-filtered significant genes (noise_p_adj < 0.01)\n")
  cat("For each normalization × comparison combination\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  norm_methods <- c("raw", "log", "std_log", "full")
  comp_types <- c("own_healthy", "family_mean", "ortholog_mean")
  sort_options <- c("n_stages", "mean_shift")
  total_plots <- 0
  
  for (norm_method in norm_methods) {
    for (comp_type in comp_types) {
      for (sort_by in sort_options) {
        p <- plot_multicancer_significance_heatmap(
          sig_data = sig_data,
          norm_method_local = norm_method,
          comp_type = comp_type,
          output_dir = output_dir,
          n_top = n_top,
          sort_by = sort_by
        )
        if (!is.null(p)) total_plots <- total_plots + 1
      }
    }
  }
  
  cat(sprintf("\n  Total plots created: %d / %d\n", 
              total_plots, length(norm_methods) * length(comp_types) * length(sort_options)))
}

#' Calculate Jaccard similarity matrix for corrected significant results
#' Uses gene-based significant sets for each normalization method.
#' A gene is included if it is a corrected outlier in at least one stage,
#' so overlap is compared across normalization methods on the union of all stages.
#'
#' @param sig_data Data frame from load_significant_genes()
#' @param project_id Project ID
#' @param comp_type Comparison type ("own_healthy", "family_mean", "ortholog_mean")
#' @param method_order Normalization methods to compare
#' @return Square matrix of Jaccard similarities, or NULL if no rows for this comparison
calculate_corrected_jaccard_matrix <- function(sig_data, project_id, comp_type,
                                               method_order = c("raw", "log", "std_log", "full")) {
  comp_data <- sig_data %>%
    filter(
      cancer_id == project_id,
      comparison == comp_type,
      is_outlier_stage
    )

  if (nrow(comp_data) == 0) {
    return(NULL)
  }

  method_sets <- lapply(method_order, function(method) {
    comp_data %>%
      filter(norm_method == method) %>%
      pull(gene_id) %>%
      unique()
  })
  names(method_sets) <- method_order

  jmat <- matrix(NA_real_,
                 nrow = length(method_order),
                 ncol = length(method_order),
                 dimnames = list(method_order, method_order))

  for (i in seq_along(method_order)) {
    for (j in seq_along(method_order)) {
      set_i <- method_sets[[i]]
      set_j <- method_sets[[j]]

      if (length(set_i) == 0 && length(set_j) == 0) {
        jmat[i, j] <- NA_real_
      } else if (length(union(set_i, set_j)) == 0) {
        jmat[i, j] <- NA_real_
      } else {
        jmat[i, j] <- length(intersect(set_i, set_j)) / length(union(set_i, set_j))
      }
    }
  }

  return(jmat)
}

#' Create combined Jaccard heatmap for corrected significant results for one cancer project
#' Shows overlap between normalization methods across all comparison types.
#' Missing combinations (e.g. no significant results for raw or own_healthy) are shown as NA.
#'
#' @param sig_data Data frame from load_significant_genes()
#' @param project_id Project ID
#' @param output_dir Output directory for saving plots
#' @return Combined patchwork/ggplot object
plot_corrected_jaccard_heatmap <- function(sig_data, project_id, output_dir) {

  comp_types <- c("own_healthy", "family_mean", "ortholog_mean")
  method_order <- c("raw", "log", "std_log", "full")

  comp_labels_full <- c(
    "own_healthy" = "Own Healthy",
    "family_mean" = "Family Mean",
    "ortholog_mean" = "Ortholog Mean"
  )

  norm_labels <- c(
    "raw" = "raw",
    "log" = "log",
    "std_log" = "std_log",
    "full" = "full"
  )

  jaccard_matrices <- list()

  for (comp in comp_types) {
    jmat <- calculate_corrected_jaccard_matrix(sig_data, project_id, comp, method_order = method_order)
    if (!is.null(jmat)) {
      jaccard_matrices[[comp]] <- jmat
    }
  }

  if (length(jaccard_matrices) == 0) {
    cat(sprintf("    No corrected Jaccard matrices calculated for %s\n", project_id))
    return(NULL)
  }

  plots <- list()

  for (comp in names(jaccard_matrices)) {
    jmat <- jaccard_matrices[[comp]]

    plot_data <- expand.grid(
      Method1 = rownames(jmat),
      Method2 = colnames(jmat),
      stringsAsFactors = FALSE
    )
    plot_data$Similarity <- as.vector(jmat)
    plot_data$Method1_label <- factor(norm_labels[plot_data$Method1], levels = norm_labels[method_order])
    plot_data$Method2_label <- factor(norm_labels[plot_data$Method2], levels = norm_labels[method_order])
    plot_data$label <- ifelse(is.na(plot_data$Similarity), "NA", sprintf("%.2f", plot_data$Similarity))

    p <- ggplot(plot_data, aes(x = Method1_label, y = Method2_label, fill = Similarity)) +
      geom_tile(color = "white", size = 0.5) +
      geom_text(aes(label = label), size = 3) +
      scale_fill_gradient2(
        low = "white", mid = "steelblue", high = "darkred",
        midpoint = 0.5, limits = c(0, 1), na.value = "grey90",
        name = "Jaccard\nSimilarity"
      ) +
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

  if (length(plots) == 3) {
    combined <- (plots[[1]] | plots[[2]] | plots[[3]]) +
      plot_annotation(
        title = sprintf("Corrected Jaccard Similarity: Normalization Methods - %s", project_id),
        subtitle = sprintf("Gene-level outlier overlap across all stages (noise_p_adj < %.2f)", Q_CUTOFF),
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
        title = sprintf("Corrected Jaccard Similarity: Normalization Methods - %s", project_id),
        subtitle = sprintf("Gene-level outlier overlap across all stages (noise_p_adj < %.2f)", Q_CUTOFF),
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
        title = sprintf("Corrected Jaccard Similarity: Normalization Methods - %s", project_id),
        subtitle = sprintf("Gene-level outlier overlap across all stages (noise_p_adj < %.2f)", Q_CUTOFF),
        theme = theme(
          plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10)
        )
      )
    width <- 6
    height <- 5
  }

  plot_dir <- file.path(output_dir, "jaccard_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  filename <- sprintf("jaccard_corrected_combined_%s.png", project_id)
  ggsave(file.path(plot_dir, filename), combined, width = width, height = height, dpi = 300)
  cat(sprintf("    Saved corrected Jaccard heatmap: %s\n", filename))

  for (comp in names(jaccard_matrices)) {
    csv_filename <- sprintf("jaccard_corrected_matrix_%s_%s.csv", project_id, comp)
    write.csv(jaccard_matrices[[comp]], file.path(plot_dir, csv_filename), row.names = TRUE)
  }

  return(combined)
}

# ==================== MAIN FUNCTION ====================

create_significant_gene_plots <- function(output_dir) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("SIGNIFICANT GENE TRAJECTORY PLOTS\n")
  cat(sprintf("Filter: noise_p_adj < %.2f\n", Q_CUTOFF))
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Load data
  cat("Loading all genes with adjusted p-values...\n")
  all_data <- load_all_genes()
  
  cat("\nLoading significant genes...\n")
  sig_data <- load_significant_genes()
  
  if (is.null(sig_data) || nrow(sig_data) == 0) {
    stop("No significant genes found!")
  }
  
  cat(sprintf("Total rows (all): %d\n", nrow(all_data)))
  cat(sprintf("Total stage-specific outlier rows: %d\n", sum(sig_data$is_outlier_stage, na.rm = TRUE)))
  cat(sprintf("Unique genes with >=1 outlier stage: %d\n\n", n_distinct(sig_data$gene_id[sig_data$is_outlier_stage])))
  
  plots_dir <- file.path(output_dir, "significant_gene_plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===== COMPARISON PLOTS: SIGNIFICANT vs ALL =====
  cat("Creating significant vs all comparison plots (log2 scale)...\n")
  
  projects <- unique(all_data$cancer_id)
  norm_methods <- unique(all_data$norm_method)
  
  # Individual project plots
  for (project_id in projects) {
    for (norm_method_local in norm_methods) {
      for (comp_type in COMP_TYPES) {
        plot_signed_diff_comparison(
          all_data = all_data,
          sig_data = sig_data,
          project_id = project_id,
          norm_method = norm_method_local,
          comp_type = comp_type,
          output_dir = plots_dir
        )
      }
    }
  }
  
  # Combined plots (all projects)
  for (norm_method_local in norm_methods) {
    for (comp_type in COMP_TYPES) {
      plot_signed_diff_comparison(
        all_data = all_data,
        sig_data = sig_data,
        project_id = "all",
        norm_method = norm_method_local,
        comp_type = comp_type,
        output_dir = plots_dir
      )
    }
    
    plot_signed_diff_by_norm(
      all_data = all_data,
      sig_data = sig_data,
      norm_method = norm_method_local,
      output_dir = plots_dir
    )
  }
  
  # ===== STANDARD SIGNIFICANT GENE PLOTS =====
  cat("\nCreating standard significant gene plots...\n")
  
  all_stats <- data.frame()
  
  for (project_id in projects) {
    color_manager <- create_or_load_color_manager(project_id)

    for (norm_method_local in norm_methods) {
      
      sig_proj_data <- sig_data %>%
        filter(cancer_id == project_id, norm_method == norm_method_local)
      
      if (sum(sig_proj_data$is_outlier_stage, na.rm = TRUE) == 0) next
      
      proj_dir <- file.path(plots_dir, norm_method_local, project_id)
      dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
      
      for (comp_type in COMP_TYPES) {
        n_genes <- sig_proj_data %>% 
          filter(comparison == comp_type, is_outlier_stage) %>% 
          pull(gene_id) %>% 
          unique() %>% 
          length()
        
        if (n_genes == 0) next
        
        # 1. Distance distribution
        plot_significant_distance_distribution(
          sig_proj_data, project_id, norm_method_local, comp_type,
          output_file = file.path(proj_dir, paste0("01_distance_distribution_", comp_type, ".png"))
        )
        
        # 2. Top trajectories (n=5)
        plot_significant_top_trajectories(
          sig_proj_data, project_id, norm_method_local, comp_type, n_top = 5,
          output_file = file.path(proj_dir, paste0("02_top_trajectories_", comp_type, ".png")),
          color_manager = color_manager
        )
        
        # 2b. Capped heatmap (Top 20)
        if (n_genes >= 10) {
          plot_significant_heatmap_capped(
            sig_proj_data, project_id, norm_method_local, comp_type, n_top = 20,
            output_file = file.path(proj_dir, paste0("02b_heatmap_capped_20_", comp_type, ".png"))
          )
        }
        
        # 2c. Uncapped heatmap (ALL genes)
        plot_significant_heatmap_uncapped(
          sig_proj_data, project_id, norm_method_local, comp_type,
          output_file = file.path(proj_dir, paste0("02c_heatmap_uncapped_all_", comp_type, ".png"))
        )
        
        # 3. Stage heatmap
        plot_significant_stage_heatmap(
          sig_proj_data, project_id, norm_method_local, comp_type,
          output_file = file.path(proj_dir, paste0("03_stage_heatmap_", comp_type, ".png"))
        )
        
        # 4. Direction pie
        plot_significant_direction_pie(
          sig_proj_data, project_id, norm_method_local, comp_type,
          output_file = file.path(proj_dir, paste0("04_direction_pie_", comp_type, ".png"))
        )
        
        # 5. Distribution
        plot_significant_distribution(
          sig_proj_data, project_id, norm_method_local, comp_type,
          output_file = file.path(proj_dir, paste0("05_distribution_", comp_type, ".png"))
        )
      }
      
      # 6. Comparison barplot
      plot_significant_comparison_barplot(
        sig_proj_data, project_id, norm_method_local,
        output_file = file.path(proj_dir, "06_comparison_barplot.png")
      )
      
      # Collect statistics
      stats <- sig_proj_data %>%
        filter(is_outlier_stage) %>%
        group_by(comparison, direction) %>%
        summarise(
          n_genes = n_distinct(gene_id),
          .groups = "drop"
        ) %>%
        mutate(
          project = project_id,
          norm = norm_method_local
        )
      all_stats <- bind_rows(all_stats, stats)
    }

    plot_corrected_jaccard_heatmap(
      sig_data = sig_data,
      project_id = project_id,
      output_dir = plots_dir
    )

    if (!is.null(color_manager) && is.function(color_manager$save_map)) {
      color_manager$save_map(file.path(dirname(GENE_OUTPUT_DIR), project_id))
    }
  }
  
  # Save statistics
  write.csv(all_stats, file.path(plots_dir, "significant_genes_summary.csv"), row.names = FALSE)
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SIGNIFICANT GENE PLOTS COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("Results in: %s\n", plots_dir))
}

# ==================== MAIN EXECUTION ====================

if (sys.nframe() == 0) {
  output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create significant gene plots (individual projects)
  create_significant_gene_plots(output_dir)
  
  # Load significant genes for multi-cancer heatmaps
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("LOADING SIGNIFICANT GENES FOR MULTI-CANCER HEATMAPS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  sig_data <- load_significant_genes()
  
  # Create multi-cancer heatmaps
  if (!is.null(sig_data) && nrow(sig_data) > 0) {
    create_all_multicancer_heatmaps(
      sig_data = sig_data,
      output_dir = output_dir,
      n_top = 30
    )
  } else {
    cat("No significant genes found, skipping multi-cancer heatmaps\n")
  }
}