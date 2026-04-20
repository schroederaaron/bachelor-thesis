#!/usr/bin/env Rscript
# compare_tensoromics_with_de.R
# Compares TensorOmics outliers with differentially expressed genes from edgeR/limma
# - Central all_genes_pvalues.rds file
# - Noise threshold flexible (0.01, 0.05)
# - Top N comparisons (capped and uncapped versions)
# - Strict stage separation
# - Heatmap export with NA values for missing data
# - Distribution plots of signed differences

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)
library(reshape2)
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
  "std_log" = "log(mean(gene-wise scaled TPM))",
  "full" = "log(mean(quantile(gene-wise scaled TPM)))"
)

# All three comparison types
COMP_TYPES <- c("own_healthy", "family_mean", "ortholog_mean")
COMP_DISPLAY <- c(
  "own_healthy" = "Gene vs own healthy",
  "family_mean" = "Gene vs family mean",
  "ortholog_mean" = "Gene vs ortholog mean"
)

# Parameter combinations to test
PARAM_COMBINATIONS <- list(
  list(noise = 0.01, name = "n01"),
  list(noise = 0.05, name = "n05")
)

# Top N for comparisons (for capped version)
TOP_N_VALUES <- c(500, 1000)

# ==================== DATA LOADING ====================

#' Load TensorOmics data from central all_genes_pvalues file
load_tensoromics_data <- function(project_id) {
  sig_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")
  data_file <- file.path(sig_dir, "all_genes_pvalues.rds")
  
  if (!file.exists(data_file)) {
    return(NULL)
  }
  
  all_data <- readRDS(data_file)
  project_data <- all_data %>% filter(cancer_id == project_id)
  
  if (nrow(project_data) == 0) {
    return(NULL)
  }
  
  return(project_data)
}

#' Extract significant genes with noise threshold
extract_significant_genes <- function(data, norm_method, comp_type, stage, noise_alpha = 0.01) {
  
  # Filter by normalization, comparison type, and stage
  filtered <- data %>%
    filter(
      norm_method == !!norm_method,
      comparison == !!comp_type,
      stage == !!stage
    )
  
  if (nrow(filtered) == 0) {
    return(NULL)
  }
  
  # Apply noise correction and filter
  significant <- filtered %>%
    filter(noise_p_value_adj < noise_alpha)
  
  return(significant)
}

#' Load edgeR results (stage-specific)
load_edger_results <- function(project_id, stage) {
  file_path <- file.path(DE_OUTPUT_DIR, project_id, stage, "edgeR_results.csv")
  
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  edger_data <- read.csv(file_path)
  
  # Significant genes (FDR < 0.01)
  significant <- edger_data %>%
    filter(FDR < 0.01) %>%
    arrange(PValue)
  
  return(significant)
}

#' Load limma/voom results (stage-specific)
load_limma_results <- function(project_id, stage) {
  file_path <- file.path(DE_OUTPUT_DIR, project_id, stage, "voom_results.csv")
  
  if (!file.exists(file_path)) {
    return(NULL)
  }
  
  limma_data <- read.csv(file_path)
  
  # Significant genes (adj.P.Val < 0.01)
  significant <- limma_data %>%
    filter(adj.P.Val < 0.01) %>%
    arrange(P.Value)
  
  return(significant)
}

# ==================== HELPER FUNCTIONS ====================

#' Calculate percentage of TO genes found in DE
percentage_in_de <- function(to_genes, de_genes) {
  if (length(to_genes) == 0) return(NA_real_)
  in_de <- intersect(to_genes, de_genes)
  return(length(in_de) / length(to_genes) * 100)
}

# ==================== DISTRIBUTION PLOTS ====================

#' Create distribution plots of signed differences
create_distribution_plots <- function(to_all_data, to_significant, de_edger, de_limma,
                                       project_id, norm_method, stage, comp_type,
                                       param_name, output_dir) {
  
  # Extract signed differences
  all_distances <- abs(to_all_data$signed_difference)
  to_sig_distances <- if (!is.null(to_significant) && nrow(to_significant) > 0) {
    abs(to_significant$signed_difference)
  } else {
    numeric(0)
  }
  
  # DE gene IDs
  de_edger_ids <- if (!is.null(de_edger)) de_edger$gene_id else character(0)
  de_limma_ids <- if (!is.null(de_limma)) de_limma$gene_id else character(0)
  
  # TO genes in DE methods
  to_all_ids <- to_all_data$gene_id
  to_sig_ids <- if (!is.null(to_significant)) to_significant$gene_id else character(0)
  
  # Distances for genes in DE methods
  edger_distances <- if (length(de_edger_ids) > 0) {
    to_all_data %>% 
      filter(gene_id %in% de_edger_ids) %>%
      pull(signed_difference) %>%
      abs()
  } else {
    numeric(0)
  }
  
  limma_distances <- if (length(de_limma_ids) > 0) {
    to_all_data %>% 
      filter(gene_id %in% de_limma_ids) %>%
      pull(signed_difference) %>%
      abs()
  } else {
    numeric(0)
  }
  
  # Overlaps: TO significant & in edgeR
  to_sig_in_edger_distances <- if (length(to_sig_ids) > 0 && length(de_edger_ids) > 0) {
    intersect_ids <- intersect(to_sig_ids, de_edger_ids)
    if (length(intersect_ids) > 0) {
      to_all_data %>%
        filter(gene_id %in% intersect_ids) %>%
        pull(signed_difference) %>%
        abs()
    } else {
      numeric(0)
    }
  } else {
    numeric(0)
  }
  
  # Overlaps: TO significant & in limma
  to_sig_in_limma_distances <- if (length(to_sig_ids) > 0 && length(de_limma_ids) > 0) {
    intersect_ids <- intersect(to_sig_ids, de_limma_ids)
    if (length(intersect_ids) > 0) {
      to_all_data %>%
        filter(gene_id %in% intersect_ids) %>%
        pull(signed_difference) %>%
        abs()
    } else {
      numeric(0)
    }
  } else {
    numeric(0)
  }
  
  # Common x-axis limits for all plots
  all_values <- c(all_distances, to_sig_distances, edger_distances, limma_distances,
                  to_sig_in_edger_distances, to_sig_in_limma_distances)
  x_limits <- if (length(all_values) > 0) {
    quantile(all_values, probs = c(0, 0.99), na.rm = TRUE)
  } else {
    c(0, 1)
  }
  
  # Create plots
  # Plot 1: TO All + TO Significant + edgeR
  p1 <- ggplot() +
    geom_density(data = data.frame(distance = all_distances), 
                 aes(x = distance, fill = "TO All"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = to_sig_distances), 
                 aes(x = distance, fill = "TO Significant"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = edger_distances), 
                 aes(x = distance, fill = "DEG edgeR"), alpha = 0.5, adjust = 1.5) +
    scale_fill_manual(values = c("TO All" = "gray70", 
                                  "TO Significant" = "steelblue", 
                                  "DEG edgeR" = "darkorange"),
                      name = NULL) +
    labs(title = "1. TO All + TO Significant + DEG edgeR",
         x = "|signed difference|", y = "Density") +
    xlim(x_limits) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(size = 10, face = "bold"))
  
  # Plot 2: TO All + TO Significant + limma
  p2 <- ggplot() +
    geom_density(data = data.frame(distance = all_distances), 
                 aes(x = distance, fill = "TO All"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = to_sig_distances), 
                 aes(x = distance, fill = "TO Significant"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = limma_distances), 
                 aes(x = distance, fill = "DEG limma"), alpha = 0.5, adjust = 1.5) +
    scale_fill_manual(values = c("TO All" = "gray70", 
                                  "TO Significant" = "steelblue", 
                                  "DEG limma" = "forestgreen"),
                      name = NULL) +
    labs(title = "2. TO All + TO Significant + DEG limma",
         x = "|signed difference|", y = "Density") +
    xlim(x_limits) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(size = 10, face = "bold"))
  
  # Plot 3: TO All + TO Significant only
  p3 <- ggplot() +
    geom_density(data = data.frame(distance = all_distances), 
                 aes(x = distance, fill = "TO All"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = to_sig_distances), 
                 aes(x = distance, fill = "TO Significant"), alpha = 0.5, adjust = 1.5) +
    scale_fill_manual(values = c("TO All" = "gray70", "TO Significant" = "steelblue"),
                      name = NULL) +
    labs(title = "3. TO All + TO Significant",
         x = "|signed difference|", y = "Density") +
    xlim(x_limits) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(size = 10, face = "bold"))
  
  # Plot 4: TO All + edgeR only
  p4 <- ggplot() +
    geom_density(data = data.frame(distance = all_distances), 
                 aes(x = distance, fill = "TO All"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = edger_distances), 
                 aes(x = distance, fill = "DEG edgeR"), alpha = 0.5, adjust = 1.5) +
    scale_fill_manual(values = c("TO All" = "gray70", "DEG edgeR" = "darkorange"),
                      name = NULL) +
    labs(title = "4. TO All + DEG edgeR only",
         x = "|signed difference|", y = "Density") +
    xlim(x_limits) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(size = 10, face = "bold"))
  
  # Plot 5: TO All + limma only
  p5 <- ggplot() +
    geom_density(data = data.frame(distance = all_distances), 
                 aes(x = distance, fill = "TO All"), alpha = 0.5, adjust = 1.5) +
    geom_density(data = data.frame(distance = limma_distances), 
                 aes(x = distance, fill = "DEG limma"), alpha = 0.5, adjust = 1.5) +
    scale_fill_manual(values = c("TO All" = "gray70", "DEG limma" = "forestgreen"),
                      name = NULL) +
    labs(title = "5. TO All + DEG limma only",
         x = "|signed difference|", y = "Density") +
    xlim(x_limits) +
    theme_minimal() +
    theme(legend.position = "right", plot.title = element_text(size = 10, face = "bold"))
  
  # Combine into 2x3 grid
  combined_plot <- (p1 + p2) / (p3 + p4) / (p5 + plot_spacer()) +
    plot_annotation(
      title = sprintf("%s | %s | %s | %s | %s", 
                      project_id, NORM_DISPLAY[norm_method], stage, 
                      COMP_DISPLAY[comp_type], param_name),
      subtitle = "Distribution of |signed differences|",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                    plot.subtitle = element_text(hjust = 0.5, size = 10))
    )
  
  # Save
  plot_filename <- sprintf("dist_%s_%s_%s_%s_%s.png",
                           gsub(" ", "_", project_id),
                           norm_method,
                           gsub(" ", "_", stage),
                           comp_type,
                           param_name)
  
  ggsave(file.path(output_dir, "distribution_plots", plot_filename),
         combined_plot, width = 14, height = 12, dpi = 300)
  
  return(combined_plot)
}

# ==================== HEATMAP LABEL HELPERS ====================

make_heatmap_title <- function(project_id, noise_alpha, top_n = NULL, uncapped = FALSE, metric = "percentage") {
  scope_label <- if (uncapped) {
    "All significant TensorOmics genes"
  } else {
    sprintf("Top %d significant TensorOmics genes", top_n)
  }

  metric_label <- if (metric == "count") {
    "Counts by stage"
  } else {
    "Overlap with DE methods"
  }

  sprintf("%s | Noise threshold %.2f | %s | %s",
          project_id, noise_alpha, scope_label, metric_label)
}

make_heatmap_subtitle <- function(metric = "percentage") {
  if (metric == "count") {
    "Cell values show the number of significant TensorOmics genes per stage. Gray = no data"
  } else {
    "Cell values show the percentage of significant TensorOmics genes also identified by edgeR or limma. Gray = no data"
  }
}

# ==================== HEATMAP FUNCTIONS ====================

#' Create heatmap for a project and parameter combination
create_project_heatmap <- function(stats_data, project_id, params, output_dir, top_n = 500) {
  
  noise_alpha <- params$noise
  param_name <- params$name
  
  # Define all possible rows (comparison type + normalization)
  row_names <- c(
    "own_vs_healthy_raw", "own_vs_healthy_log", "own_vs_healthy_std_log", "own_vs_healthy_full",
    "family_mean_raw", "family_mean_log", "family_mean_std_log", "family_mean_full",
    "ortholog_mean_raw", "ortholog_mean_log", "ortholog_mean_std_log", "ortholog_mean_full"
  )
  
  # Define all possible columns (stage + DE method)
  col_names <- c(
    "edgeR - Stage I", "limma - Stage I",
    "edgeR - Stage II", "limma - Stage II",
    "edgeR - Stage III", "limma - Stage III",
    "edgeR - Stage IV", "limma - Stage IV"
  )
  
  # Create complete grid with NA values
  complete_grid <- expand.grid(
    row = row_names,
    col = col_names,
    stringsAsFactors = FALSE
  )
  complete_grid$percentage <- NA_real_
  complete_grid$n_genes <- NA_integer_
  
  # Fill in existing data
  if (!is.null(stats_data[[as.character(noise_alpha)]])) {
    
    for (comp_type in names(stats_data[[as.character(noise_alpha)]])) {
      for (norm in names(stats_data[[as.character(noise_alpha)]][[comp_type]])) {
        
        stats <- stats_data[[as.character(noise_alpha)]][[comp_type]][[norm]]
        
        # Create row name
        if (comp_type == "own_healthy") {
          row_base <- "own_vs_healthy"
        } else {
          row_base <- comp_type
        }
        row_name <- paste0(row_base, "_", norm)
        
        # For each stage
        for (stage in names(stats)) {
          for (de_method in c("edgeR", "limma")) {
            pct_col <- ifelse(de_method == "edgeR", "pct_in_edger", "pct_in_limma")
            pct_value <- stats[[stage]][[pct_col]]
            n_to <- stats[[stage]][["n_to"]]
            
            if (!is.null(pct_value) && !is.na(pct_value)) {
              col_name <- paste0(de_method, " - ", stage)
              
              # Find the correct row in complete_grid and update
              idx <- which(complete_grid$row == row_name & complete_grid$col == col_name)
              if (length(idx) > 0) {
                complete_grid$percentage[idx] <- pct_value
                complete_grid$n_genes[idx] <- n_to
              }
            }
          }
        }
      }
    }
  }
  
  # Convert to factor with defined order
  complete_grid$row <- factor(complete_grid$row, levels = row_names)
  complete_grid$col <- factor(complete_grid$col, levels = col_names)
  
  # Remove rows that are completely NA
  rows_with_data <- complete_grid %>%
    group_by(row) %>%
    summarise(has_data = any(!is.na(percentage))) %>%
    filter(has_data) %>%
    pull(row)
  
  plot_data <- complete_grid %>% filter(row %in% rows_with_data)
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # ===== PERCENTAGE HEATMAP =====
  p_pct <- ggplot(plot_data, aes(x = col, y = row, fill = percentage)) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(percentage), "NA", sprintf("%.0f%%", percentage))), 
              size = 3, color = ifelse(is.na(plot_data$percentage), "gray50", "black")) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100), 
                        na.value = "lightgray", name = "% in DE") +
    labs(
      title = make_heatmap_title(project_id, noise_alpha, top_n = top_n, uncapped = FALSE, metric = "percentage"),
      subtitle = make_heatmap_subtitle("percentage"),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )
  
  ggsave(
    file.path(output_dir, sprintf("heatmap_%s_%s_top%d.png", 
                                   project_id, param_name, top_n)),
    p_pct, width = 12, height = 10, dpi = 300
  )
  
  # ===== GENE COUNT HEATMAP =====
  count_plot_data <- plot_data %>%
    mutate(stage_only = gsub("^(edgeR|limma) - ", "", as.character(col))) %>%
    group_by(row, stage_only) %>%
    summarise(
      n_genes = dplyr::first(n_genes[!is.na(n_genes)]),
      .groups = "drop"
    )

  count_col_names <- c("Stage I", "Stage II", "Stage III", "Stage IV")
  count_plot_data$stage_only <- factor(count_plot_data$stage_only, levels = count_col_names)

  p_n <- ggplot(count_plot_data, aes(x = stage_only, y = row, fill = n_genes)) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(n_genes), "NA", n_genes)), size = 3,
              color = ifelse(is.na(count_plot_data$n_genes), "gray50", "black")) +
    scale_fill_gradient(low = "white", high = "darkgreen", 
                        na.value = "lightgray", name = "TO gene count") +
    labs(
      title = make_heatmap_title(project_id, noise_alpha, top_n = top_n, uncapped = FALSE, metric = "count"),
      subtitle = make_heatmap_subtitle("count"),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )
  
  ggsave(
    file.path(output_dir, sprintf("heatmap_%s_%s_top%d_counts.png", 
                                   project_id, param_name, top_n)),
    p_n, width = 12, height = 10, dpi = 300
  )
  
  # ===== SAVE AS CSV (with NA) =====
  csv_data <- plot_data %>%
    select(row, col, percentage, n_genes) %>%
    pivot_wider(
      id_cols = row,
      names_from = col,
      values_from = c(percentage, n_genes),
      names_sep = "_"
    )
  
  write.csv(csv_data, 
            file.path(output_dir, sprintf("table_%s_%s_top%d.csv", 
                                           project_id, param_name, top_n)),
            row.names = FALSE)
  
  return(list(
    percentage_plot = p_pct,
    count_plot = p_n,
    data = plot_data
  ))
}

#' Create heatmaps for all projects (capped and uncapped versions)
create_all_heatmaps <- function(all_stats_list, output_dir, top_n_values = c(500, 1000)) {
  
  cat("\nCreating heatmaps...\n")
  
  all_heatmaps <- list()
  
  for (project_id in names(all_stats_list)) {
    cat(sprintf("  %s\n", project_id))
    
    project_heatmaps <- list()
    
    for (params in PARAM_COMBINATIONS) {
      
      # Extract data for this project and parameters
      noise_alpha <- as.character(params$noise)
      
      if (!is.null(all_stats_list[[project_id]][[noise_alpha]])) {
        
        for (top_n in top_n_values) {
          
          # Create data structure for this top N
          topN_data <- list()
          topN_data[[noise_alpha]] <- list()
          
          for (comp_type in names(all_stats_list[[project_id]][[noise_alpha]])) {
            topN_data[[noise_alpha]][[comp_type]] <- list()
            
            for (norm in names(all_stats_list[[project_id]][[noise_alpha]][[comp_type]])) {
              if (!is.null(all_stats_list[[project_id]][[noise_alpha]][[comp_type]][[norm]]$stage_topN[[as.character(top_n)]])) {
                topN_data[[noise_alpha]][[comp_type]][[norm]] <- 
                  all_stats_list[[project_id]][[noise_alpha]][[comp_type]][[norm]]$stage_topN[[as.character(top_n)]]
              }
            }
          }
          
          # Create heatmap for these parameters
          hm <- create_project_heatmap(
            stats_data = topN_data,
            project_id = project_id,
            params = params,
            output_dir = output_dir,
            top_n = top_n
          )
          
          if (!is.null(hm)) {
            key <- sprintf("%s_top%d", params$name, top_n)
            project_heatmaps[[key]] <- hm
          }
        }
      }
    }
    
    all_heatmaps[[project_id]] <- project_heatmaps
  }
  
  return(all_heatmaps)
}

# ==================== UNCAPPED HEATMAP FUNCTIONS ====================

#' Create uncapped heatmap for a project (all genes, no top N limit)
create_project_heatmap_uncapped <- function(stats_data, project_id, params, output_dir) {
  
  noise_alpha <- params$noise
  param_name <- params$name
  
  # Define all possible rows (comparison type + normalization)
  row_names <- c(
    "own_vs_healthy_raw", "own_vs_healthy_log", "own_vs_healthy_std_log", "own_vs_healthy_full",
    "family_mean_raw", "family_mean_log", "family_mean_std_log", "family_mean_full",
    "ortholog_mean_raw", "ortholog_mean_log", "ortholog_mean_std_log", "ortholog_mean_full"
  )
  
  # Define all possible columns (stage + DE method)
  col_names <- c(
    "edgeR - Stage I", "limma - Stage I",
    "edgeR - Stage II", "limma - Stage II",
    "edgeR - Stage III", "limma - Stage III",
    "edgeR - Stage IV", "limma - Stage IV"
  )
  
  # Create complete grid with NA values
  complete_grid <- expand.grid(
    row = row_names,
    col = col_names,
    stringsAsFactors = FALSE
  )
  complete_grid$percentage <- NA_real_
  complete_grid$n_genes <- NA_integer_
  
  # Fill in existing data
  if (!is.null(stats_data[[as.character(noise_alpha)]])) {
    
    for (comp_type in names(stats_data[[as.character(noise_alpha)]])) {
      for (norm in names(stats_data[[as.character(noise_alpha)]][[comp_type]])) {
        
        stats <- stats_data[[as.character(noise_alpha)]][[comp_type]][[norm]]
        
        # Create row name
        if (comp_type == "own_healthy") {
          row_base <- "own_vs_healthy"
        } else {
          row_base <- comp_type
        }
        row_name <- paste0(row_base, "_", norm)
        
        # For each stage
        for (stage in names(stats$stage_all)) {
          for (de_method in c("edgeR", "limma")) {
            pct_col <- ifelse(de_method == "edgeR", "pct_in_edger", "pct_in_limma")
            pct_value <- stats$stage_all[[stage]][[pct_col]]
            n_to <- stats$stage_all[[stage]][["n_to"]]
            
            if (!is.null(pct_value) && !is.na(pct_value)) {
              col_name <- paste0(de_method, " - ", stage)
              
              idx <- which(complete_grid$row == row_name & complete_grid$col == col_name)
              if (length(idx) > 0) {
                complete_grid$percentage[idx] <- pct_value
                complete_grid$n_genes[idx] <- n_to
              }
            }
          }
        }
      }
    }
  }
  
  # Convert to factor with defined order
  complete_grid$row <- factor(complete_grid$row, levels = row_names)
  complete_grid$col <- factor(complete_grid$col, levels = col_names)
  
  # Remove rows that are completely NA
  rows_with_data <- complete_grid %>%
    group_by(row) %>%
    summarise(has_data = any(!is.na(percentage))) %>%
    filter(has_data) %>%
    pull(row)
  
  plot_data <- complete_grid %>% filter(row %in% rows_with_data)
  
  if (nrow(plot_data) == 0) {
    return(NULL)
  }
  
  # ===== PERCENTAGE HEATMAP (UNCAPPED) =====
  p_pct <- ggplot(plot_data, aes(x = col, y = row, fill = percentage)) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(percentage), "NA", sprintf("%.0f%%", percentage))), 
              size = 3, color = ifelse(is.na(plot_data$percentage), "gray50", "black")) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 100), 
                        na.value = "lightgray", name = "% in DE") +
    labs(
      title = make_heatmap_title(project_id, noise_alpha, uncapped = TRUE, metric = "percentage"),
      subtitle = make_heatmap_subtitle("percentage"),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )
  
  ggsave(
    file.path(output_dir, sprintf("heatmap_%s_%s_uncapped_all.png", 
                                   project_id, param_name)),
    p_pct, width = 12, height = 10, dpi = 300
  )
  
  # ===== GENE COUNT HEATMAP (UNCAPPED) =====
  count_plot_data <- plot_data %>%
    mutate(stage_only = gsub("^(edgeR|limma) - ", "", as.character(col))) %>%
    group_by(row, stage_only) %>%
    summarise(
      n_genes = dplyr::first(n_genes[!is.na(n_genes)]),
      .groups = "drop"
    )

  count_col_names <- c("Stage I", "Stage II", "Stage III", "Stage IV")
  count_plot_data$stage_only <- factor(count_plot_data$stage_only, levels = count_col_names)

  p_n <- ggplot(count_plot_data, aes(x = stage_only, y = row, fill = n_genes)) +
    geom_tile() +
    geom_text(aes(label = ifelse(is.na(n_genes), "NA", n_genes)), size = 3,
              color = ifelse(is.na(count_plot_data$n_genes), "gray50", "black")) +
    scale_fill_gradient(low = "white", high = "darkgreen", 
                        na.value = "lightgray", name = "TO gene count") +
    labs(
      title = make_heatmap_title(project_id, noise_alpha, uncapped = TRUE, metric = "count"),
      subtitle = make_heatmap_subtitle("count"),
      x = "",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9),
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "right"
    )
  
  ggsave(
    file.path(output_dir, sprintf("heatmap_%s_%s_uncapped_all_counts.png", 
                                   project_id, param_name)),
    p_n, width = 12, height = 10, dpi = 300
  )
  
  # ===== SAVE AS CSV (with NA) =====
  csv_data <- plot_data %>%
    select(row, col, percentage, n_genes) %>%
    pivot_wider(
      id_cols = row,
      names_from = col,
      values_from = c(percentage, n_genes),
      names_sep = "_"
    )
  
  write.csv(csv_data, 
            file.path(output_dir, sprintf("table_%s_%s_uncapped_all.csv", 
                                           project_id, param_name)),
            row.names = FALSE)
  
  return(list(
    percentage_plot = p_pct,
    count_plot = p_n,
    data = plot_data
  ))
}

#' Create uncapped heatmaps for all projects
create_uncapped_heatmaps <- function(all_stats_list, output_dir) {
  
  cat("\nCreating uncapped heatmaps (all genes)...\n")
  
  all_heatmaps <- list()
  
  for (project_id in names(all_stats_list)) {
    cat(sprintf("  %s\n", project_id))
    
    project_heatmaps <- list()
    
    for (params in PARAM_COMBINATIONS) {
      
      hm <- create_project_heatmap_uncapped(
        stats_data = all_stats_list[[project_id]],
        project_id = project_id,
        params = params,
        output_dir = output_dir
      )
      
      if (!is.null(hm)) {
        project_heatmaps[[params$name]] <- hm
      }
    }
    
    all_heatmaps[[project_id]] <- project_heatmaps
  }
  
  return(all_heatmaps)
}

# ==================== MAIN FUNCTION ====================

compare_tensoromics_with_de <- function(project_id, output_dir, all_stats_list = NULL) {
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("COMPARISON: %s\n", project_id))
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Load TensorOmics data (once for all normalizations)
  to_all_data <- load_tensoromics_data(project_id)
  
  if (is.null(to_all_data)) {
    cat("  No TensorOmics data. Skipping project.\n")
    return(all_stats_list)
  }
  
  proj_out_dir <- file.path(output_dir, project_id)
  dir.create(proj_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Directory for distribution plots
  dist_plot_dir <- file.path(output_dir, "distribution_plots")
  dir.create(dist_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize stats list for this project
  if (is.null(all_stats_list)) {
    all_stats_list <- list()
  }
  if (is.null(all_stats_list[[project_id]])) {
    all_stats_list[[project_id]] <- list()
  }
  
  # For each comparison type
  for (comp_type in COMP_TYPES) {
    comp_display <- COMP_DISPLAY[comp_type]
    cat(sprintf("\n>>> Comparison type: %s\n", comp_display))
    
    # For each stage
    for (stage in STAGES) {
      cat(sprintf("\n  Stage: %s\n", stage))
      
      # Load DE results for THIS stage
      edger_all <- load_edger_results(project_id, stage)
      limma_all <- load_limma_results(project_id, stage)
      
      if (is.null(edger_all) || is.null(limma_all)) {
        cat("    DE data incomplete. Skipping.\n")
        next
      }
      
      cat(sprintf("    edgeR: %d significant genes\n", nrow(edger_all)))
      cat(sprintf("    limma: %d significant genes\n", nrow(limma_all)))
      
      # For each parameter combination
      for (params in PARAM_COMBINATIONS) {
        noise_alpha <- params$noise
        param_name <- params$name
        
        cat(sprintf("\n    --- Parameters: %s (noise=%.2f) ---\n", 
                    param_name, noise_alpha))
        
        # Collect TensorOmics data for all normalizations
        to_data_list <- list()
        
        for (norm in NORM_METHODS) {
          to_data <- extract_significant_genes(
            data = to_all_data,
            norm_method = norm,
            comp_type = comp_type,
            stage = stage,
            noise_alpha = noise_alpha
          )
          
          if (!is.null(to_data) && nrow(to_data) > 0) {
            to_data_list[[NORM_DISPLAY[norm]]] <- to_data
            cat(sprintf("      %s: %d significant genes\n", 
                        NORM_DISPLAY[norm], nrow(to_data)))
          }
        }
        
        if (length(to_data_list) == 0) {
          cat("      No TensorOmics data\n")
          next
        }
        
        # ===== ALL GENES STATISTICS =====
        cat("\n      --- All genes statistics ---\n")
        
        # Initialize stats for this combination
        if (is.null(all_stats_list[[project_id]][[as.character(noise_alpha)]])) {
          all_stats_list[[project_id]][[as.character(noise_alpha)]] <- list()
        }
        if (is.null(all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]])) {
          all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]] <- list()
        }
        
        for (i in seq_along(to_data_list)) {
          norm_name <- names(to_data_list)[i]
          to_genes <- to_data_list[[i]]$gene_id
          norm_key <- names(NORM_DISPLAY)[which(NORM_DISPLAY == norm_name)]
          
          # Calculate percentages
          pct_edger <- percentage_in_de(to_genes, edger_all$gene_id)
          pct_limma <- percentage_in_de(to_genes, limma_all$gene_id)
          
          cat(sprintf("        %s: %d TO genes, %.1f%% in edgeR, %.1f%% in limma\n", 
                      norm_name, length(to_genes), pct_edger, pct_limma))
          
          # Store for heatmap (all genes)
          if (is.null(all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]])) {
            all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]] <- list()
          }
          
          if (is.null(all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_all)) {
            all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_all <- list()
          }
          
          all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_all[[stage]] <- list(
            pct_in_edger = pct_edger,
            pct_in_limma = pct_limma,
            n_to = length(to_genes)
          )
        }
        
        # ===== TOP N STATISTICS (for multiple N values) =====
        for (top_n in TOP_N_VALUES) {
          cat(sprintf("\n      --- Top %d statistics ---\n", top_n))
          
          for (i in seq_along(to_data_list)) {
            norm_name <- names(to_data_list)[i]
            to_data <- to_data_list[[i]]
            norm_key <- names(NORM_DISPLAY)[which(NORM_DISPLAY == norm_name)]
            
            # Top N genes by distance_p_value
            n_to_top <- min(top_n, nrow(to_data))
            to_genes_top <- to_data %>%
              head(n_to_top) %>%
              pull(gene_id)
            
            # Percentages for Top N
            pct_edger_top <- percentage_in_de(to_genes_top, edger_all$gene_id)
            pct_limma_top <- percentage_in_de(to_genes_top, limma_all$gene_id)
            
            cat(sprintf("        %s Top %d: %d genes, %.1f%% in edgeR, %.1f%% in limma\n", 
                        norm_name, top_n, length(to_genes_top), pct_edger_top, pct_limma_top))
            
            # Store for heatmap (Top N)
            if (is.null(all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_topN)) {
              all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_topN <- list()
            }
            if (is.null(all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_topN[[as.character(top_n)]])) {
              all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_topN[[as.character(top_n)]] <- list()
            }
            
            all_stats_list[[project_id]][[as.character(noise_alpha)]][[comp_type]][[norm_key]]$stage_topN[[as.character(top_n)]][[stage]] <- list(
              pct_in_edger = pct_edger_top,
              pct_in_limma = pct_limma_top,
              n_to = length(to_genes_top)
            )
          }
        }
        
        # ===== CREATE DISTRIBUTION PLOTS (per normalization) =====
        cat("\n      --- Distribution plots ---\n")
        
        for (i in seq_along(to_data_list)) {
          norm_name <- names(to_data_list)[i]
          to_significant <- to_data_list[[i]]
          norm_key <- names(NORM_DISPLAY)[which(NORM_DISPLAY == norm_name)]
          
          # Filter TO All data for this combination
          to_all_filtered <- to_all_data %>%
            filter(
              norm_method == norm_key,
              comparison == comp_type,
              stage == stage
            )
          
          if (nrow(to_all_filtered) == 0) {
            next
          }
          
          # Create distribution plot
          create_distribution_plots(
            to_all_data = to_all_filtered,
            to_significant = to_significant,
            de_edger = edger_all,
            de_limma = limma_all,
            project_id = project_id,
            norm_method = norm_key,
            stage = stage,
            comp_type = comp_type,
            param_name = param_name,
            output_dir = output_dir
          )
        }
      }
    }
  }
  
  return(all_stats_list)
}

# ==================== MAIN EXECUTION ====================

if (sys.nframe() == 0) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("COMPARISON: TENSOROMICS vs EDGER/LIMMA-VOOM\n")
  cat("Parameter combinations:\n")
  for (params in PARAM_COMBINATIONS) {
    cat(sprintf("  - %s: noise=%.2f\n", params$name, params$noise))
  }
  cat(sprintf("Top-N values: %s\n", paste(TOP_N_VALUES, collapse = ", ")))
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Output directory
  output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "comparison_de")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Collect statistics across all projects
  all_stats_list <- list()
  
  # Process all projects
  for (cancer_name in names(CANCER_TYPES)) {
    project_id <- CANCER_TYPES[cancer_name]
    all_stats_list <- compare_tensoromics_with_de(project_id, output_dir, all_stats_list)
  }
  
  # ===== CREATE CAPPED HEATMAPS FOR TOP N VALUES =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("CREATING CAPPED HEATMAPS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  heatmaps_capped <- create_all_heatmaps(all_stats_list, output_dir, top_n_values = TOP_N_VALUES)
  
  # ===== CREATE UNCAPPED HEATMAPS (ALL GENES) =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("CREATING UNCAPPED HEATMAPS (ALL GENES)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  heatmaps_uncapped <- create_uncapped_heatmaps(all_stats_list, output_dir)
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("COMPARISON COMPLETE\n")
  cat(sprintf("Results in: %s\n", output_dir))
  cat(sprintf("Distribution plots in: %s\n", file.path(output_dir, "distribution_plots")))
  cat(paste(rep("=", 80), collapse = ""), "\n")
}