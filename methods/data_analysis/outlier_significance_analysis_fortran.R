#!/usr/bin/env Rscript
# outlier_significance_analysis.R
# Outlier Significance Analysis mit kNN-basierten Noise-p-Werten
# STAGE‑WEISE IMPLEMENTIERUNG mit unabhängigen Nachbarschaften für Tumor und Gesunde

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(patchwork)
library(data.table)
library(parallel)

source("config.R")
source("utils.R")
source("tensoromics_functions.R")

# Load Fortran shared library
dyn.load("noise_model.so")

# ==================== KONFIGURATION ====================

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

NORM_METHODS <- c("raw", "std_log", "full", "log")
NORM_DISPLAY <- c(
  "raw" = "mean(TPM)",
  "log" = "log(mean(TPM))",
  "std_log" = "log(mean(gene-wise scaled TPM))",
  "full" = "log(mean(quantile(gene-wise scaled TPM)))"
)

# Alle drei Vergleichstypen
COMP_TYPES <- c("own_healthy", "family_mean", "ortholog_mean")
COMP_DISPLAY <- c(
  "own_healthy" = "Gene vs own healthy",
  "family_mean" = "Gene vs family mean",
  "ortholog_mean" = "Gene vs ortholog mean"
)

MIN_FAMILY_SIZES <- c(2, 3, 4)

# Für Plot 4: Speichere Perzentile für globale Visualisierung (optional)
STORE_NULL_PERCENTILES <- TRUE

# ==================== DATEN LADEN ====================

load_all_gene_results <- function() {
  all_results <- list()
  
  for (cancer_name in names(CANCER_TYPES)) {
    project_id <- CANCER_TYPES[cancer_name]
    cat(sprintf("\nLade %s (%s)...\n", cancer_name, project_id))
    
    for (norm_method in NORM_METHODS) {
      norm_suffix <- get_norm_suffix(norm_method)
      
      # Lade Gene Results (normalisiert)
      gene_file <- file.path(
        get_norm_output_dir(GENE_OUTPUT_DIR, norm_method),
        project_id,
        paste0(project_id, "_gene_results", norm_suffix, ".rds")
      )
      
      if (!file.exists(gene_file)) {
        cat(sprintf("  Datei nicht gefunden: %s\n", gene_file))
        next
      }
      
      gene_results <- readRDS(gene_file)
      
      # Lade Family LOESS Stats
      family_file <- file.path(
        get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method),
        project_id,
        paste0(project_id, "_family_loess_stats", norm_suffix, ".rds")
      )
      
      if (!file.exists(family_file)) {
        cat(sprintf("  Family Stats nicht gefunden: %s\n", family_file))
        next
      }
      
      family_stats <- readRDS(family_file)
      
      # Lade gesunde Rohdaten (für Noise) – einmal pro Normmethode
      healthy_data_raw <- load_stage_data(
        project_id, 
        STAGES[1], 
        "healthy", 
        use_constant_healthy = TRUE, 
        norm_method = "raw", 
        apply_mean = FALSE,
        normalize = FALSE   # Rohdaten ohne Normalisierung
      )
      
      if (is.null(healthy_data_raw) || !is.matrix(healthy_data_raw$expression_vectors)) {
        cat(sprintf("  Healthy RAW Data nicht gefunden für %s\n", norm_method))
        next
      }
      
      # Preprocess healthy replicates (einmal pro Normmethode)
      cat("    Preprocessing healthy replicates for norm:", norm_method, "\n")
      healthy_preproc <- preprocess_replicates(healthy_data_raw$expression_vectors, norm_method)
      
      all_results[[paste(project_id, norm_method, sep = "_")]] <- list(
        gene_results = gene_results,
        family_stats = family_stats,
        healthy_replicates_raw = healthy_data_raw$expression_vectors,
        healthy_preproc = healthy_preproc,
        gene_to_fam = family_stats$gene_to_fam,
        n_families = family_stats$n_families,
        cancer_type = cancer_name,
        cancer_id = project_id,
        norm_method = norm_method,
        norm_display = NORM_DISPLAY[norm_method]
      )
    }
  }
  
  cat(sprintf("\n Geladen: %d Datensätze\n", length(all_results)))
  return(all_results)
}

extract_distance_pvalues <- function(gene_results) {
  
  n_genes <- nrow(gene_results)
  n_stages <- length(STAGES)
  n_comparisons <- length(COMP_TYPES)
  
  distance_pvalues <- array(NA, dim = c(n_genes, n_stages, n_comparisons),
                           dimnames = list(
                             gene_results$gene_id,
                             STAGES,
                             COMP_TYPES
                           ))
  
  for (s in 1:n_stages) {
    stage <- STAGES[s]
    stage_col <- gsub(" ", "_", stage)
    
    for (comp_type in COMP_TYPES) {
      p_col <- switch(comp_type,
        "own_healthy" = paste0("p_value_distance_own_", stage_col),
        "family_mean" = paste0("p_value_distance_fam_", stage_col),
        "ortholog_mean" = paste0("p_value_distance_orth_", stage_col)
      )
      
      if (p_col %in% colnames(gene_results)) {
        distance_pvalues[, s, comp_type] <- gene_results[[p_col]]
      }
    }
  }
  
  return(distance_pvalues)
}

# ==================== NEUE HELFER FÜR NOISE MODELL ====================

#' Preprocess replicate matrix into pre‑log space, means, residuals, and sorted order
preprocess_replicates <- function(expr_matrix, norm_method) {
  # expr_matrix: samples × genes (raw TPM values)
  # Returns list with:
  #   prelog: matrix (samples × genes) of normalized values before log transformation
  #   means: vector of per‑gene means in prelog space
  #   residuals: list of vectors (residuals per gene) = prelog - mean
  #   replicate_counts: number of replicates per gene (non‑NA)
  #   sorted_order: order of genes by mean (increasing)
  #   means_sorted: sorted means
  #   residuals_sorted: residuals in same order as sorted_order

  # Apply normalisation up to pre‑log step
  if (norm_method == "raw") {
    prelog <- expr_matrix
  } else if (norm_method == "log") {
    prelog <- expr_matrix   # log will be applied later, residuals are on raw scale
  } else if (norm_method == "std_log") {
    # Scale gene‑wise (no quantile)
    original_names <- colnames(expr_matrix)
    expr_t <- t(expr_matrix)
    scaled <- tox_normalize_by_std_dev(expr_t)
    prelog <- t(scaled)
    colnames(prelog) <- original_names
  } else if (norm_method == "full") {
    original_names <- colnames(expr_matrix)
    expr_t <- t(expr_matrix)
    scaled <- tox_normalize_by_std_dev(expr_t)
    quant <- tox_quantile_normalization(scaled)
    prelog <- t(quant)
    colnames(prelog) <- original_names
  } else {
    stop("Unknown norm_method")
  }

  # Get total number of samples
  n_samples <- nrow(prelog)
  
  # Compute means (across replicates, ignoring NAs)
  means <- colMeans(prelog, na.rm = TRUE)
  names(means) <- colnames(prelog)
  
  # Count valid replicates per gene (non-NA values)
  replicate_counts <- colSums(!is.na(prelog))
  names(replicate_counts) <- colnames(prelog)
  
  # Compute residuals per gene
  residuals <- vector("list", ncol(prelog))
  for (g in seq_len(ncol(prelog))) {
    vec <- prelog[, g]
    vec_clean <- vec[!is.na(vec)]
    if (length(vec_clean) > 0) {
      m <- means[g]
      resid <- vec_clean - m
      residuals[[g]] <- resid
    } else {
      residuals[[g]] <- numeric(0)
    }
  }
  names(residuals) <- colnames(prelog)

  # Keep only genes with at least one valid replicate and non-NA mean
  valid_genes <- !is.na(means) & replicate_counts > 0
  
  # Filter to valid genes only
  means_valid <- means[valid_genes]
  replicate_counts_valid <- replicate_counts[valid_genes]
  residuals_valid <- residuals[valid_genes]
  
  # Sort by mean
  sorted_order <- order(means_valid, na.last = NA)
  means_sorted <- means_valid[sorted_order]
  replicate_counts_sorted <- replicate_counts_valid[sorted_order]
  residuals_sorted <- residuals_valid[sorted_order]
  original_indices <- which(valid_genes)[sorted_order]
  list(
    prelog = prelog,
    means = means,
    residuals = residuals,
    replicate_counts = replicate_counts,  # All genes, with 0 for those with no data
    sorted_order = original_indices,
    means_sorted = means_sorted,
    replicate_counts_sorted = replicate_counts_sorted,
    residuals_sorted = residuals_sorted
  )
}

#' R wrapper for Fortran noise model
compute_noise_pvalues <- function(
    cancer_preproc, healthy_preproc, gene_results, family_stats, gene_to_fam,
    stage, norm_method,
    B = 10000,
    k_start = 32, k_step = 16, k_max = 1024, tau = 0.15,
    max_pool_size = 50000
) {
  
  # Prepare input arrays
  n_genes <- nrow(gene_results)
  n_families <- length(family_stats$gene_family_means_all)
  
  # Cancer data - keep as samples × genes (column-major matches Fortran)
  cancer_means <- cancer_preproc$means
  cancer_replicates <- cancer_preproc$prelog  # samples × genes
  cancer_n_genes <- ncol(cancer_replicates)
  cancer_n_samples <- nrow(cancer_replicates)
  
  # Healthy data - keep as samples × genes
  healthy_means <- healthy_preproc$means
  healthy_replicates <- healthy_preproc$prelog  # samples × genes
  healthy_n_genes <- ncol(healthy_replicates)
  healthy_n_samples <- nrow(healthy_replicates)
  
  # Observed differences
  stage_col <- gsub(" ", "_", stage)
  obs_own <- gene_results[[paste0("shift_vs_own_healthy_unscaled_", stage_col)]]
  obs_fam <- gene_results[[paste0("shift_vs_family_mean_unscaled_", stage_col)]]
  obs_orth <- gene_results[[paste0("shift_vs_ortholog_mean_unscaled_", stage_col)]]

  valid_genes_own <- as.integer(!is.na(obs_own))
  valid_genes_fam <- as.integer(!is.na(obs_fam))
  valid_genes_orth <- as.integer(!is.na(obs_orth))

  obs_own[is.na(obs_own)] <- 0
  obs_fam[is.na(obs_fam)] <- 0
  obs_orth[is.na(obs_orth)] <- 0
  
  # Family data
  family_means <- family_stats$gene_family_means_all
  ortholog_means <- family_stats$gene_family_means_ortholog
  
  # Calculate family sizes (number of genes per family)
  family_sizes <- tabulate(gene_to_fam, nbins = n_families)
  
  # Count orthologs (genes with is_ortholog = TRUE)
  is_ortholog_sum <- sum(gene_results$is_ortholog, na.rm = TRUE)
  
  # Convert norm_method to integer
  norm_method_int <- switch(norm_method,
    "raw" = 0,
    "log" = 1,
    "std_log" = 2,
    "full" = 3,
    0  # default
  )
  
  # Output arrays
  pvalues_own <- numeric(n_genes)
  pvalues_fam <- numeric(n_genes)
  pvalues_orth <- numeric(n_genes)
  n_success <- integer(1)
  neighborhood_size_own <- integer(n_genes)
  neighborhood_size_fam <- integer(n_genes)
  neighborhood_size_orth <- integer(n_genes)
  neighborhood_size_cancer <- integer(n_genes)
  
  # Call Fortran subroutine
  result <- .Fortran("compute_noise_pvalues_fortran",
    cancer_means = as.double(cancer_means),
    cancer_replicates = as.double(cancer_replicates),
    cancer_n_genes = as.integer(cancer_n_genes),
    cancer_n_samples = as.integer(cancer_n_samples),
    healthy_means = as.double(healthy_means),
    healthy_replicates = as.double(healthy_replicates),
    healthy_n_genes = as.integer(healthy_n_genes),
    healthy_n_samples = as.integer(healthy_n_samples),
    obs_own = as.double(obs_own),
    obs_fam = as.double(obs_fam),
    obs_orth = as.double(obs_orth),
    family_means = as.double(family_means),
    ortholog_means = as.double(ortholog_means),
    valid_genes_own = as.integer(valid_genes_own),
    valid_genes_fam = as.integer(valid_genes_fam),
    valid_genes_orth = as.integer(valid_genes_orth),
    family_sizes = as.integer(family_sizes),
    is_ortholog_sum = as.integer(is_ortholog_sum),
    gene_to_family = as.integer(gene_to_fam),
    n_genes = as.integer(n_genes),
    n_families = as.integer(n_families),
    norm_method = as.integer(norm_method_int),
    B = as.integer(B),
    k_start = as.integer(k_start),
    k_step = as.integer(k_step),
    k_max = as.integer(k_max),
    tau = as.double(tau),
    pvalues_own = as.double(pvalues_own),
    pvalues_fam = as.double(pvalues_fam),
    pvalues_orth = as.double(pvalues_orth),
    n_success = as.integer(n_success),
    max_pool_size = as.integer(max_pool_size),
    neighborhood_size_own = as.integer(neighborhood_size_own),
    neighborhood_size_fam = as.integer(neighborhood_size_fam),
    neighborhood_size_orth = as.integer(neighborhood_size_orth),
    neighborhood_size_cancer = as.integer(neighborhood_size_cancer),
    PACKAGE = "noise_model"
  )

    # replace invalid value indicators with NA
    pvalues_own <- result$pvalues_own
    pvalues_fam <- result$pvalues_fam
    pvalues_orth <- result$pvalues_orth

    neighborhood_size_cancer <- result$neighborhood_size_cancer
    neighborhood_size_fam <- result$neighborhood_size_fam
    neighborhood_size_orth <- result$neighborhood_size_orth
    neighborhood_size_own <- result$neighborhood_size_own

    neighborhood_size_cancer[neighborhood_size_cancer < 0] <- NA
    neighborhood_size_fam[neighborhood_size_fam < 0] <- NA
    neighborhood_size_orth[neighborhood_size_orth < 0] <- NA
    neighborhood_size_own[neighborhood_size_own < 0] <- NA

    pvalues_own[pvalues_own < 0 | pvalues_own > 1] <- NA
    pvalues_fam[pvalues_fam < 0 | pvalues_fam > 1] <- NA
    pvalues_orth[pvalues_orth < 0 | pvalues_orth > 1] <- NA

    # Now create the matrix
    pvalues_matrix <- cbind(
        own_healthy = pvalues_own,
        family_mean = pvalues_fam,
        ortholog_mean = pvalues_orth
    )
    rownames(pvalues_matrix) <- gene_results$gene_id
    
    neighborhood_matrix <- cbind(
      neighborhood_size_cancer = neighborhood_size_cancer,
      neighborhood_size_fam = neighborhood_size_fam,
      neighborhood_size_orth = neighborhood_size_orth,
      neighborhood_size_own = neighborhood_size_own
    )
    rownames(neighborhood_matrix) <- gene_results$gene_id
  
  list(
    pvalues = pvalues_matrix,
    n_success = result$n_success,
    neighborhoods = neighborhood_matrix,
    null_sample = numeric(0)  # For compatibility with existing code
  )
}

# ==================== GESAMTTABELLE ERSTELLEN (MIT STAGE-WEISER VERARBEITUNG) ====================

process_single_combination <- function(cancer_type, project_id, norm_method, stage, 
                                       gene_results, family_stats, healthy_preproc,
                                       gene_to_fam, n_families, norm_display,
                                       output_dir) {
  
  # Load cancer data for this stage
  cancer_data_raw <- load_stage_data(
    project_id = project_id,
    stage = stage,
    data_type = "cancer",
    use_constant_healthy = FALSE,
    norm_method = "raw",
    apply_mean = FALSE,
    normalize = FALSE
  )
  
  if (is.null(cancer_data_raw) || !is.matrix(cancer_data_raw$expression_vectors)) {
    return(NULL)
  }
  
  # Preprocess cancer replicates
  cancer_preproc <- preprocess_replicates(cancer_data_raw$expression_vectors, norm_method)
  
  # Compute noise p-values
  noise_result <- compute_noise_pvalues(
    cancer_preproc = cancer_preproc,
    healthy_preproc = healthy_preproc,
    gene_results = gene_results,
    family_stats = family_stats,
    gene_to_fam = gene_to_fam,
    stage = stage,
    norm_method = norm_method,
    B = 10000,
    k_start = 32, 
    k_step = 16, 
    k_max = 1024, 
    tau = 0.15,
    max_pool_size = 50000
  )
  
  # Extract distance p-values
  distance_p <- extract_distance_pvalues(gene_results)
  family_n_genes <- tabulate(gene_to_fam, nbins = n_families)
  n_genes <- nrow(gene_results)
  
  # Pre-allocate
  results_list <- list()
  row_idx <- 1
  
  # Get neighborhood size matrix
  neighborhood_mat <- noise_result$neighborhoods
  
  for (g in 1:n_genes) {
    gene_id <- gene_results$gene_id[g]
    for (comp_type in COMP_TYPES) {
      dist_p_val <- distance_p[g, stage, comp_type]
      noise_p_val <- noise_result$pvalues[g, comp_type]
      
      if (is.na(dist_p_val) || is.na(noise_p_val)) next
      
      # Scaled with families sds
      shift_col <- switch(comp_type,
                          "own_healthy" = paste0("shift_vs_own_healthy_", gsub(" ", "_", stage)),
                          "family_mean" = paste0("shift_vs_family_mean_", gsub(" ", "_", stage)),
                          "ortholog_mean" = paste0("shift_vs_ortholog_mean_", gsub(" ", "_", stage)))
      
      if (!shift_col %in% colnames(gene_results)) next
      
      shift <- gene_results[g, shift_col]
      direction <- ifelse(is.na(shift), NA, ifelse(shift > 0, "up", "down"))
      
      fam_id <- gene_to_fam[g]
      family_size <- if (!is.na(fam_id) && fam_id > 0) family_n_genes[fam_id] else NA
      family_mean_val <- if (!is.na(fam_id) && fam_id > 0) 
        family_stats$gene_family_means_all[fam_id] else NA
      
      # Get neighborhood sizes for this gene
      n_cancer <- neighborhood_mat[g, "neighborhood_size_cancer"]
      n_healthy <- switch(comp_type,
                          "own_healthy" = neighborhood_mat[g, "neighborhood_size_own"],
                          "family_mean" = neighborhood_mat[g, "neighborhood_size_fam"],
                          "ortholog_mean" = neighborhood_mat[g, "neighborhood_size_orth"])
      
      results_list[[row_idx]] <- data.frame(
        gene_id = gene_id,
        cancer_type = cancer_type,
        cancer_id = project_id,
        stage = stage,
        normalization = norm_display,
        norm_method = norm_method,
        comparison = comp_type,
        direction = direction,
        signed_difference = shift,
        distance_p_value = dist_p_val,
        noise_p_value = noise_p_val,
        family_size = family_size,
        family_mean = family_mean_val,
        neighborhood_size_cancer = n_cancer,
        neighborhood_size_healthy = n_healthy,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }
  
  # Combine results for this combination
  if (length(results_list) == 0) {
    return(NULL)
  }
  
  result_df <- do.call(rbind, results_list)
  
  # ===== SAVE INTERMEDIATE RESULTS =====
  safe_stage <- gsub(" ", "_", stage)
  filename <- sprintf("%s_%s_%s_results.rds", 
                      project_id, norm_method, safe_stage)
  intermediate_file <- file.path(output_dir, "intermediate", filename)
  
  dir.create(dirname(intermediate_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(result_df, intermediate_file)
  
  csv_file <- file.path(output_dir, "intermediate", 
                        gsub(".rds", ".csv", filename))
  write.csv(result_df, csv_file, row.names = FALSE)
  
  done_file <- file.path(output_dir, "intermediate", "done", 
                         sprintf("%s_%s_%s.done", 
                                 project_id, norm_method, safe_stage))
  dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste(
    "Completed:", Sys.time(),
    "\nGenes:", nrow(gene_results),
    "\nSuccessfully processed:", noise_result$n_success,
    "\nRows in output:", nrow(result_df)
  ), done_file)
  
  return(result_df)
}

create_all_genes_table <- function(all_results, n_cores = 64, output_dir) {
  
  cat("\n>>> Erstelle Gesamttabelle (stage-weise) mit paralleler Verarbeitung...\n")
  cat(sprintf("    Verwende %d Kerne für parallele Kombinationen\n", n_cores))
  cat(sprintf("    Intermediate Ergebnisse werden gespeichert in: %s/intermediate/\n", output_dir))
  
  # Create all combination tasks
  tasks <- list()
  task_id <- 1
  
  for (res_name in names(all_results)) {
    res <- all_results[[res_name]]
    
    for (stage in STAGES) {
      tasks[[task_id]] <- list(
        cancer_type = res$cancer_type,
        project_id = res$cancer_id,
        norm_method = res$norm_method,
        norm_display = res$norm_display,
        stage = stage,
        gene_results = res$gene_results,
        family_stats = res$family_stats,
        healthy_preproc = res$healthy_preproc,
        gene_to_fam = res$gene_to_fam,
        n_families = res$n_families,
        output_dir = output_dir  # Pass output directory
      )
      task_id <- task_id + 1
    }
  }
  
  cat(sprintf("    Gesamt: %d Kombinationen zu verarbeiten\n", length(tasks)))

  combination_results <- mclapply(tasks, function(task) {
    tryCatch({
    process_single_combination(
        cancer_type = task$cancer_type,
        project_id = task$project_id,
        norm_method = task$norm_method,
        norm_display = task$norm_display,
        stage = task$stage,
        gene_results = task$gene_results,
        family_stats = task$family_stats,
        healthy_preproc = task$healthy_preproc,
        gene_to_fam = task$gene_to_fam,
        n_families = task$n_families,
        output_dir = task$output_dir
    )
    }, error = function(e) {
    cat(sprintf("Error processing %s %s %s: %s\n", 
                task$cancer_type, task$norm_method, task$stage, e$message))
    return(NULL)
    })
  }, mc.cores = n_cores, mc.preschedule = TRUE)
  
  # Combine all results
  cat("\n>>> Kombiniere Ergebnisse...\n")
  valid_results <- combination_results[!sapply(combination_results, is.null)]
  cat(sprintf("    Erfolgreich verarbeitet: %d / %d Kombinationen\n", 
              length(valid_results), length(tasks)))
  
  if (length(valid_results) > 0) {
    final_result <- rbindlist(valid_results, fill = TRUE)
    cat(sprintf("    Finale Tabellengröße: %d Zeilen\n", nrow(final_result)))
    return(final_result)
  } else {
    stop("Keine erfolgreichen Ergebnisse!")
  }
}
# ==================== PLOT FUNKTIONEN ====================

plot_jaccard_heatmap <- function(data, title = "", subtitle = "") {
  
  norm_methods <- c("raw", "log", "std_log", "full")
  norm_sets <- list()
  
  for (norm in norm_methods) {
    norm_data <- data %>% filter(norm_method == norm)
    norm_sets[[norm]] <- unique(paste(norm_data$gene_id, norm_data$cancer_id, norm_data$stage))
  }
  
  jaccard_mat <- matrix(NA, 4, 4)
  rownames(jaccard_mat) <- norm_methods
  colnames(jaccard_mat) <- norm_methods
  
  for (i in 1:4) {
    for (j in i:4) {
      if (length(norm_sets[[i]]) > 0 && length(norm_sets[[j]]) > 0) {
        intersection <- length(intersect(norm_sets[[i]], norm_sets[[j]]))
        union <- length(union(norm_sets[[i]], norm_sets[[j]]))
        jaccard_mat[i, j] <- intersection / union
        jaccard_mat[j, i] <- jaccard_mat[i, j]
      }
    }
  }
  
  plot_data <- as.data.frame(as.table(jaccard_mat))
  colnames(plot_data) <- c("Method1", "Method2", "Jaccard")
  plot_data <- plot_data %>% filter(!is.na(Jaccard))
  
  ggplot(plot_data, aes(x = Method1, y = Method2, fill = Jaccard)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 5) +
    scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
    labs(title = title, subtitle = subtitle, x = "", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# ==================== HAUPTSCHLEIFE ====================

analyze_outlier_significance <- function(output_dir = NULL, n_cores = 64) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("OUTLIER SIGNIFICANCE ANALYSIS (parallel combinations)\n")
  cat(sprintf("Using %d cores for parallel combination processing\n", n_cores))
  cat("Parameters: k_start=32, k_step=16, k_max=1024, tau=0.15, B=10000\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Load all data (this is still sequential, but only once)
  all_results <- load_all_gene_results()
  if (length(all_results) == 0) stop("Keine Daten gefunden!")
  
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create results table in parallel
  cat("\n>>> Erstelle Gesamttabelle (parallel)...\n")
  all_genes_table <- create_all_genes_table(all_results, n_cores = n_cores, output_dir = output_dir)
  
  all_genes_table <- all_genes_table %>%
    group_by(cancer_id, normalization, stage, comparison) %>%
    mutate(
      noise_p_value_adj = p.adjust(noise_p_value, method = "BH")
    ) %>%
    ungroup()

  # Save results
  saveRDS(all_genes_table, file.path(output_dir, "all_genes_pvalues.rds"))
  write.csv(all_genes_table, file.path(output_dir, "all_genes_pvalues.csv"), row.names = FALSE)
  
  # ===== 3. BASELINE AUSGABE =====
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("UNKORRIGIERTE WERTE (BASELINE)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  cat("\n--- Summary distance_p_values ---\n")
  print(summary(all_genes_table$distance_p_value))
  cat("\n--- Summary noise_p_values ---\n")
  print(summary(all_genes_table$noise_p_value))
  
  # ===== 4. SELEKTION: NOISE-P-WERT < 0.05 =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SELEKTION: noise_p_adj < 0.01 (sortiert nach distance_p_value)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  ALPHA_NOISE <- 0.01
  
  significant_by_norm <- list()
  
  for (current_norm in unique(all_genes_table$norm_method)) {
    
    cat(sprintf("\n>>> Verarbeite Normalisierung: %s <<<\n", current_norm))
    
    norm_data <- all_genes_table %>% 
      filter(norm_method == current_norm)
    
    # ===== DIAGNOSE PRO VERGLEICHSTYP =====
    cat(sprintf("\n  DIAGNOSE %s:\n", current_norm))
    
    for (comp_type in COMP_TYPES) {
      comp_data <- norm_data %>% filter(comparison == comp_type)
      
      cat(sprintf("\n    --- %s ---\n", comp_type))
      cat(sprintf("      Noise p (unkorr) < 0.01: %d (%.1f%%)\n", 
                  sum(comp_data$noise_p_value < 0.01, na.rm = TRUE),
                  100 * mean(comp_data$noise_p_value < 0.01, na.rm = TRUE)))
      cat(sprintf("      Noise p (adj) < 0.01: %d (%.1f%%)\n\n", 
                  sum(comp_data$noise_p_value_adj < 0.01, na.rm = TRUE),
                  100 * mean(comp_data$noise_p_value_adj < 0.01, na.rm = TRUE)))
      cat("      Summary of noise_p_value: \n")
      print(summary(comp_data$noise_p_value))
      cat("      Summary of noise_p_value_adj: \n")
      print(summary(comp_data$noise_p_value_adj))
    }
    
    # ===== SELEKTION: Nur noise_p_adj < 0.01 =====
    significant <- norm_data %>%
      filter(noise_p_value_adj < ALPHA_NOISE) %>%
      arrange(distance_p_value)
    
    cat(sprintf("\n  → Signifikante Einträge (noise_p_adj < 0.01): %d\n", nrow(significant)))
    
    if (nrow(significant) > 0) {
      cat("\n  Top 10 nach distance_p_value:\n")
      print(significant %>%
        select(gene_id, stage, comparison, distance_p_value, noise_p_value_adj) %>%
        head(10))
      
      filename <- sprintf("significant_%s.csv", 
                         gsub("[^a-zA-Z0-9]", "_", current_norm))
      write.csv(significant, file.path(output_dir, filename), row.names = FALSE)
      
      significant_by_norm[[current_norm]] <- significant
    }
  }
  
  # ===== 5. GESAMTTABELLE ÜBER ALLE NORMALISIERUNGEN =====
  cat("\n>>> Gesamttabelle über alle Normalisierungen <<<\n")
  
  if (length(significant_by_norm) > 0) {
    all_significant <- rbindlist(significant_by_norm, fill = TRUE)
    all_significant <- all_significant %>% arrange(distance_p_value)
    
    cat(sprintf("\n  Gesamt: %d signifikante Einträge (noise_p_adj < 0.01)\n", nrow(all_significant)))
    
    # Statistik pro Normalisierung
    cat("\n  Statistik pro Normalisierung:\n")
    norm_stats <- all_significant %>%
      group_by(norm_method) %>%
      summarise(
        n = n(),
        median_distance_p = median(distance_p_value, na.rm = TRUE),
        q95_distance_p = quantile(distance_p_value, 0.95, na.rm = TRUE),
        min_distance_p = min(distance_p_value, na.rm = TRUE)
      ) %>%
      arrange(desc(n))
    print(norm_stats)
    
    # Statistik pro Vergleichstyp
    cat("\n  Statistik pro Vergleichstyp:\n")
    comp_stats <- all_significant %>%
      group_by(comparison) %>%
      summarise(
        n = n(),
        median_distance_p = median(distance_p_value, na.rm = TRUE),
        q95_distance_p = quantile(distance_p_value, 0.95, na.rm = TRUE)
      ) %>%
      arrange(desc(n))
    print(comp_stats)
    
    # Statistik pro Krebsart
    cat("\n  Statistik pro Krebsart (Top 10):\n")
    cancer_stats <- all_significant %>%
      group_by(cancer_type) %>%
      summarise(
        n = n(),
        median_distance_p = median(distance_p_value, na.rm = TRUE)
      ) %>%
      arrange(desc(n)) %>%
      head(10)
    print(cancer_stats)
    
    write.csv(all_significant, 
              file.path(output_dir, "significant_all.csv"), 
              row.names = FALSE)
    saveRDS(all_significant, 
            file.path(output_dir, "significant_all.rds"))
  } else {
    cat("\n  Keine signifikanten Einträge gefunden!\n")
    all_significant <- data.frame()
  }
  
  # ===== 6. JACCARD PLOTS =====
  cat("\n>>> Jaccard Plots (basierend auf noise_p_adj < 0.01)\n")
  
  plots_dir <- file.path(output_dir, "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (nrow(all_significant) > 0) {
    for (min_size in MIN_FAMILY_SIZES) {
      for (comp_type in COMP_TYPES) {
        
        plot_data <- all_significant %>%
          filter(
            comparison == comp_type, 
            family_size >= min_size
          )
        
        if (nrow(plot_data) > 0) {
          p_adj <- plot_jaccard_heatmap(
            plot_data,
            title = sprintf("%s (min %d) - KORRIGIERT", comp_type, min_size),
            subtitle = sprintf("noise_p_adj < 0.01, n=%d", nrow(plot_data))
          )
          
          filename_adj <- sprintf("jaccard_ADJ_%s_min%d.png", comp_type, min_size)
          ggsave(file.path(plots_dir, filename_adj), p_adj, 
                 width = 7, height = 6, dpi = 300)
          
          # Rohdaten für Vergleich
          raw_data <- all_genes_table %>%
            filter(
              comparison == comp_type,
              family_size >= min_size,
              noise_p_value < 0.01
            )
          
          if (nrow(raw_data) > 0 && nrow(plot_data) > 0) {
            norm_methods <- c("raw", "log", "std_log", "full")
            
            raw_sets <- list()
            for (norm in norm_methods) {
              norm_data <- raw_data %>% filter(norm_method == norm)
              raw_sets[[norm]] <- unique(paste(norm_data$gene_id, 
                                               norm_data$cancer_id, 
                                               norm_data$stage))
            }
            
            adj_sets <- list()
            for (norm in norm_methods) {
              norm_data <- plot_data %>% filter(norm_method == norm)
              adj_sets[[norm]] <- unique(paste(norm_data$gene_id, 
                                               norm_data$cancer_id, 
                                               norm_data$stage))
            }
            
            comparison_mat <- matrix(NA, 4, 4)
            rownames(comparison_mat) <- norm_methods
            colnames(comparison_mat) <- norm_methods
            
            for (i in 1:4) {
              for (j in 1:4) {
                if (length(raw_sets[[i]]) > 0 && length(adj_sets[[j]]) > 0) {
                  intersection <- length(intersect(raw_sets[[i]], adj_sets[[j]]))
                  union <- length(union(raw_sets[[i]], adj_sets[[j]]))
                  comparison_mat[i, j] <- intersection / union
                }
              }
            }
            
            plot_data_comp <- as.data.frame(as.table(comparison_mat))
            colnames(plot_data_comp) <- c("Raw_Norm", "Adj_Norm", "Jaccard")
            plot_data_comp <- plot_data_comp %>% filter(!is.na(Jaccard))
            
            p_comp <- ggplot(plot_data_comp, aes(x = Adj_Norm, y = Raw_Norm, fill = Jaccard)) +
              geom_tile() +
              geom_text(aes(label = sprintf("%.2f", Jaccard)), size = 4) +
              scale_fill_gradient(low = "white", high = "steelblue", limits = c(0, 1)) +
              labs(title = sprintf("%s (min %d) - Raw data vs. noise_p_adj < 0.01", 
                                   comp_type, min_size),
                   x = "noise_p_adj < 0.01", y = "Raw data (noise_p < 0.01)") +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            
            filename_comp <- sprintf("jaccard_COMP_%s_min%d.png", comp_type, min_size)
            ggsave(file.path(plots_dir, filename_comp), p_comp, 
                   width = 7, height = 6, dpi = 300)
          }
        }
      }
    }
  }
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("ANALYSE ABGESCHLOSSEN\n")
  cat(sprintf("Ergebnisse in: %s\n", output_dir))
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  return(invisible(list(
    data = all_genes_table,
    significant = all_significant
  )))
}

# ==================== AUSFÜHRUNG ====================

if (sys.nframe() == 0) {
  set.seed(42)
  base_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")
  
  cat("\n", paste(rep("#", 80), collapse = ""), "\n")
  cat("### ADAPTIVE NOISE MODEL ANALYSIS (STAGE-WISE) ###\n")
  cat("Parameters: k_start=32, k_step=16, k_max=1024, tau=0.15, B=10000\n")
  cat(paste(rep("#", 80), collapse = ""), "\n")
  
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- analyze_outlier_significance(output_dir = base_dir)
  
  saveRDS(list(
    timestamp = Sys.time(),
    n_significant = nrow(results$significant),
    parameters = list(
      k_start = 32,
      k_step = 16,
      k_max = 1024,
      tau = 0.15,
      B = 10000
    )
  ), file.path(base_dir, "metadata.rds"))
}