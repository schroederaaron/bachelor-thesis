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
  "std_log" = "log(gene-wise scaled TPM)",
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
        normalize = FALSE
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

#' Adaptive gathering of residuals from genes with similar mean expression

#' @param target_mean The mean expression to find neighbors around
#' @param means_sorted Sorted vector of gene mean expressions
#' @param residuals_sorted List of residual vectors in same order as means_sorted
#' @param k_start Initial minimum number of residuals to collect
#' @param k_step Minimum number of new residuals to add per expansion step
#' @param k_max Maximum number of residuals to collect
#' @param tau Stopping criterion for relative change in mean absolute residual
#' @return Vector of pooled residuals
gather_residuals <- function(target_mean, means_sorted, residuals_sorted,
                             k_start = 32, k_step = 16, k_max = 1024,
                             tau = 0.15) {
  
  # Step 1: Find closest index using binary search
  n_genes <- length(means_sorted)
  if (n_genes == 0) return(numeric(0))
  
  # Binary search to find insertion point
  pos <- findInterval(target_mean, means_sorted)
  
  # Determine the closest index
  if (pos == 0) {
    pos <- 1
  } else if (pos == n_genes) {
    pos <- n_genes
  } else {
    # Check which neighbor is closer
    if (abs(means_sorted[pos] - target_mean) > abs(means_sorted[pos + 1] - target_mean)) {
      pos <- pos + 1
    }
  }
  
  # Track which genes we've included
  left <- pos
  right <- pos
  included_indices <- pos
  
  # Function to collect all residuals from a set of gene indices
  get_pool_from_indices <- function(indices) {
    all_res <- c()
    for (idx in indices) {
      all_res <- c(all_res, residuals_sorted[[idx]])
    }
    all_res
  }
  
  # Start with the closest gene
  pool <- residuals_sorted[[pos]]
  current_size <- length(pool)
  
  # Keep track of neighbors we need to add in order of increasing distance
  left_candidate <- pos - 1
  right_candidate <- pos + 1
  
  # Precompute distances to target for efficient neighbor selection
  # Since we'll be adding neighbors one by one, we can compute distances on the fly
  
  # Step 2: Initial expansion to reach at least k_start residuals
  while (current_size < k_start && (left_candidate >= 1 || right_candidate <= n_genes)) {
    # Determine which candidate is closer in mean expression space
    left_dist <- if (left_candidate >= 1) abs(means_sorted[left_candidate] - target_mean) else Inf
    right_dist <- if (right_candidate <= n_genes) abs(means_sorted[right_candidate] - target_mean) else Inf
    
    if (left_dist <= right_dist && left_candidate >= 1) {
      # Add the left neighbor
      pool <- c(pool, residuals_sorted[[left_candidate]])
      current_size <- length(pool)
      left_candidate <- left_candidate - 1
    } else if (right_candidate <= n_genes) {
      # Add the right neighbor
      pool <- c(pool, residuals_sorted[[right_candidate]])
      current_size <- length(pool)
      right_candidate <- right_candidate + 1
    } else {
      # No more candidates
      break
    }
  }
  
  # If still not enough residuals, return what we have (or empty if too few)
  if (current_size < 10) {
    return(numeric(0))
  }
  
  # Compute initial mean absolute residual for adaptive growth
  S_old <- mean(abs(pool))
  if (S_old == 0) {
    return(pool)  # No variability, can't expand meaningfully
  }
  
  # Step 3: Adaptive growth with k_step increments
  while (current_size < k_max && (left_candidate >= 1 || right_candidate <= n_genes)) {
    
    # We will try to add at least k_step new residuals
    new_pool <- pool
    added_this_round <- 0
    genes_added <- 0
    
    # Keep adding neighbors (one at a time, alternating by distance) until we've added at least k_step residuals
    # or run out of neighbors
    while (added_this_round < k_step && (left_candidate >= 1 || right_candidate <= n_genes)) {
      
      left_dist <- if (left_candidate >= 1) abs(means_sorted[left_candidate] - target_mean) else Inf
      right_dist <- if (right_candidate <= n_genes) abs(means_sorted[right_candidate] - target_mean) else Inf
      
      if (left_dist <= right_dist && left_candidate >= 1) {
        # Add the left neighbor
        new_resids <- residuals_sorted[[left_candidate]]
        new_pool <- c(new_pool, new_resids)
        added_this_round <- added_this_round + length(new_resids)
        genes_added <- genes_added + 1
        left_candidate <- left_candidate - 1
      } else if (right_candidate <= n_genes) {
        # Add the right neighbor
        new_resids <- residuals_sorted[[right_candidate]]
        new_pool <- c(new_pool, new_resids)
        added_this_round <- added_this_round + length(new_resids)
        genes_added <- genes_added + 1
        right_candidate <- right_candidate + 1
      } else {
        # No more candidates
        break
      }
    }
    
    # If no genes were added, break
    if (genes_added == 0) {
      break
    }
    
    # Compute new mean absolute residual
    S_new <- mean(abs(new_pool))
    rel_change <- (S_new - S_old) / S_old
    
    # Check stopping condition
    if (rel_change > tau) {
      # Don't accept this expansion
      break
    }
    
    # Accept expansion
    pool <- new_pool
    current_size <- length(pool)
    S_old <- S_new
  }
  
  return(pool)
}

# ==================== NOISE P-WERT BERECHNUNG (STAGE-WEISE) ====================

compute_noise_pvalues_stage <- function(
    cancer_preproc, healthy_preproc, gene_results, family_stats, gene_to_fam,
    stage, norm_method,
    B = 10000,
    k_start = 32, k_step = 16, k_max = 1024, tau = 0.15, output_dir = NULL, cancer_id
) {
  set.seed(42)
  
  # Stage-specific column names
  stage_col <- gsub(" ", "_", stage)
  own_col <- paste0("shift_vs_own_healthy_", stage_col)
  fam_col <- paste0("shift_vs_family_mean_", stage_col)
  orth_col <- paste0("shift_vs_ortholog_mean_", stage_col)
  
  if (!(own_col %in% colnames(gene_results))) {
    stop("Missing column: ", own_col)
  }
  
  n_genes <- nrow(gene_results)
  noise_pvalues <- array(NA, dim = c(n_genes, 3),
                         dimnames = list(gene_results$gene_id,
                                         c("own_healthy", "family_mean", "ortholog_mean")))
  
  # ========== PRE-COMPUTE CONSTANTS AND LOOKUPS ==========
  
  # Convert to vectors for faster access
  gene_ids <- gene_results$gene_id
  gene_to_fam_int <- as.integer(gene_to_fam)
  
  # Create index mappings (name to position)
  cancer_gene_names <- names(cancer_preproc$means)
  healthy_gene_names <- names(healthy_preproc$means)
  
  gene_to_cancer_idx <- match(gene_ids, cancer_gene_names)
  gene_to_healthy_idx <- match(gene_ids, healthy_gene_names)
  
  # Precompute family sizes for O(1) lookup
  family_sizes <- table(gene_to_fam)
  names(family_sizes) <- names(family_sizes)
  
  # Precompute ortholog sum (constant across genes)
  is_ortholog_sum <- sum(gene_results$is_ortholog, na.rm = TRUE)
  
  # Extract vectors from family_stats for faster access
  family_means_all <- family_stats$gene_family_means_all
  ortholog_means_all <- family_stats$gene_family_means_ortholog
  family_sds_all <- family_stats$gene_smoothed_sds_all
  family_sds_orth <- family_stats$gene_smoothed_sds_ortholog
  
  # Extract healthy data for faster access
  healthy_means <- healthy_preproc$means
  healthy_replicate_counts <- healthy_preproc$replicate_counts
  healthy_means_sorted <- healthy_preproc$means_sorted
  healthy_residuals_sorted <- healthy_preproc$residuals_sorted
  
  # Extract cancer data for faster access
  cancer_means <- cancer_preproc$means
  cancer_replicate_counts <- cancer_preproc$replicate_counts
  cancer_means_sorted <- cancer_preproc$means_sorted
  cancer_residuals_sorted <- cancer_preproc$residuals_sorted
  
  # Extract observed differences as vectors for direct access
  obs_own_vec <- gene_results[[own_col]]
  obs_fam_vec <- gene_results[[fam_col]]
  obs_orth_vec <- gene_results[[orth_col]]
  
  # ========== CACHES FOR RESIDUAL POOLS ==========
  # Note: These need to be accessible within process_gene
  # We'll use environments to allow modification from within the function
  family_pool_cache <- new.env()
  ortholog_pool_cache <- new.env()
  
  # ========== PROCESS GENE FUNCTION ==========
  process_gene <- function(g) {
    # Wrap everything in tryCatch to catch errors
    tryCatch({
      # Get indices
      cancer_idx <- gene_to_cancer_idx[g]
      if (is.na(cancer_idx)) {
        return(list(success = FALSE, reason = "no_gene_in_cancer"))
      }
      
      healthy_idx <- gene_to_healthy_idx[g]
      if (is.na(healthy_idx)) {
        return(list(success = FALSE, reason = "no_gene_in_healthy"))
      }
      
      # Get cancer values using index
      mu_cancer <- cancer_means[cancer_idx]
      r_cancer <- cancer_replicate_counts[cancer_idx]
      
      if (is.na(r_cancer) || r_cancer == 0) {
        return(list(success = FALSE, reason = "zero_replicates_cancer"))
      }
      
      # Get healthy values using index
      mu_healthy <- healthy_means[healthy_idx]
      r_healthy <- healthy_replicate_counts[healthy_idx]
      
      if (is.na(r_healthy) || r_healthy == 0) {
        return(list(success = FALSE, reason = "zero_replicates_healthy"))
      }
      
      # Family info using direct indexing
      fam_id <- gene_to_fam_int[g]
      if (is.na(fam_id) || fam_id == 0 || fam_id > length(family_means_all)) {
        return(list(success = FALSE, reason = "no_family"))
      }
      
      family_mean <- family_means_all[fam_id]
      ortholog_mean <- ortholog_means_all[fam_id]
      family_sd_all <- family_sds_all[fam_id]
      family_sd_orth <- family_sds_orth[fam_id]
      family_size <- family_sizes[as.character(fam_id)]
      
      # Quick validation
      if (is.na(family_sd_all) || is.na(family_sd_orth) || 
          is.null(family_size) || family_size == 0 || is_ortholog_sum == 0) {
        return(list(success = FALSE, reason = "invalid_family_data"))
      }
      
      # ========== COLLECT RESIDUAL POOLS ==========
      
      # Cancer residuals (per gene - cannot cache)
      cancer_resid_pool <- gather_residuals(mu_cancer,
                                            cancer_means_sorted,
                                            cancer_residuals_sorted,
                                            k_start, k_step, k_max, tau)
      if (length(cancer_resid_pool) < 10) {
        return(list(success = FALSE, reason = "cancer_resid_pool"))
      }
      
      # Healthy own residuals (per gene - cannot cache)
      healthy_resid_pool_own <- gather_residuals(mu_healthy,
                                                 healthy_means_sorted,
                                                 healthy_residuals_sorted,
                                                 k_start, k_step, k_max, tau)
      if (length(healthy_resid_pool_own) < 10) {
        return(list(success = FALSE, reason = "healthy_resid_pool_own"))
      }
      
      # Family pool with caching
      fam_key <- as.character(fam_id)
      if (!exists(fam_key, envir = family_pool_cache)) {
        family_pool_cache[[fam_key]] <- gather_residuals(family_mean,
                                                          healthy_means_sorted,
                                                          healthy_residuals_sorted,
                                                          k_start, k_step, k_max, tau)
      }
      healthy_resid_pool_fam <- family_pool_cache[[fam_key]]
      if (length(healthy_resid_pool_fam) < 10) {
        return(list(success = FALSE, reason = "healthy_resid_pool_fam"))
      }
      
      # Ortholog pool with caching
      if (!exists(fam_key, envir = ortholog_pool_cache)) {
        ortholog_pool_cache[[fam_key]] <- gather_residuals(ortholog_mean,
                                                            healthy_means_sorted,
                                                            healthy_residuals_sorted,
                                                            k_start, k_step, k_max, tau)
      }
      healthy_resid_pool_orth <- ortholog_pool_cache[[fam_key]]
      if (length(healthy_resid_pool_orth) < 10) {
        return(list(success = FALSE, reason = "healthy_resid_pool_orth"))
      }
      
      # Observed differences (direct vector access)
      obs_own <- obs_own_vec[g]
      obs_fam <- obs_fam_vec[g]
      obs_orth <- obs_orth_vec[g]
      
      if (is.na(obs_own) || is.na(obs_fam) || is.na(obs_orth)) {
        return(list(success = FALSE, reason = "missing_obs"))
      }
      
      # ========== VECTORIZED P-VALUE COMPUTATION ==========
      
      compute_pvalue_vec <- function(mu_c, r_c, resid_c,
                                     mu_h, r_h, resid_h,
                                     obs, sd_factor) {
        # Set seed for reproducibility (unique per gene)
        set.seed(42 + g * 1000)
        
        # Vectorized sampling - sample all residuals at once
        # Cancer side: sample r_c * B residuals, reshape, then rowMeans
        cancer_samples <- matrix(sample(resid_c, B * r_c, replace = TRUE), 
                                 nrow = B, ncol = r_c)
        eta_c <- rowMeans(cancer_samples)
        
        # Healthy side
        healthy_samples <- matrix(sample(resid_h, B * r_h, replace = TRUE), 
                                  nrow = B, ncol = r_h)
        eta_h <- rowMeans(healthy_samples)
        
        # Perturbed means
        x_c <- mu_c + eta_c
        x_h <- mu_h + eta_h
        
        # Calculate null distances based on norm_method
        if (norm_method == "raw") {
          null_dists <- abs(x_c - x_h)
        } else {
          # Log-transformation (ensure positivity)
          x_c <- pmax(x_c, 0) + 1
          x_h <- pmax(x_h, 0) + 1
          null_dists <- abs(log2(x_c) - log2(x_h))
        }
        
        # Scale by family SD
        if (!is.na(sd_factor) && sd_factor > 0) {
          null_dists <- null_dists / sd_factor
        }
        
        # Empirical p-value
        p <- (1 + sum(null_dists >= abs(obs))) / (1 + B)
        return(p)
      }
      
      # Compute p-values for all three comparisons
      p_own <- compute_pvalue_vec(mu_cancer, r_cancer, cancer_resid_pool,
                                  mu_healthy, r_healthy, healthy_resid_pool_own,
                                  obs_own, family_sd_all)
      
      p_fam <- compute_pvalue_vec(mu_cancer, r_cancer, cancer_resid_pool,
                                  family_mean, family_size, healthy_resid_pool_fam,
                                  obs_fam, family_sd_all)
      
      p_orth <- compute_pvalue_vec(mu_cancer, r_cancer, cancer_resid_pool,
                                   ortholog_mean, is_ortholog_sum, healthy_resid_pool_orth,
                                   obs_orth, family_sd_orth)
      
      return(list(
        success = TRUE,
        pvalues = c(own = p_own, fam = p_fam, orth = p_orth),
        reason = "success"
      ))
      
    }, error = function(e) {
      # Return error information if something goes wrong
      return(list(success = FALSE, reason = paste("error:", e$message)))
    })
  }
  
  # ========== PARALLEL PROCESSING ==========
  
  # Determine number of cores (use up to 32)
  n_cores <- min(24, detectCores())
  cat(sprintf("      Using %d cores for parallel processing\n", n_cores))
  
  # Process in chunks to manage memory
  chunk_size <- 1000
  gene_indices <- 1:n_genes
  chunks <- split(gene_indices, ceiling(seq_along(gene_indices) / chunk_size))
  cat(sprintf("      Processing %d chunks\n", length(chunks)))
  
  # Initialize counters
  skip_no_gene_in_cancer <- 0
  skip_no_gene_in_healthy <- 0
  skip_zero_replicates_cancer <- 0
  skip_zero_replicates_healthy <- 0
  skip_no_family <- 0
  skip_invalid_family <- 0
  skip_cancer_resid_pool <- 0
  skip_healthy_resid_pool_own <- 0
  skip_healthy_resid_pool_fam <- 0
  skip_healthy_resid_pool_orth <- 0
  skip_missing_obs <- 0
  skip_errors <- 0
  processed_genes <- 0
  
  for (chunk_idx in seq_along(chunks)) {
    chunk_genes <- chunks[[chunk_idx]]
    cat(sprintf("        Chunk %d/%d: processing %d genes...\n", 
                chunk_idx, length(chunks), length(chunk_genes)))
    
    # Process chunk in parallel
    chunk_results <- mclapply(chunk_genes, process_gene, 
                              mc.cores = n_cores, 
                              mc.preschedule = TRUE)
    
    # Check if results are valid
    if (is.null(chunk_results) || !is.list(chunk_results)) {
      cat("        WARNING: chunk_results is NULL or not a list\n")
      next
    }
    
    # Collect results
    for (local_idx in seq_along(chunk_genes)) {
      g <- chunk_genes[local_idx]
      result <- chunk_results[[local_idx]]
      
      # Skip if result is NULL or not a list
      if (is.null(result) || !is.list(result)) {
        skip_errors <- skip_errors + 1
        next
      }
      
      if (result$success) {
        noise_pvalues[g, "own_healthy"] <- result$pvalues["own"]
        noise_pvalues[g, "family_mean"] <- result$pvalues["fam"]
        noise_pvalues[g, "ortholog_mean"] <- result$pvalues["orth"]
        processed_genes <- processed_genes + 1
      } else {
        # Update skip counters based on reason
        reason <- result$reason
        if (grepl("error:", reason)) {
          skip_errors <- skip_errors + 1
          if (skip_errors <= 5) {
            cat(sprintf("          Error for gene %d: %s\n", g, reason))
          }
        } else if (reason == "no_gene_in_cancer") skip_no_gene_in_cancer <- skip_no_gene_in_cancer + 1
        else if (reason == "no_gene_in_healthy") skip_no_gene_in_healthy <- skip_no_gene_in_healthy + 1
        else if (reason == "zero_replicates_cancer") skip_zero_replicates_cancer <- skip_zero_replicates_cancer + 1
        else if (reason == "zero_replicates_healthy") skip_zero_replicates_healthy <- skip_zero_replicates_healthy + 1
        else if (reason == "no_family") skip_no_family <- skip_no_family + 1
        else if (reason == "invalid_family_data") skip_invalid_family <- skip_invalid_family + 1
        else if (reason == "cancer_resid_pool") skip_cancer_resid_pool <- skip_cancer_resid_pool + 1
        else if (reason == "healthy_resid_pool_own") skip_healthy_resid_pool_own <- skip_healthy_resid_pool_own + 1
        else if (reason == "healthy_resid_pool_fam") skip_healthy_resid_pool_fam <- skip_healthy_resid_pool_fam + 1
        else if (reason == "healthy_resid_pool_orth") skip_healthy_resid_pool_orth <- skip_healthy_resid_pool_orth + 1
        else if (reason == "missing_obs") skip_missing_obs <- skip_missing_obs + 1
        else skip_errors <- skip_errors + 1
      }
    }
    
    # Clean up after each chunk
    gc()
  }
  
  # ========== OUTPUT STATISTICS ==========
  cat(sprintf("      Filter statistics:\n"))
  cat(sprintf("        - No gene in cancer: %d\n", skip_no_gene_in_cancer))
  cat(sprintf("        - No gene in healthy: %d\n", skip_no_gene_in_healthy))
  cat(sprintf("        - Zero replicates (cancer): %d\n", skip_zero_replicates_cancer))
  cat(sprintf("        - Zero replicates (healthy): %d\n", skip_zero_replicates_healthy))
  cat(sprintf("        - No family info: %d\n", skip_no_family))
  cat(sprintf("        - Invalid family data: %d\n", skip_invalid_family))
  cat(sprintf("        - Cancer resid pool <10: %d\n", skip_cancer_resid_pool))
  cat(sprintf("        - Healthy resid pool (own) <10: %d\n", skip_healthy_resid_pool_own))
  cat(sprintf("        - Healthy resid pool (fam) <10: %d\n", skip_healthy_resid_pool_fam))
  cat(sprintf("        - Healthy resid pool (orth) <10: %d\n", skip_healthy_resid_pool_orth))
  cat(sprintf("        - Missing observed differences: %d\n", skip_missing_obs))
  cat(sprintf("        - Errors: %d\n", skip_errors))
  cat(sprintf("        - Successfully processed: %d\n", processed_genes))

  if (!is.null(output_dir)) {
    safe_stage <- gsub(" ", "_", stage)
    result_file <- file.path(output_dir, "intermediate", 
                              sprintf("%s_%s_%s_pvalues.rds", 
                                      cancer_id, norm_method, safe_stage))
    
    # Also create a marker file for completion
    done_file <- file.path(output_dir, "intermediate", "done",
                            sprintf("%s_%s_%s.done", 
                                    cancer_id, norm_method, safe_stage))
    dir.create(dirname(done_file), recursive = TRUE, showWarnings = FALSE)
    
    # Save the p-values
    saveRDS(list(
      pvalues = noise_pvalues,
      statistics = list(
        processed_genes = processed_genes,
        total_genes = n_genes,
        skip_statistics = list(
          no_gene_in_cancer = skip_no_gene_in_cancer,
          no_gene_in_healthy = skip_no_gene_in_healthy,
          zero_replicates_cancer = skip_zero_replicates_cancer,
          zero_replicates_healthy = skip_zero_replicates_healthy,
          no_family = skip_no_family,
          invalid_family = skip_invalid_family,
          cancer_resid_pool = skip_cancer_resid_pool,
          healthy_resid_pool_own = skip_healthy_resid_pool_own,
          healthy_resid_pool_fam = skip_healthy_resid_pool_fam,
          healthy_resid_pool_orth = skip_healthy_resid_pool_orth,
          missing_obs = skip_missing_obs,
          errors = skip_errors
        )
      ),
      timestamp = Sys.time()
    ), result_file)
  }
  
  list(pvalues = noise_pvalues, null_sample = numeric(0))
}

# ==================== GESAMTTABELLE ERSTELLEN (MIT STAGE-WEISER VERARBEITUNG) ====================

create_all_genes_table <- function(all_results, output_dir = NULL) {
  
  cat("\n>>> Erstelle Gesamttabelle (stage-weise) ...\n")
  
  # Pre-allocate Liste
  total_estimated <- 0
  for (res_name in names(all_results)) {
    res <- all_results[[res_name]]
    n_genes <- nrow(res$gene_results)
    total_estimated <- total_estimated + (n_genes * length(STAGES) * length(COMP_TYPES))
  }
  buffer <- 1.1
  alloc_rows <- ceiling(total_estimated * buffer)
  all_genes_list <- vector("list", alloc_rows)
  global_idx <- 1
  
  # Für Null-Distanzen
  null_dist_list <- vector("list", length(all_results) * length(STAGES))
  null_idx <- 1
  
  # Schleife über alle Datensätze
  for (res_name in names(all_results)) {
    res <- all_results[[res_name]]
    cat(sprintf("\n  Verarbeite %s (%s, %s)...\n", 
                res$cancer_type, res$norm_method, res$norm_display))
    
    family_n_genes <- tabulate(res$gene_to_fam, nbins = res$n_families)
    
    # Extrahiere Distance p-Werte (werden nicht mehr benötigt? Wir behalten sie für Vergleich)
    distance_p <- extract_distance_pvalues(res$gene_results)
    
    # Gesunde Preprocessed Daten (bereits vorhanden)
    healthy_preproc <- res$healthy_preproc
    
    # Für jedes Stadium separat rechnen
    for (stage in STAGES) {
      cat(sprintf("    Stadium: %s\n", stage))
      
      # Lade Krebs-Replikate für dieses Stadium (Rohdaten)
      cancer_data_raw <- load_stage_data(
        project_id = res$cancer_id,
        stage = stage,
        data_type = "cancer",
        use_constant_healthy = FALSE,
        norm_method = "raw",
        apply_mean = FALSE,
        normalize = FALSE
      )
      
      if (is.null(cancer_data_raw) || !is.matrix(cancer_data_raw$expression_vectors)) {
        cat("      Keine Krebsdaten für dieses Stadium.\n")
        next
      }
      
      # Preprocess Krebs-Replikate mit der gleichen Normmethode
      cancer_preproc <- preprocess_replicates(cancer_data_raw$expression_vectors, res$norm_method)
      
      # Noise p-Werte für dieses Stadium berechnen
      noise_result <- compute_noise_pvalues_stage(
        cancer_preproc = cancer_preproc,
        healthy_preproc = healthy_preproc,
        gene_results = res$gene_results,
        family_stats = res$family_stats,
        gene_to_fam = res$gene_to_fam,
        stage = stage,
        norm_method = res$norm_method,
        B = 10000,
        k_start = 32, k_step = 16, k_max = 1024, tau = 0.15,
        output_dir = output_dir,
        cancer_id = res$cancer_id
      )
      
      noise_p <- noise_result$pvalues
      
      # Null-Distanzen sammeln (optional)
      if (length(noise_result$null_sample) > 0) {
        null_dist_list[[null_idx]] <- noise_result$null_sample
        null_idx <- null_idx + 1
      }
      
      n_genes <- nrow(res$gene_results)
      rows_added <- 0
      
      # Daten für dieses Stadium vorbereiten (Pre‑allocate Vektoren)
      gene_ids_all <- character(n_genes * length(COMP_TYPES))
      cancer_types_all <- character(n_genes * length(COMP_TYPES))
      cancer_ids_all <- character(n_genes * length(COMP_TYPES))
      stages_all <- character(n_genes * length(COMP_TYPES))
      normalizations_all <- character(n_genes * length(COMP_TYPES))
      norm_methods_all <- character(n_genes * length(COMP_TYPES))
      comparisons_all <- character(n_genes * length(COMP_TYPES))
      directions_all <- character(n_genes * length(COMP_TYPES))
      signed_diffs_all <- numeric(n_genes * length(COMP_TYPES))
      dist_p_all <- numeric(n_genes * length(COMP_TYPES))
      noise_p_all <- numeric(n_genes * length(COMP_TYPES))
      family_sizes_all <- integer(n_genes * length(COMP_TYPES))
      family_means_all <- numeric(n_genes * length(COMP_TYPES))
      
      local_idx <- 1
      
      for (g in 1:n_genes) {
        gene_id <- res$gene_results$gene_id[g]
        for (comp_type in COMP_TYPES) {
          dist_p_val <- distance_p[g, stage, comp_type]
          noise_p_val <- noise_p[g, comp_type]
          
          if (is.na(dist_p_val) || is.na(noise_p_val)) next
          
          shift_col <- switch(comp_type,
                              "own_healthy" = paste0("shift_vs_own_healthy_", gsub(" ", "_", stage)),
                              "family_mean" = paste0("shift_vs_family_mean_", gsub(" ", "_", stage)),
                              "ortholog_mean" = paste0("shift_vs_ortholog_mean_", gsub(" ", "_", stage)))
          
          if (!shift_col %in% colnames(res$gene_results)) next
          
          shift <- res$gene_results[g, shift_col]
          direction <- ifelse(is.na(shift), NA, ifelse(shift > 0, "up", "down"))
          
          fam_id <- res$gene_to_fam[g]
          family_size <- if (!is.na(fam_id) && fam_id > 0) family_n_genes[fam_id] else NA
          family_mean_val <- if (!is.na(fam_id) && fam_id > 0) 
            res$family_stats$gene_family_means_all[fam_id] else NA
          
          gene_ids_all[local_idx] <- gene_id
          cancer_types_all[local_idx] <- res$cancer_type
          cancer_ids_all[local_idx] <- res$cancer_id
          stages_all[local_idx] <- stage
          normalizations_all[local_idx] <- res$norm_display
          norm_methods_all[local_idx] <- res$norm_method
          comparisons_all[local_idx] <- comp_type
          directions_all[local_idx] <- direction
          signed_diffs_all[local_idx] <- shift
          dist_p_all[local_idx] <- dist_p_val
          noise_p_all[local_idx] <- noise_p_val
          family_sizes_all[local_idx] <- family_size
          family_means_all[local_idx] <- family_mean_val
          
          local_idx <- local_idx + 1
          rows_added <- rows_added + 1
        }
      }
      
      if (rows_added > 0) {
        stage_df <- data.frame(
          gene_id = gene_ids_all[1:rows_added],
          cancer_type = cancer_types_all[1:rows_added],
          cancer_id = cancer_ids_all[1:rows_added],
          stage = stages_all[1:rows_added],
          normalization = normalizations_all[1:rows_added],
          norm_method = norm_methods_all[1:rows_added],
          comparison = comparisons_all[1:rows_added],
          direction = directions_all[1:rows_added],
          signed_difference = signed_diffs_all[1:rows_added],
          distance_p_value = dist_p_all[1:rows_added],
          noise_p_value = noise_p_all[1:rows_added],
          family_size = family_sizes_all[1:rows_added],
          family_mean = family_means_all[1:rows_added],
          stringsAsFactors = FALSE
        )
        
        all_genes_list[[global_idx]] <- stage_df
        global_idx <- global_idx + 1
      }
      
      cat(sprintf("      → %d Zeilen für Stadium %s\n", rows_added, stage))
    }
  }
  
  cat(sprintf("\n>>> Gesamt: %d Zeilen generiert\n", global_idx - 1))
  
  # Kombiniere alle Teile
  all_genes_list <- all_genes_list[1:(global_idx - 1)]
  result <- rbindlist(all_genes_list)
  
  # Null-Distanzen speichern (optional)
  null_dist_list <- null_dist_list[1:(null_idx - 1)]
  all_null_distances <- unlist(null_dist_list)
  null_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "diagnostic_plots", "null_distances")
  dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(all_null_distances, file.path(null_dir, "global_null_distances.rds"))
  
  cat(sprintf("\n FINALE TABELLE: %d Zeilen\n", nrow(result)))
  return(result)
}

# ==================== PLOT FUNKTIONEN (unverändert) ====================

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

analyze_outlier_significance <- function(output_dir = NULL) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("OUTLIER SIGNIFICANCE ANALYSIS (stage‑wise adaptive kNN, B=10000)\n")
  cat("Parameters: k_start=32, k_step=16, k_max=1024, tau=0.15\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  all_results <- load_all_gene_results()
  if (length(all_results) == 0) stop("Keine Daten gefunden!")
  
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===== 1. GESAMTTABELLE =====
  cat("\n>>> Erstelle Gesamttabelle...\n")
  all_genes_table <- create_all_genes_table(all_results, output_dir = output_dir)
  
  # ===== 2. SPEICHERE ROHDATEN =====
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
  cat("SELEKTION: noise_p_adj < 0.05 (sortiert nach distance_p_value)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  ALPHA_NOISE <- 0.05
  
  significant_by_norm <- list()
  
  for (current_norm in unique(all_genes_table$norm_method)) {
    
    cat(sprintf("\n>>> Verarbeite Normalisierung: %s <<<\n", current_norm))
    
    norm_data <- all_genes_table %>% 
      filter(norm_method == current_norm)
    
    # BH-Korrektur für noise_p_value (pro cancer_id, norm, stage, comp)
    norm_data <- norm_data %>%
      group_by(cancer_id, normalization, stage, comparison) %>%
      mutate(
        noise_p_value_adj = p.adjust(noise_p_value, method = "BH")
      ) %>%
      ungroup()
    
    # ===== DIAGNOSE PRO VERGLEICHSTYP =====
    cat(sprintf("\n  DIAGNOSE %s:\n", current_norm))
    
    for (comp_type in COMP_TYPES) {
      comp_data <- norm_data %>% filter(comparison == comp_type)
      
      cat(sprintf("\n    --- %s ---\n", comp_type))
      cat(sprintf("      Noise p (unkorr) < 0.05: %d (%.1f%%)\n", 
                  sum(comp_data$noise_p_value < 0.05, na.rm = TRUE),
                  100 * mean(comp_data$noise_p_value < 0.05, na.rm = TRUE)))
      cat(sprintf("      Noise p (adj) < 0.05: %d (%.1f%%)\n\n", 
                  sum(comp_data$noise_p_value_adj < 0.05, na.rm = TRUE),
                  100 * mean(comp_data$noise_p_value_adj < 0.05, na.rm = TRUE)))
      cat("      Summary of noise_p_value: \n")
      print(summary(comp_data$noise_p_value))
      cat("      Summary of noise_p_value_adj: \n")
      print(summary(comp_data$noise_p_value_adj))
    }
    
    # ===== SELEKTION: Nur noise_p_adj < 0.05 =====
    significant <- norm_data %>%
      filter(noise_p_value_adj < ALPHA_NOISE) %>%
      arrange(distance_p_value)
    
    cat(sprintf("\n  → Signifikante Einträge (noise_p_adj < 0.05): %d\n", nrow(significant)))
    
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
    
    cat(sprintf("\n  Gesamt: %d signifikante Einträge (noise_p_adj < 0.05)\n", nrow(all_significant)))
    
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
  cat("\n>>> Jaccard Plots (basierend auf noise_p_adj < 0.05)\n")
  
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
            subtitle = sprintf("noise_p_adj < 0.05, n=%d", nrow(plot_data))
          )
          
          filename_adj <- sprintf("jaccard_ADJ_%s_min%d.png", comp_type, min_size)
          ggsave(file.path(plots_dir, filename_adj), p_adj, 
                 width = 7, height = 6, dpi = 300)
          
          # Rohdaten für Vergleich
          raw_data <- all_genes_table %>%
            filter(
              comparison == comp_type,
              family_size >= min_size,
              distance_p_value < 0.05,
              noise_p_value < 0.05
            )
          
          if (nrow(raw_data) > 0 && nrow(plot_data) > 0) {
            norm_methods <- c("raw", "log", "std_log")
            
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
            
            comparison_mat <- matrix(NA, 3, 3)
            rownames(comparison_mat) <- norm_methods
            colnames(comparison_mat) <- norm_methods
            
            for (i in 1:3) {
              for (j in 1:3) {
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
              labs(title = sprintf("%s (min %d) - Rohdaten vs. noise_p_adj < 0.05", 
                                   comp_type, min_size),
                   x = "noise_p_adj < 0.05", y = "Rohdaten (dist_p<0.05 & noise_p<0.05)") +
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
  base_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_stagewise")
  
  cat("\n", paste(rep("#", 80), collapse = ""), "\n")
  cat("### ADAPTIVE NOISE MODEL ANALYSIS (STAGE‑WISE) ###\n")
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