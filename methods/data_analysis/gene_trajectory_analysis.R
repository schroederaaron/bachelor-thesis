library(dplyr)

source("config.R")
source("utils.R")

#' Analyze gene trajectories for one normalization method
#' @param project_id Project ID
#' @param ortholog_proteins Vector of ortholog protein IDs
#' @param norm_method Normalization method: "raw", "full", "std_log", or "log"
#' @param use_constant_healthy Whether to use constant healthy reference (default TRUE)
#' @return List with gene analysis results
analyze_gene_trajectories_norm <- function(project_id, ortholog_proteins, 
                                           norm_method = "raw",
                                           use_constant_healthy = TRUE) {
  
  norm_display <- get_norm_display(norm_method)
  norm_suffix <- get_norm_suffix(norm_method)
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("GENE TRAJECTORY ANALYSIS: %s [%s]\n", project_id, norm_display))
  cat(sprintf("%s\n\n", paste(rep("=", 60), collapse = "")))
  
  # ==================== LOAD FAMILY LOESS STATS ====================
  # These contain the smoothed within-family gene SDs needed for normalization
  
  family_output_dir <- get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method)
  family_loess_file <- file.path(family_output_dir, project_id, 
                                paste0(project_id, "_family_loess_stats", norm_suffix, ".rds"))
  
  if (!file.exists(family_loess_file)) {
    cat(sprintf("  Error: Family LOESS stats not found for [%s].\n", norm_display))
    cat(sprintf("  Expected file: %s\n", family_loess_file))
    cat(sprintf("  Run family analysis with the SAME normalization method first!\n"))
    return(NULL)
  }
  
  # Load family LOESS stats
  family_stats <- readRDS(family_loess_file)
  
  # Extract the smoothed within-family gene SDs
  # These are what we use to normalize gene-level shifts
  gene_smoothed_sds_all <- family_stats$gene_smoothed_sds_all
  gene_smoothed_sds_ortholog <- family_stats$gene_smoothed_sds_ortholog
  gene_family_means_all <- family_stats$gene_family_means_all
  gene_family_means_ortholog <- family_stats$gene_family_means_ortholog
  gene_to_fam_ref <- family_stats$gene_to_fam
  
  cat(sprintf("  Loaded family LOESS stats: %d families\n", 
              length(gene_smoothed_sds_all)))
  cat(sprintf("  Families with valid smoothed SDs (all genes): %d\n", 
              sum(!is.na(gene_smoothed_sds_all))))
  cat(sprintf("  Families with valid smoothed SDs (orthologs): %d\n", 
              sum(!is.na(gene_smoothed_sds_ortholog))))
  
  # ==================== CREATE OUTPUT DIRECTORIES ====================
  
  output_dir <- get_norm_output_dir(GENE_OUTPUT_DIR, norm_method)
  proj_dir <- file.path(output_dir, project_id)
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ==================== LOAD CONSTANT HEALTHY REFERENCE ====================
  
  healthy_constant <- NULL
  healthy_constant_means <- NULL
  healthy_constant_gene_ids <- NULL
  
  if (use_constant_healthy) {
    cat("  Loading constant healthy reference...\n")
    healthy_constant_data <- load_stage_data(project_id, STAGES[1], "healthy", 
                                            use_constant_healthy = TRUE, 
                                            norm_method = norm_method)
    if (!is.null(healthy_constant_data)) {
      # With our updated load_stage_data, expression_vectors is now a 1-row matrix
      # containing gene means
      if (is.matrix(healthy_constant_data$expression_vectors)) {
        healthy_constant_means <- as.vector(healthy_constant_data$expression_vectors)
        names(healthy_constant_means) <- colnames(healthy_constant_data$expression_vectors)
      } else {
        healthy_constant_means <- healthy_constant_data$expression_vectors
      }
      healthy_constant_gene_ids <- healthy_constant_data$gene_ids
      cat(sprintf("  ✓ Loaded constant healthy reference (%d genes) [%s]\n", 
                  length(healthy_constant_means), norm_display))
    } else {
      cat("  Warning: Could not load constant healthy reference. Using stage-specific healthy data.\n")
      use_constant_healthy <- FALSE
    }
  }
  
  # ==================== INITIALIZE STORAGE ====================
  
  gene_results_all <- list()
  
  # Store stage-specific shifts for each comparison type
  stage_shifts_gene <- list(
    gene_vs_own_healthy = list(),      # cancer gene vs its own healthy counterpart
    gene_vs_family_mean = list(),      # cancer gene vs healthy family mean (all genes)
    gene_vs_family_ortholog_mean = list()  # cancer gene vs healthy family ortholog mean
  )
  
  # For storing gene-level data across stages
  ref_genes <- NULL
  n_genes <- 0
  n_families <- length(gene_smoothed_sds_all)
  
  # ==================== PROCESS EACH STAGE ====================
  
  for (stage in STAGES) {
    cat(sprintf("\n  Processing %s...\n", stage))
    
    # Load cancer data with same normalization
    cancer_data <- load_stage_data(project_id, stage, "cancer", 
                                  use_constant_healthy = FALSE, 
                                  norm_method = norm_method)
    if (is.null(cancer_data)) {
      cat(sprintf("    Warning: Missing cancer data for %s\n", stage))
      next
    }
    
    # Extract cancer gene means (now a vector from the 1-row matrix)
    if (is.matrix(cancer_data$expression_vectors)) {
      cancer_gene_means <- as.vector(cancer_data$expression_vectors)
      names(cancer_gene_means) <- colnames(cancer_data$expression_vectors)
    } else {
      cancer_gene_means <- cancer_data$expression_vectors
    }
    gene_ids <- cancer_data$gene_ids
    gene_to_fam <- cancer_data$gene_to_fam
    
    # Verify gene order matches reference
    if (is.null(ref_genes)) {
      ref_genes <- gene_ids
      n_genes <- length(ref_genes)
    }
    
    # Get healthy gene means (either constant or stage-specific)
    if (use_constant_healthy && !is.null(healthy_constant_means)) {
      # Use constant healthy reference, aligning genes to current order
      healthy_gene_means <- healthy_constant_means[gene_ids]
      cat("    Using constant healthy reference\n")
    } else {
      # Load stage-specific healthy data
      healthy_data <- load_stage_data(project_id, stage, "healthy", 
                                     use_constant_healthy = FALSE, 
                                     norm_method = norm_method)
      if (is.null(healthy_data)) {
        cat(sprintf("    Warning: Missing healthy data for %s\n", stage))
        next
      }
      
      # Extract healthy gene means
      if (is.matrix(healthy_data$expression_vectors)) {
        healthy_gene_means <- as.vector(healthy_data$expression_vectors)
        names(healthy_gene_means) <- colnames(healthy_data$expression_vectors)
      } else {
        healthy_gene_means <- healthy_data$expression_vectors
      }
    }
    
    # Identify ortholog genes
    is_ortholog <- gene_ids %in% ortholog_proteins
        
    # ==================== COMPUTE GENE SHIFTS ====================
    # Using the smoothed within-family gene SDs from family analysis
    
    cat(sprintf("    Computing shifts for %d genes...\n", length(gene_ids)))
    
    # Initialize shift vectors
    shifts_own <- rep(NA, length(gene_ids))
    shifts_fam <- rep(NA, length(gene_ids))
    shifts_orth <- rep(NA, length(gene_ids))
    shifts_own_unscaled <- rep(NA, length(gene_ids))
    shifts_fam_unscaled <- rep(NA, length(gene_ids))
    shifts_orth_unscaled <- rep(NA, length(gene_ids))
    
    # For tracking which families have valid SDs
    valid_sd_all <- !is.na(gene_smoothed_sds_all) & gene_smoothed_sds_all > 0
    valid_sd_orth <- !is.na(gene_smoothed_sds_ortholog) & gene_smoothed_sds_ortholog > 0
    
    cat(sprintf("    Families with valid SD (all genes): %d\n", sum(valid_sd_all)))
    cat(sprintf("    Families with valid SD (orthologs): %d\n", sum(valid_sd_orth)))
    
    for (i in 1:length(gene_ids)) {
      fam_id <- gene_to_fam[i]
      
      # 1. Gene vs its own healthy counterpart
      # Normalized by within-family gene SD (smoothed, all genes)
      if (!is.na(cancer_gene_means[i]) && !is.na(healthy_gene_means[i]) &&
          fam_id <= length(gene_smoothed_sds_all) &&
          !is.na(gene_smoothed_sds_all[fam_id]) && gene_smoothed_sds_all[fam_id] > 0) {
        shifts_own[i] <- (cancer_gene_means[i] - healthy_gene_means[i]) / gene_smoothed_sds_all[fam_id]
        shifts_own_unscaled[i] <- cancer_gene_means[i] - healthy_gene_means[i]
      }
      
      # 2. Gene vs healthy family mean (all genes)
      # Normalized by within-family gene SD (smoothed, all genes)
      if (!is.na(cancer_gene_means[i]) && 
          fam_id <= length(gene_family_means_all) && !is.na(gene_family_means_all[fam_id]) &&
          fam_id <= length(gene_smoothed_sds_all) &&
          !is.na(gene_smoothed_sds_all[fam_id]) && gene_smoothed_sds_all[fam_id] > 0) {
        shifts_fam[i] <- (cancer_gene_means[i] - gene_family_means_all[fam_id]) / gene_smoothed_sds_all[fam_id]
        shifts_fam_unscaled[i] <- cancer_gene_means[i] - gene_family_means_all[fam_id]
      }
      
      # 3. Gene vs healthy family ortholog mean
      # Normalized by within-family ortholog gene SD (smoothed)
      if (!is.na(cancer_gene_means[i]) &&
          fam_id <= length(gene_family_means_ortholog) && !is.na(gene_family_means_ortholog[fam_id]) &&
          fam_id <= length(gene_smoothed_sds_ortholog) &&
          !is.na(gene_smoothed_sds_ortholog[fam_id]) && gene_smoothed_sds_ortholog[fam_id] > 0) {
        shifts_orth[i] <- (cancer_gene_means[i] - gene_family_means_ortholog[fam_id]) / gene_smoothed_sds_ortholog[fam_id]
        shifts_orth_unscaled[i] <- cancer_gene_means[i] - gene_family_means_ortholog[fam_id]
      }
    }
    
    # Store stage shifts
    stage_shifts_gene$gene_vs_own_healthy[[stage]] <- shifts_own
    stage_shifts_gene$gene_vs_own_healthy_unscaled[[stage]] <- shifts_own_unscaled
    stage_shifts_gene$gene_vs_family_mean[[stage]] <- shifts_fam
    stage_shifts_gene$gene_vs_family_mean_unscaled[[stage]] <- shifts_fam_unscaled
    stage_shifts_gene$gene_vs_family_ortholog_mean[[stage]] <- shifts_orth
    stage_shifts_gene$gene_vs_family_ortholog_mean_unscaled[[stage]] <- shifts_orth_unscaled

    # Create detailed comparison list for saving
    comparisons <- vector("list", length(gene_ids))
    for (i in 1:length(gene_ids)) {
      comparisons[[i]] <- list(
        gene_id = gene_ids[i],
        family_id = gene_to_fam[i],
        is_ortholog = is_ortholog[i],
        shift_gene_vs_own_healthy = shifts_own[i],
        shift_gene_vs_own_healthy_unscaled = shifts_own_unscaled[i],
        shift_gene_vs_family_mean = shifts_fam[i],
        shift_gene_vs_family_mean_unscaled = shifts_fam_unscaled[i],
        shift_gene_vs_family_ortholog_mean = shifts_orth[i],
        shift_gene_vs_family_ortholog_mean_unscaled = shifts_orth_unscaled[i],
        cancer_expression = cancer_gene_means[i],
        healthy_expression = healthy_gene_means[i]
      )
    }
    
    # Store stage results
    stage_result <- list(
      stage = stage,
      norm_method = norm_method,
      norm_display = norm_display,
      comparisons = comparisons,
      n_genes = length(gene_ids),
      n_families = n_families,
      gene_ids = gene_ids,
      gene_to_fam = gene_to_fam,
      is_ortholog = is_ortholog,
      shifts_own = shifts_own,
      shifts_own_unscaled = shifts_own_unscaled,
      shifts_fam = shifts_fam,
      shifts_fam_unscaled = shifts_fam_unscaled,
      shifts_orth = shifts_orth,
      shifts_orth_unscaled = shifts_orth_unscaled,
      gene_smoothed_sds_all = gene_smoothed_sds_all,
      gene_smoothed_sds_ortholog = gene_smoothed_sds_ortholog
    )
    
    gene_results_all[[stage]] <- stage_result
    
    # Save intermediate results
    save_intermediate_results(stage_result, project_id, 
                             paste0("gene_results_", gsub(" ", "_", stage)), 
                             GENE_OUTPUT_DIR, norm_method)
  }
  
  # Check if we have data for at least 2 stages
  if (length(gene_results_all) < 2) {
    cat("  Insufficient data for gene trajectory analysis\n")
    return(NULL)
  }
  
  # ==================== CREATE STAGE SHIFT MATRICES ====================
  
  cat("\n  Creating stage shift matrices...\n")
  
  # Initialize matrices for each comparison type
  stage_shifts_matrices <- list(
    gene_vs_own_healthy = matrix(NA, nrow = n_genes, ncol = length(STAGES)),
    gene_vs_family_mean = matrix(NA, nrow = n_genes, ncol = length(STAGES)),
    gene_vs_family_ortholog_mean = matrix(NA, nrow = n_genes, ncol = length(STAGES))
  )
  
  # Set column names
  for (comp_type in names(stage_shifts_matrices)) {
    colnames(stage_shifts_matrices[[comp_type]]) <- STAGES
    rownames(stage_shifts_matrices[[comp_type]]) <- ref_genes
  }
  
  # Fill matrices
  for (stage in STAGES) {
    if (stage %in% names(gene_results_all)) {
      stage_result <- gene_results_all[[stage]]
      
      # Match genes to reference ordering
      gene_match <- match(ref_genes, stage_result$gene_ids)
      
      for (comp_type in names(stage_shifts_matrices)) {
        # Extract the appropriate shift vector
        shift_vec <- switch(comp_type,
          gene_vs_own_healthy = stage_result$shifts_own,
          gene_vs_family_mean = stage_result$shifts_fam,
          gene_vs_family_ortholog_mean = stage_result$shifts_orth
        )
        
        # Place in matrix with correct gene order
        if (!is.null(shift_vec) && length(shift_vec) >= max(gene_match, na.rm = TRUE)) {
          stage_shifts_matrices[[comp_type]][, stage] <- shift_vec[gene_match]
        }
      }
    }
  }

  # After building stage_shifts_matrices, add:
  cat("\n=== Stage shift matrix valid counts ===\n")
  for (comp_type in names(stage_shifts_matrices)) {
    cat("\n", comp_type, ":\n")
    for (stage in STAGES) {
      n_valid <- sum(!is.na(stage_shifts_matrices[[comp_type]][, stage]))
      cat("  ", stage, ": ", n_valid, " valid genes\n")
    }
  }
  cat("======================================\n")
  
  # ==================== COMPUTE GENE TRAJECTORY METRICS ====================
  
  cat("  Computing gene trajectory metrics...\n")
  
  # Initialize results data frame
  gene_trajectory_results <- data.frame(
    gene_id = ref_genes,
    family_id = gene_to_fam_ref[match(ref_genes, names(gene_to_fam_ref))],
    is_ortholog = gene_results_all[[names(gene_results_all)[1]]]$is_ortholog[
      match(ref_genes, gene_results_all[[names(gene_results_all)[1]]]$gene_ids)
    ],
    
    # Total distances (sum of absolute shifts)
    total_distance_gene_vs_own = numeric(n_genes),
    total_distance_gene_vs_fam = numeric(n_genes),
    total_distance_gene_vs_orth = numeric(n_genes),
    
    # Scaled total distances (will be computed later)
    total_distance_scaled_own = numeric(n_genes),
    total_distance_scaled_fam = numeric(n_genes),
    total_distance_scaled_orth = numeric(n_genes),
    
    # Outlier flags
    is_outlier_own = logical(n_genes),
    is_outlier_fam = logical(n_genes),
    is_outlier_orth = logical(n_genes),
    
    # Number of stages with valid data
    n_stages_own = integer(n_genes),
    n_stages_fam = integer(n_genes),
    n_stages_orth = integer(n_genes),
    
    stringsAsFactors = FALSE
  )
  
  # Compute total distances (sum of absolute shifts)
  for (i in 1:n_genes) {
    for (comp_type in c("own", "fam", "orth")) {
      matrix_name <- switch(comp_type,
        own = "gene_vs_own_healthy",
        fam = "gene_vs_family_mean",
        orth = "gene_vs_family_ortholog_mean"
      )
      
      shifts <- stage_shifts_matrices[[matrix_name]][i, ]
      valid_shifts <- shifts[!is.na(shifts)]
      
      if (length(valid_shifts) > 0) {
        gene_trajectory_results[i, paste0("total_distance_gene_vs_", comp_type)] <- sum(abs(valid_shifts))
        gene_trajectory_results[i, paste0("n_stages_", comp_type)] <- length(valid_shifts)
      }
    }
  }
  
  # ==================== SCALED TOTAL DISTANCE OUTLIERS ====================
  
  cat("  Computing scaled total distance outliers...\n")

  # Gene vs own
  dist_col <- "total_distance_gene_vs_own"
  distances <- gene_trajectory_results[[dist_col]]

  outlier_result <- compute_scaled_total_outliers(
    total_distances = distances,
    percentile = OUTLIER_PERCENTILE
  )

  gene_trajectory_results$total_distance_scaled_own <- outlier_result$scaled_distances
  gene_trajectory_results$is_outlier_own <- outlier_result$is_outlier
  gene_trajectory_results$p_value_total_own <- outlier_result$p_values

  # Gene vs family mean
  dist_col <- "total_distance_gene_vs_fam"
  distances <- gene_trajectory_results[[dist_col]]

  outlier_result <- compute_scaled_total_outliers(
    total_distances = distances,
    percentile = OUTLIER_PERCENTILE
  )

  gene_trajectory_results$total_distance_scaled_fam <- outlier_result$scaled_distances
  gene_trajectory_results$is_outlier_fam <- outlier_result$is_outlier
  gene_trajectory_results$p_value_total_fam <- outlier_result$p_values 

  # Gene vs ortholog mean
  dist_col <- "total_distance_gene_vs_orth"
  distances <- gene_trajectory_results[[dist_col]]

  outlier_result <- compute_scaled_total_outliers(
    total_distances = distances,
    percentile = OUTLIER_PERCENTILE
  )

  gene_trajectory_results$total_distance_scaled_orth <- outlier_result$scaled_distances
  gene_trajectory_results$is_outlier_orth <- outlier_result$is_outlier
  gene_trajectory_results$p_value_total_orth <- outlier_result$p_values
  
  # ==================== GENE VELOCITIES ====================
  
  cat("  Computing gene velocities and detecting outliers...\n")

  gene_velocities <- list()
  gene_velocity_outliers <- list()
  gene_velocity_pvalues <- list()

  for (comp_type in names(stage_shifts_matrices)) {
    vel_result <- compute_velocity_outliers(
      shifts_matrix = stage_shifts_matrices[[comp_type]],
      stages = STAGES,
      percentile = OUTLIER_PERCENTILE
    )
    
    gene_velocities[[comp_type]] <- vel_result$velocities
    gene_velocity_outliers[[comp_type]] <- vel_result$is_outlier
    gene_velocity_pvalues[[comp_type]] <- vel_result$p_values
  }
  
  # ==================== STAGE-SPECIFIC GENE OUTLIERS ====================
  
  cat("  Detecting stage-specific gene outliers...\n")

  stage_outliers_gene <- list()

  for (comp_type in names(stage_shifts_matrices)) {
    stage_outliers_gene[[comp_type]] <- list()
    
    # The distances are already normalized using family SDs
    cat(sprintf("    Processing %s\n", comp_type))
    
    for (stage in STAGES) {
      # Get absolute shifts for this stage
      stage_shifts <- abs(stage_shifts_matrices[[comp_type]][, stage])
      valid_idx <- !is.na(stage_shifts)
      n_valid <- sum(valid_idx)
      
      cat(sprintf("      %s: %d valid genes\n", stage, n_valid))
      
      if (n_valid >= 10) {
        
        # Identify outliers using the distance values
        outlier_result <- tox_identify_outliers(
          rdi = stage_shifts[valid_idx],
          percentile = OUTLIER_PERCENTILE
        )
        
        # Map back to full gene vector
        is_outlier <- rep(FALSE, n_genes)
        is_outlier[valid_idx] <- outlier_result$is_outlier
        # Store results
        stage_outliers_gene[[comp_type]][[stage]] <- list(
          is_outlier = is_outlier,
          threshold = outlier_result$threshold,
          rdi_values = stage_shifts[valid_idx],  # Store RDI values for reference
          n_valid = n_valid,
          method = "rdi_based",
          p_values = outlier_result$p_values,
          valid_idx = valid_idx
        )
        
        cat(sprintf("        Found %d outliers (%.1f%%)\n", 
                    sum(is_outlier), sum(is_outlier)/n_valid*100))
        
      } else {
        # Not enough valid genes
        cat(sprintf("        WARNING: Only %d valid genes, skipping\n", n_valid))
        stage_outliers_gene[[comp_type]][[stage]] <- list(
          is_outlier = rep(FALSE, n_genes),
          threshold = NA,
          rdi_values = rep(NA, n_genes),
          n_valid = n_valid,
          method = "insufficient_data",
          p_values = rep(NA, n_genes),
          valid_idx = valid_idx
        )
      }
    }
  }
  
  # ==================== SAVE GENE RESULTS ====================

  # ===== 1. PER-STAGE SHIFT SPALTEN HINZUFÜGEN =====
  cat("  Adding per-stage shift columns to results...\n")

  for (stage in STAGES) {
    stage_col <- gsub(" ", "_", stage)
    
    # Family mean shifts
    gene_trajectory_results[[paste0("shift_vs_family_mean_", stage_col)]] <- 
      stage_shifts_matrices[["gene_vs_family_mean"]][, stage]
    gene_trajectory_results[[paste0("shift_vs_family_mean_unscaled_", stage_col)]] <-
      gene_results_all[[stage]]$shifts_fam_unscaled
    
    # Ortholog mean shifts
    gene_trajectory_results[[paste0("shift_vs_ortholog_mean_", stage_col)]] <- 
      stage_shifts_matrices[["gene_vs_family_ortholog_mean"]][, stage]
    gene_trajectory_results[[paste0("shift_vs_ortholog_mean_unscaled_", stage_col)]] <-
      gene_results_all[[stage]]$shifts_orth_unscaled

    # Own healthy shifts
    gene_trajectory_results[[paste0("shift_vs_own_healthy_", stage_col)]] <- 
      stage_shifts_matrices[["gene_vs_own_healthy"]][, stage]
    gene_trajectory_results[[paste0("shift_vs_own_healthy_unscaled_", stage_col)]] <-
      gene_results_all[[stage]]$shifts_own_unscaled
  }

  # ===== 2. DISTANCE P-WERTE AUS STAGE_OUTLIERS EXTRAHIEREN =====
  cat("  Adding distance p-values to results...\n")

  # Initialisiere p-Wert Spalten mit NA
  for (stage in STAGES) {
    stage_col <- gsub(" ", "_", stage)
    gene_trajectory_results[[paste0("p_value_distance_own_", stage_col)]] <- NA
    gene_trajectory_results[[paste0("p_value_distance_fam_", stage_col)]] <- NA
    gene_trajectory_results[[paste0("p_value_distance_orth_", stage_col)]] <- NA
  }

  # Fülle aus stage_outliers_gene
  if (!is.null(stage_outliers_gene)) {
    for (stage in STAGES) {
      stage_col <- gsub(" ", "_", stage)
      
      # Family mean
      if (!is.null(stage_outliers_gene[["gene_vs_family_mean"]][[stage]]$p_values)) {
        p_vals <- stage_outliers_gene[["gene_vs_family_mean"]][[stage]]$p_values
        valid_idx <- stage_outliers_gene[["gene_vs_family_mean"]][[stage]]$valid_idx
        p_full <- rep(NA, n_genes)
        p_full[valid_idx] <- p_vals
        gene_trajectory_results[[paste0("p_value_distance_fam_", stage_col)]] <- p_full
      }
      
      # Ortholog mean
      if (!is.null(stage_outliers_gene[["gene_vs_family_ortholog_mean"]][[stage]]$p_values)) {
        p_vals <- stage_outliers_gene[["gene_vs_family_ortholog_mean"]][[stage]]$p_values
        valid_idx <- stage_outliers_gene[["gene_vs_family_ortholog_mean"]][[stage]]$valid_idx
        p_full <- rep(NA, n_genes)
        p_full[valid_idx] <- p_vals
        gene_trajectory_results[[paste0("p_value_distance_orth_", stage_col)]] <- p_full
      }

      if (!is.null(stage_outliers_gene[["gene_vs_own_healthy"]][[stage]]$p_values)) {
        p_vals <- stage_outliers_gene[["gene_vs_own_healthy"]][[stage]]$p_values
        valid_idx <- stage_outliers_gene[["gene_vs_own_healthy"]][[stage]]$valid_idx
        p_full <- rep(NA, n_genes)
        p_full[valid_idx] <- p_vals
        gene_trajectory_results[[paste0("p_value_distance_own_", stage_col)]] <- p_full
      }
    }
  }

  # ===== 3. FAMILY INFORMATION HINZUFÜGEN =====
  cat("  Adding family information to results...\n")

  # Family means (für ME-log-P plots)
  gene_trajectory_results$family_mean <- family_stats$gene_family_means_all[gene_trajectory_results$family_id]

  # Family sizes (für Filter)
  gene_trajectory_results$family_n_genes <- family_stats$family_n_genes_all[gene_trajectory_results$family_id]

  cat("  Saving gene trajectory results...\n")
  
  # Add normalization info
  gene_trajectory_results$norm_method <- norm_method
  gene_trajectory_results$norm_display <- norm_display
  
  # Save complete gene trajectory results
  save_final_results(gene_trajectory_results, project_id, "gene_results", 
                    GENE_OUTPUT_DIR, norm_method)
  
  # Save stage shifts
  stage_shifts_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    shifts = stage_shifts_matrices,
    stages = STAGES,
    gene_ids = ref_genes,
    gene_to_fam = gene_trajectory_results$family_id,
    comparison_types = names(stage_shifts_matrices)
  )
  
  save_final_results(stage_shifts_data, project_id, "gene_shifts", 
                    GENE_OUTPUT_DIR, norm_method)
  
  # Save velocities
  velocities_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    velocities = gene_velocities,
    outliers = gene_velocity_outliers,
    stages = STAGES,
    comparison_types = names(gene_velocities)
  )
  
  save_final_results(velocities_data, project_id, "gene_velocities", 
                    GENE_OUTPUT_DIR, norm_method)
  
  # Save outlier information
  outliers_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    total_distance_outliers = list(
      gene_vs_own = gene_trajectory_results[, c("gene_id", "total_distance_scaled_own", "is_outlier_own", "p_value_total_own")],
      gene_vs_fam = gene_trajectory_results[, c("gene_id", "total_distance_scaled_fam", "is_outlier_fam", "p_value_total_fam")],
      gene_vs_orth = gene_trajectory_results[, c("gene_id", "total_distance_scaled_orth", "is_outlier_orth", "p_value_total_orth")]
    ),
    stage_specific_outliers = stage_outliers_gene,
    velocity_outliers = gene_velocity_outliers,
    velocity_pvalues = gene_velocity_pvalues,
    family_sds = list(
      all_genes = gene_smoothed_sds_all,
      orthologs = gene_smoothed_sds_ortholog
    )
  )
  
  save_final_results(outliers_data, project_id, "gene_outliers", 
                    GENE_OUTPUT_DIR, norm_method)
  
  # ==================== CREATE GENE SUMMARY ====================
  
  cat("  Creating gene analysis summary...\n")
  
  n_ortholog_genes <- sum(gene_trajectory_results$is_ortholog, na.rm = TRUE)
  
  summary_info <- list(
    project_id = project_id,
    norm_method = norm_method,
    norm_display = norm_display,
    n_genes = n_genes,
    n_ortholog_genes = n_ortholog_genes,
    n_stages_analyzed = length(gene_results_all),
    use_constant_healthy = use_constant_healthy
  )
  
  comp_types <- c("own", "fam", "orth")
  comp_labels <- c("Gene vs Own", "Gene vs Family Mean", "Gene vs Ortholog Mean")
  
  for (i in seq_along(comp_types)) {
    comp_type <- comp_types[i]
    outlier_col <- paste0("is_outlier_", comp_type)
    
    if (outlier_col %in% names(gene_trajectory_results)) {
      n_outliers <- sum(gene_trajectory_results[[outlier_col]], na.rm = TRUE)
      summary_info[[paste0("n_outliers_", comp_type)]] <- n_outliers
      summary_info[[paste0("percent_outliers_", comp_type)]] <- 
        round(n_outliers / n_genes * 100, 1)
    }
  }
  
  # Count velocity outliers
  total_velocity_outliers <- 0
  for (comp_type in names(gene_velocity_outliers)) {
    for (trans in names(gene_velocity_outliers[[comp_type]])) {
      total_velocity_outliers <- total_velocity_outliers + 
        sum(gene_velocity_outliers[[comp_type]][[trans]], na.rm = TRUE)
    }
  }
  summary_info$n_velocity_outliers <- total_velocity_outliers
  
  # Save summary
  save_final_results(summary_info, project_id, "gene_summary", 
                    GENE_OUTPUT_DIR, norm_method)
  
  # Text summary
  summary_file <- file.path(proj_dir, paste0(project_id, "_gene_summary", norm_suffix, ".txt"))
  sink(summary_file)
  cat(sprintf("GENE TRAJECTORY ANALYSIS - %s\n", project_id))
  cat(sprintf("Normalization: %s\n", norm_display))
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  cat(sprintf("Total genes: %d\n", summary_info$n_genes))
  cat(sprintf("Ortholog genes: %d (%.1f%%)\n", 
              summary_info$n_ortholog_genes,
              summary_info$n_ortholog_genes / summary_info$n_genes * 100))
  cat(sprintf("Stages analyzed: %d/%d\n", summary_info$n_stages_analyzed, length(STAGES)))
  cat(sprintf("Constant healthy reference: %s\n", ifelse(use_constant_healthy, "Yes", "No")))
  cat(sprintf("\n=== TOTAL DISTANCE OUTLIERS (scaled by within-family gene SD) ===\n"))
  
  for (i in seq_along(comp_types)) {
    comp_type <- comp_types[i]
    comp_label <- comp_labels[i]
    
    if (paste0("n_outliers_", comp_type) %in% names(summary_info)) {
      cat(sprintf("\n%s:\n", comp_label))
      cat(sprintf("  Outliers: %d (%.1f%%)\n",
                  summary_info[[paste0("n_outliers_", comp_type)]],
                  summary_info[[paste0("percent_outliers_", comp_type)]]))
    }
  }
  
  cat(sprintf("\n=== VELOCITY OUTLIERS ===\n"))
  cat(sprintf("  Total outliers across all transitions: %d\n", 
              summary_info$n_velocity_outliers))
  sink()
  
  cat(sprintf("\n  ✓ Gene analysis complete for %s [%s]\n", project_id, norm_display))
  cat(sprintf("  ✓ Results saved to %s\n", proj_dir))
  
  return(list(
    project_id = project_id,
    norm_method = norm_method,
    norm_display = norm_display,
    gene_trajectory_results = gene_trajectory_results,
    stage_shifts = stage_shifts_data,
    velocities = velocities_data,
    outliers = outliers_data,
    summary = summary_info
  ))
}

#' Run gene trajectory analysis for all normalization methods
#' @param ortholog_proteins Vector of ortholog protein IDs
#' @param norm_methods Vector of normalization methods to run
#' @param use_constant_healthy Whether to use constant healthy reference (default TRUE)
#' @return List of results for all projects and normalization methods
run_gene_analysis_all <- function(ortholog_proteins, 
                                  norm_methods = c("raw", "full", "std_log"),
                                  use_constant_healthy = TRUE) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("GENE TRAJECTORY ANALYSIS - ALL NORMALIZATION METHODS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  cat("NOTE: This analysis requires family-level results with the SAME\n")
  cat("      normalization method to be generated first!\n\n")
  
  all_results <- list()
  
  for (norm_method in norm_methods) {
    norm_display <- get_norm_display(norm_method)
    cat(sprintf("\n>>> RUNNING WITH NORMALIZATION: %s <<<\n", norm_display))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    for (project_id in PROJECTS) {
      cat(sprintf("\nProcessing: %s\n", project_id))
      
      tryCatch({
        results <- analyze_gene_trajectories_norm(
          project_id, ortholog_proteins, norm_method, use_constant_healthy
        )
        if (!is.null(results)) {
          if (!project_id %in% names(all_results)) {
            all_results[[project_id]] <- list()
          }
          all_results[[project_id]][[norm_method]] <- results
          cat(sprintf("  ✓ Successfully analyzed %s [%s]\n", project_id, norm_display))
        }
      }, error = function(e) {
        cat(sprintf("  ✗ Error analyzing %s [%s]: %s\n", 
                    project_id, norm_display, e$message))
        print(e)
      })
    }
  }
  
  # Create combined summaries for each normalization method
  for (norm_method in norm_methods) {
    norm_results <- list()
    for (project_id in names(all_results)) {
      if (norm_method %in% names(all_results[[project_id]])) {
        norm_results[[project_id]] <- all_results[[project_id]][[norm_method]]
      }
    }
    if (length(norm_results) > 0) {
      create_gene_combined_summary(norm_results, norm_method)
    }
  }
  
  # Create cross-normalization comparison
  create_gene_norm_comparison(all_results)
  
  return(all_results)
}

#' Create combined summary for one normalization method
#' @param all_results List of results for one normalization method
#' @param norm_method Normalization method
create_gene_combined_summary <- function(all_results, norm_method) {
  
  norm_display <- get_norm_display(norm_method)
  norm_suffix <- get_norm_suffix(norm_method)
  output_dir <- get_norm_output_dir(GENE_OUTPUT_DIR, norm_method)
  
  summary_data <- data.frame()
  
  for (project_id in names(all_results)) {
    results <- all_results[[project_id]]
    
    summary_row <- data.frame(
      Project = project_id,
      Normalization = norm_display,
      Genes = results$summary$n_genes,
      Ortholog_Genes = results$summary$n_ortholog_genes,
      Stages_Analyzed = results$summary$n_stages_analyzed,
      Constant_Healthy = results$summary$use_constant_healthy,
      Velocity_Outliers = results$summary$n_velocity_outliers,
      stringsAsFactors = FALSE
    )
    
    comp_types <- c("own", "fam", "orth")
    comp_names <- c("Gene_vs_Own", "Gene_vs_Family", "Gene_vs_Ortholog")
    
    for (i in seq_along(comp_types)) {
      comp_type <- comp_types[i]
      comp_name <- comp_names[i]
      
      n_outliers_col <- paste0("n_outliers_", comp_type)
      percent_col <- paste0("percent_outliers_", comp_type)
      
      if (n_outliers_col %in% names(results$summary)) {
        summary_row[[paste0("Outliers_", comp_name)]] <- results$summary[[n_outliers_col]]
        summary_row[[paste0("Percent_Outliers_", comp_name)]] <- results$summary[[percent_col]]
      } else {
        summary_row[[paste0("Outliers_", comp_name)]] <- NA
        summary_row[[paste0("Percent_Outliers_", comp_name)]] <- NA
      }
    }
    
    summary_data <- rbind(summary_data, summary_row)
  }
  
  # Save
  saveRDS(summary_data, file.path(output_dir, paste0("project_summary", norm_suffix, ".rds")))
  write.csv(summary_data,
            file.path(output_dir, paste0("project_summary", norm_suffix, ".csv")),
            row.names = FALSE)
  
  cat(sprintf("\n✓ Combined gene summary saved for [%s] to %s\n", 
              norm_display, output_dir))
  
  return(summary_data)
}

#' Create cross-normalization comparison for gene analysis
#' @param all_results Nested list: all_results[[project_id]][[norm_method]]
create_gene_norm_comparison <- function(all_results) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("GENE CROSS-NORMALIZATION COMPARISON\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Create comparison directory
  comp_dir <- file.path(dirname(GENE_OUTPUT_DIR), "gene_norm_comparison")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Collect data
  comp_data <- data.frame()
  
  for (project_id in names(all_results)) {
    for (norm_method in names(all_results[[project_id]])) {
      results <- all_results[[project_id]][[norm_method]]
      norm_display <- get_norm_display(norm_method)
      
      comp_data <- rbind(comp_data, data.frame(
        Project = project_id,
        Normalization = norm_display,
        Norm_Method = norm_method,
        Genes = results$summary$n_genes,
        Outliers_Own = results$summary$n_outliers_own,
        Outliers_Own_Pct = results$summary$percent_outliers_own,
        Outliers_Fam = results$summary$n_outliers_fam,
        Outliers_Fam_Pct = results$summary$percent_outliers_fam,
        Outliers_Orth = ifelse(!is.null(results$summary$n_outliers_orth), 
                              results$summary$n_outliers_orth, NA),
        Outliers_Orth_Pct = ifelse(!is.null(results$summary$percent_outliers_orth), 
                                  results$summary$percent_outliers_orth, NA),
        Velocity_Outliers = results$summary$n_velocity_outliers,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save comparison
  saveRDS(comp_data, file.path(comp_dir, "gene_norm_comparison.rds"))
  write.csv(comp_data, file.path(comp_dir, "gene_norm_comparison.csv"), row.names = FALSE)
  
  cat(sprintf("✓ Gene normalization comparison saved to %s\n", comp_dir))
  
  return(comp_data)
}

# ==================== MAIN EXECUTION ====================

if (sys.nframe() == 0) {
  source("tensoromics_functions.R")
  ortholog_proteins <- load_ortholog_info(ORTHOLOG_FILE)

  print(paste0("Number of Orthologs: ", length(ortholog_proteins)))
  
  # Run all three normalization methods
  gene_results <- run_gene_analysis_all(
    ortholog_proteins, 
    norm_methods = c("raw", "full", "std_log", "log"),
    use_constant_healthy = TRUE
  )
  
  cat(sprintf("\n✓ Gene analysis complete. Analyzed %d projects with 3 normalization methods\n", 
              length(gene_results)))
}