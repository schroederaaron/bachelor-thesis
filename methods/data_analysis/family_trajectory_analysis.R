# family_trajectory_analysis.R
# Main family-level trajectory analysis
# WITH NORMALIZATION LOOP: runs all 3 normalization methods

library(dplyr)

source("config.R")
source("utils.R")

#' Analyze family trajectories for one normalization method
#' @param project_id Project ID
#' @param ortholog_proteins Vector of ortholog protein IDs
#' @param norm_method Normalization method: "raw", "full", or "std_log"
#' @param use_constant_healthy Whether to use constant healthy reference
#' @return List with analysis results
analyze_family_trajectories_norm <- function(project_id, ortholog_proteins, 
                                             norm_method = "raw", 
                                             use_constant_healthy = TRUE) {
  
  norm_display <- get_norm_display(norm_method)
  norm_suffix <- get_norm_suffix(norm_method)
  
  cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
  cat(sprintf("FAMILY TRAJECTORY ANALYSIS: %s [%s]\n", project_id, norm_display))
  cat(sprintf("%s\n\n", paste(rep("=", 60), collapse = "")))
  
  # Create project directory (normalization-specific)
  output_dir <- get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method)
  proj_dir <- file.path(output_dir, project_id)
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Initialize storage
  cancer_results_all <- list()
  healthy_results_all <- list()
  healthy_constant <- NULL
  
  cancer_means_all <- list()
  cancer_means_ortholog <- list()
  healthy_means_all <- list()
  healthy_means_ortholog <- list()
  healthy_gene_sds_all <- list()
  healthy_gene_sds_ortholog <- list()
  healthy_n_genes_all <- list()
  healthy_n_genes_ortholog <- list()
  
  shifts_all <- list()
  shifts_ortholog <- list()
  
  stage_outliers_all <- list()
  stage_outliers_ortholog <- list()
  stage_outlier_directions_all <- list()
  stage_outlier_directions_ortholog <- list()
  
  loess_data_all <- NULL
  loess_data_ortholog <- NULL
  
  # ===== LOAD CONSTANT HEALTHY REFERENCE =====
  if (use_constant_healthy) {
    cat("  Loading constant healthy reference...\n")
    healthy_constant_data <- load_stage_data(project_id, STAGES[1], "healthy", 
                                            use_constant_healthy = TRUE, 
                                            norm_method = norm_method)
    if (!is.null(healthy_constant_data)) {
      healthy_constant <- compute_family_means_with_orthologs(
        healthy_constant_data, 
        ortholog_proteins
      )
      cat(sprintf("  ✓ Loaded constant healthy reference (%d samples) [%s]\n", 
                  healthy_constant$n_samples, norm_display))
    } else {
      cat("  Warning: Could not load constant healthy reference. Using stage-specific healthy data.\n")
      use_constant_healthy <- FALSE
    }
  }
  
  # ===== PROCESS EACH STAGE =====
  for (stage in STAGES) {
    cat(sprintf("\n  Processing %s...\n", stage))
    
    # Load cancer data with normalization
    cancer_data <- load_stage_data(project_id, stage, "cancer", 
                                  use_constant_healthy = FALSE, 
                                  norm_method = norm_method)
    if (is.null(cancer_data)) {
      cat(sprintf("    Warning: Missing cancer data for %s\n", stage))
      next
    }
    
    # Load healthy data
    if (use_constant_healthy && !is.null(healthy_constant)) {
      healthy_results <- healthy_constant
      cat("    Using constant healthy reference\n")
    } else {
      healthy_data <- load_stage_data(project_id, stage, "healthy", 
                                     use_constant_healthy = FALSE, 
                                     norm_method = norm_method)
      if (is.null(healthy_data)) {
        cat(sprintf("    Warning: Missing healthy data for %s\n", stage))
        next
      }
      healthy_results <- compute_family_means_with_orthologs(
        healthy_data, 
        ortholog_proteins
      )
    }
    
    # Compute cancer family means
    cancer_results <- compute_family_means_with_orthologs(
      cancer_data, 
      ortholog_proteins
    )
    
    if (is.null(cancer_results) || is.null(healthy_results)) {
      cat(sprintf("    Warning: Could not compute means for %s\n", stage))
      next
    }
    
    # Store results
    cancer_results_all[[stage]] <- cancer_results
    healthy_results_all[[stage]] <- healthy_results
    
    cancer_means_all[[stage]] <- cancer_results$family_means_all
    cancer_means_ortholog[[stage]] <- cancer_results$family_means_ortholog
    healthy_means_all[[stage]] <- healthy_results$family_means_all
    healthy_means_ortholog[[stage]] <- healthy_results$family_means_ortholog
    
    healthy_gene_sds_all[[stage]] <- healthy_results$family_gene_sds_all
    healthy_gene_sds_ortholog[[stage]] <- healthy_results$family_gene_sds_ortholog
    healthy_n_genes_all[[stage]] <- healthy_results$family_n_genes_all
    healthy_n_genes_ortholog[[stage]] <- healthy_results$family_n_genes_ortholog
    
    n_families <- cancer_results$n_families
    
    # ===== ALL GENES =====
    cat("    Computing normalized shifts for all genes...\n")

    # Zuerst: Berechne die Gene Statistics für ALLE Gene
    gene_stats_all <- compute_family_gene_expression_stats(
      expr_matrix = healthy_results$gene_means,  # Vektor von Gen-Mittelwerten
      gene_to_fam = healthy_results$gene_to_fam,
      n_families = n_families,
      ortholog_set = NULL,  # Alle Gene
      min_genes_per_family = MIN_GENES_PER_ALL_FAM,
      sd_floor = NULL
    )

    all_genes_result <- compute_family_normalized_shifts(
      cancer_means = cancer_results$family_means_all,
      healthy_means = healthy_results$family_means_all,
      family_stats = gene_stats_all, 
      n_families = n_families,
      percentile = OUTLIER_PERCENTILE
    )

    shifts_all[[stage]] <- all_genes_result$normalized_shifts
    stage_outliers_all[[stage]] <- all_genes_result$family_outliers
    stage_outlier_directions_all[[stage]] <- all_genes_result$family_outlier_direction

    if (is.null(loess_data_all) && !is.null(all_genes_result$loess_x)) {
      loess_data_all <- list(
        smoothed_sds = all_genes_result$smoothed_sds,
        family_means = all_genes_result$family_means_healthy,
        loess_x = all_genes_result$loess_x,
        loess_y = all_genes_result$loess_y,
        loess_y_smoothed = gene_stats_all$smoothed_sds,
        valid_families = all_genes_result$valid_families_smooth,
        lower_threshold = all_genes_result$lower_threshold,
        upper_threshold = all_genes_result$upper_threshold
      )
    }

    # ===== ORTHOLOGS ONLY =====
    if (sum(cancer_results$ortholog_counts) >= 10) {
      cat("    Computing normalized shifts for orthologs only...\n")
      
      # Berechne Gene Statistics für ORTHOLOGS ONLY
      gene_stats_ortholog <- compute_family_gene_expression_stats(
        expr_matrix = healthy_results$gene_means,
        gene_to_fam = healthy_results$gene_to_fam,
        n_families = n_families,
        ortholog_set = healthy_results$is_ortholog,  # Nur Ortholog-Gene
        min_genes_per_family = MIN_GENES_PER_ORTH_FAM,
        sd_floor = NULL
      )
      
      ortholog_result <- compute_family_normalized_shifts(
        cancer_means = cancer_results$family_means_ortholog,
        healthy_means = healthy_results$family_means_ortholog,
        family_stats = gene_stats_ortholog,  # <- Auch hier family_stats
        n_families = n_families,
        percentile = OUTLIER_PERCENTILE
      )
      
      shifts_ortholog[[stage]] <- ortholog_result$normalized_shifts
      stage_outliers_ortholog[[stage]] <- ortholog_result$family_outliers
      stage_outlier_directions_ortholog[[stage]] <- ortholog_result$family_outlier_direction
      
      if (is.null(loess_data_ortholog) && !is.null(ortholog_result$loess_x)) {
        loess_data_ortholog <- list(
          smoothed_sds = ortholog_result$smoothed_sds,
          raw_gene_sds = gene_stats_ortholog$raw_sds,
          family_means = ortholog_result$family_means_healthy,
          loess_x = ortholog_result$loess_x,
          loess_y = ortholog_result$loess_y,
          loess_y_smoothed = gene_stats_ortholog$smoothed_sds,
          valid_families = ortholog_result$valid_families_smooth,
          lower_threshold = ortholog_result$lower_threshold,
          upper_threshold = ortholog_result$upper_threshold
        )
      }
    } else {
      cat(sprintf("    Warning: Not enough orthologs for reliable analysis\n"))
      shifts_ortholog[[stage]] <- rep(NA, n_families)
      stage_outliers_ortholog[[stage]] <- rep(FALSE, n_families)
      stage_outlier_directions_ortholog[[stage]] <- rep(NA, n_families)
    }
    
    # Save intermediate results with normalization suffix
    stage_results <- list(
      stage = stage,
      norm_method = norm_method,
      norm_display = norm_display,
      cancer = cancer_results,
      healthy = healthy_results,
      shifts_all = shifts_all[[stage]],
      shifts_ortholog = shifts_ortholog[[stage]],
      outliers_all = stage_outliers_all[[stage]],
      outliers_ortholog = stage_outliers_ortholog[[stage]],
      outlier_directions_all = stage_outlier_directions_all[[stage]],
      outlier_directions_ortholog = stage_outlier_directions_ortholog[[stage]],
      loess_data_all = loess_data_all,
      loess_data_ortholog = loess_data_ortholog
    )
    
    save_intermediate_results(stage_results, project_id, 
                             paste0("stage_", gsub(" ", "_", stage)), 
                             FAMILY_OUTPUT_DIR, norm_method)
  }
  
  # Check if we have enough data
  if (length(cancer_means_all) < 2) {
    cat("  Insufficient data for trajectory analysis\n")
    return(NULL)
  }
  
  # ===== COMPUTE TRAJECTORY METRICS =====
  cat("\n  Computing trajectory metrics...\n")
  
  n_families <- cancer_results$n_families
  common_families <- 1:n_families
  first_stage <- names(cancer_means_all)[1]
  ortholog_counts <- cancer_results_all[[first_stage]]$ortholog_counts
  
  # Initialize results
  trajectory_results <- data.frame(
    family_id = common_families,
    ortholog_count = ortholog_counts,
    has_orthologs = ortholog_counts > 0,
    n_genes = healthy_results_all[[first_stage]]$family_n_genes_all,
    total_distance_all = numeric(n_families),
    total_distance_ortholog = numeric(n_families),
    total_distance_scaled_all = numeric(n_families),
    total_distance_scaled_ortholog = numeric(n_families),
    is_total_outlier_all = logical(n_families),
    is_total_outlier_ortholog = logical(n_families),
    n_stages_all = integer(n_families),
    n_stages_ortholog = integer(n_families),
    stringsAsFactors = FALSE
  )
  
  # Stage shifts matrices
  stage_shifts_all <- matrix(NA, nrow = n_families, ncol = length(STAGES))
  stage_shifts_ortholog <- matrix(NA, nrow = n_families, ncol = length(STAGES))
  colnames(stage_shifts_all) <- STAGES
  colnames(stage_shifts_ortholog) <- STAGES
  
  stage_outlier_matrix_all <- matrix(FALSE, nrow = n_families, ncol = length(STAGES))
  stage_outlier_matrix_ortholog <- matrix(FALSE, nrow = n_families, ncol = length(STAGES))
  colnames(stage_outlier_matrix_all) <- STAGES
  colnames(stage_outlier_matrix_ortholog) <- STAGES
  
  # Compute per-family metrics
  for (f in common_families) {
    family_shifts_all <- numeric()
    family_shifts_ortholog <- numeric()
    stage_count_all <- 0
    stage_count_ortholog <- 0
    
    for (stage in STAGES) {
      if (stage %in% names(shifts_all)) {
        shift_all <- shifts_all[[stage]][f]
        shift_ortholog <- shifts_ortholog[[stage]][f]
        
        if (!is.na(shift_all)) {
          stage_shifts_all[f, stage] <- shift_all
          family_shifts_all <- c(family_shifts_all, shift_all)
          stage_count_all <- stage_count_all + 1
          
          if (stage %in% names(stage_outliers_all)) {
            stage_outlier_matrix_all[f, stage] <- stage_outliers_all[[stage]][f]
          }
        }
        
        if (!is.na(shift_ortholog)) {
          stage_shifts_ortholog[f, stage] <- shift_ortholog
          family_shifts_ortholog <- c(family_shifts_ortholog, shift_ortholog)
          stage_count_ortholog <- stage_count_ortholog + 1
          
          if (stage %in% names(stage_outliers_ortholog)) {
            stage_outlier_matrix_ortholog[f, stage] <- stage_outliers_ortholog[[stage]][f]
          }
        }
      }
    }
    
    if (length(family_shifts_all) > 0) {
      trajectory_results$total_distance_all[f] <- sum(abs(family_shifts_all))
      trajectory_results$n_stages_all[f] <- stage_count_all
    }
    if (length(family_shifts_ortholog) > 0) {
      trajectory_results$total_distance_ortholog[f] <- sum(abs(family_shifts_ortholog))
      trajectory_results$n_stages_ortholog[f] <- stage_count_ortholog
    }
  }
  
  # ===== SCALED TOTAL DISTANCE OUTLIERS =====
  cat("  Computing scaled total distance outliers...\n")
  
  if (!is.null(loess_data_all$smoothed_sds)) {
    total_outliers_scaled_all <- compute_scaled_total_outliers(
      total_distances = trajectory_results$total_distance_all,
      percentile = OUTLIER_PERCENTILE
    )
    
    trajectory_results$total_distance_scaled_all <- total_outliers_scaled_all$scaled_distances
    trajectory_results$is_total_outlier_all <- total_outliers_scaled_all$is_outlier
  }
  
  if (!is.null(loess_data_ortholog$smoothed_sds)) {
    total_outliers_scaled_ortholog <- compute_scaled_total_outliers(
      total_distances = trajectory_results$total_distance_ortholog,
      percentile = OUTLIER_PERCENTILE
    )
    
    trajectory_results$total_distance_scaled_ortholog <- total_outliers_scaled_ortholog$scaled_distances
    trajectory_results$is_total_outlier_ortholog <- total_outliers_scaled_ortholog$is_outlier
  }
  
  # ===== VELOCITY OUTLIERS =====
  cat("  Computing velocity outliers...\n")
  
  velocity_outliers_all <- compute_velocity_outliers(
    shifts_matrix = stage_shifts_all,
    stages = STAGES,
    percentile = OUTLIER_PERCENTILE
  )
  
  velocity_outliers_ortholog <- compute_velocity_outliers(
    shifts_matrix = stage_shifts_ortholog,
    stages = STAGES,
    percentile = OUTLIER_PERCENTILE
  )
  
  # ===== SAVE RESULTS =====
  cat("  Saving results...\n")
  
  # Add normalization info to results
  trajectory_results$norm_method <- norm_method
  trajectory_results$norm_display <- norm_display
  
  # Save with normalization suffix
  save_final_results(trajectory_results, project_id, "family_results", 
                    FAMILY_OUTPUT_DIR, norm_method)
  
  # Stage shifts
  stage_shifts_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    all_genes = stage_shifts_all,
    orthologs_only = stage_shifts_ortholog,
    stages = STAGES,
    family_ids = common_families
  )
  
  save_final_results(stage_shifts_data, project_id, "family_shifts", 
                    FAMILY_OUTPUT_DIR, norm_method)
  
  # Velocities
  velocities_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    all_genes = velocity_outliers_all$velocities,
    orthologs_only = velocity_outliers_ortholog$velocities,
    outliers_all = velocity_outliers_all$is_outlier,
    outliers_ortholog = velocity_outliers_ortholog$is_outlier,
    thresholds_all = velocity_outliers_all$thresholds,
    thresholds_ortholog = velocity_outliers_ortholog$thresholds,
    transitions = names(velocity_outliers_all$velocities)
  )
  
  save_final_results(velocities_data, project_id, "family_velocities", 
                    FAMILY_OUTPUT_DIR, norm_method)
  
  # LOESS stats - CRITICAL for gene analysis
  loess_stats_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    gene_smoothed_sds_all = loess_data_all$smoothed_sds,
    gene_family_means_all = loess_data_all$family_means,
    loess_x_all = loess_data_all$loess_x,
    loess_y_all = loess_data_all$loess_y,
    loess_y_smoothed_all = loess_data_all$loess_y_smoothed,
    valid_families_all = loess_data_all$valid_families,
    gene_smoothed_sds_ortholog = if (!is.null(loess_data_ortholog)) 
      loess_data_ortholog$smoothed_sds else rep(NA, n_families),
    gene_raw_sds_ortholog = if (!is.null(loess_data_ortholog))
      loess_data_ortholog$raw_gene_sds else rep(NA, n_families),
    gene_family_means_ortholog = if (!is.null(loess_data_ortholog))
      loess_data_ortholog$family_means else rep(NA, n_families),
    loess_x_ortholog = if (!is.null(loess_data_ortholog))
      loess_data_ortholog$loess_x else NULL,
    loess_y_ortholog = if (!is.null(loess_data_ortholog))
      loess_data_ortholog$loess_y else NULL,
    loess_y_smoothed_ortholog = if (!is.null(loess_data_ortholog))
      loess_data_ortholog$loess_y_smoothed else NULL,
    valid_families_ortholog = if (!is.null(loess_data_ortholog))
      loess_data_ortholog$valid_families else NULL,
    gene_to_fam = healthy_results_all[[first_stage]]$gene_to_fam,
    n_families = n_families,
    outlier_thresholds = list(
      lower = loess_data_all$lower_threshold,
      upper = loess_data_all$upper_threshold
    )
  )
  
  save_final_results(loess_stats_data, project_id, "family_loess_stats", 
                    FAMILY_OUTPUT_DIR, norm_method)
  
  # Outliers
  outliers_data <- list(
    norm_method = norm_method,
    norm_display = norm_display,
    stage_outliers = list(
      all_genes = stage_outlier_matrix_all,
      orthologs_only = stage_outlier_matrix_ortholog,
      directions_all = stage_outlier_directions_all,
      directions_ortholog = stage_outlier_directions_ortholog
    ),
    total_outliers = list(
      all_genes = trajectory_results[, c("family_id", "total_distance_scaled_all", "is_total_outlier_all")],
      orthologs_only = trajectory_results[, c("family_id", "total_distance_scaled_ortholog", "is_total_outlier_ortholog")]
    ),
    velocity_outliers = list(
      all_genes = velocities_data$outliers_all,
      orthologs_only = velocities_data$outliers_ortholog
    )
  )
  
  save_final_results(outliers_data, project_id, "family_outliers", 
                    FAMILY_OUTPUT_DIR, norm_method)
  
  # ===== CREATE SUMMARY =====
  cat("  Creating summary...\n")
  
  families_with_orthologs <- sum(trajectory_results$has_orthologs, na.rm = TRUE)
  
  summary_info <- list(
    project_id = project_id,
    norm_method = norm_method,
    norm_display = norm_display,
    n_families = n_families,
    n_stages_analyzed = length(cancer_means_all),
    n_families_with_orthologs = families_with_orthologs,
    use_constant_healthy = use_constant_healthy,
    n_total_outliers_all = sum(trajectory_results$is_total_outlier_all, na.rm = TRUE),
    n_total_outliers_ortholog = sum(trajectory_results$is_total_outlier_ortholog, na.rm = TRUE),
    percent_total_outliers_all = round(
      sum(trajectory_results$is_total_outlier_all, na.rm = TRUE) / n_families * 100, 1),
    percent_total_outliers_ortholog = ifelse(families_with_orthologs > 0,
      round(sum(trajectory_results$is_total_outlier_ortholog, na.rm = TRUE) / 
            families_with_orthologs * 100, 1), 0),
    n_velocity_outliers_all = sum(sapply(velocity_outliers_all$is_outlier, sum, na.rm = TRUE)),
    n_velocity_outliers_ortholog = sum(sapply(velocity_outliers_ortholog$is_outlier, sum, na.rm = TRUE))
  )
  
  # Add stage-specific outlier counts
  for (stage in names(stage_outliers_all)) {
    stage_key <- gsub(" ", "_", stage)
    summary_info[[paste0("stage_", stage_key, "_outliers_all")]] <- 
      sum(stage_outliers_all[[stage]], na.rm = TRUE)
    summary_info[[paste0("stage_", stage_key, "_outliers_ortholog")]] <- 
      sum(stage_outliers_ortholog[[stage]], na.rm = TRUE)
    
    if (stage %in% names(stage_outlier_directions_all)) {
      dirs <- stage_outlier_directions_all[[stage]]
      summary_info[[paste0("stage_", stage_key, "_up_all")]] <- sum(dirs == "up", na.rm = TRUE)
      summary_info[[paste0("stage_", stage_key, "_down_all")]] <- sum(dirs == "down", na.rm = TRUE)
    }
  }
  
  save_final_results(summary_info, project_id, "family_summary", 
                    FAMILY_OUTPUT_DIR, norm_method)
  
  # Text summary
  summary_file <- file.path(get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method), 
                           project_id, paste0(project_id, "_family_summary", norm_suffix, ".txt"))
  sink(summary_file)
  cat(sprintf("FAMILY TRAJECTORY ANALYSIS - %s\n", project_id))
  cat(sprintf("Normalization: %s\n", norm_display))
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  cat(sprintf("Total families: %d\n", summary_info$n_families))
  cat(sprintf("Families with orthologs: %d (%.1f%%)\n", 
              summary_info$n_families_with_orthologs,
              summary_info$n_families_with_orthologs / summary_info$n_families * 100))
  cat(sprintf("Stages analyzed: %d/%d\n", summary_info$n_stages_analyzed, length(STAGES)))
  cat(sprintf("Constant healthy reference: %s\n", ifelse(use_constant_healthy, "Yes", "No")))
  cat(sprintf("\n=== TOTAL DISTANCE OUTLIERS (scaled by within-family gene SD) ===\n"))
  cat(sprintf("All genes: %d (%.1f%%)\n",
              summary_info$n_total_outliers_all, 
              summary_info$percent_total_outliers_all))
  cat(sprintf("Orthologs only: %d (%.1f%%)\n",
              summary_info$n_total_outliers_ortholog, 
              summary_info$percent_total_outliers_ortholog))
  cat(sprintf("\n=== VELOCITY OUTLIERS ===\n"))
  cat(sprintf("All genes: %d total outliers across transitions\n", 
              summary_info$n_velocity_outliers_all))
  cat(sprintf("Orthologs only: %d total outliers across transitions\n", 
              summary_info$n_velocity_outliers_ortholog))
  sink()
  
  cat(sprintf("\n  ✓ Family analysis complete for %s [%s]\n", project_id, norm_display))
  cat(sprintf("  ✓ Results saved to %s\n", 
              file.path(get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method), project_id)))
  
  return(list(
    project_id = project_id,
    norm_method = norm_method,
    norm_display = norm_display,
    trajectory_results = trajectory_results,
    stage_shifts = stage_shifts_data,
    velocities = velocities_data,
    loess_stats = loess_stats_data,
    outliers = outliers_data,
    summary = summary_info
  ))
}

#' Run family trajectory analysis for all normalization methods
#' @param ortholog_proteins Vector of ortholog protein IDs
#' @param norm_methods Vector of normalization methods to run
#' @param use_constant_healthy Whether to use constant healthy reference
#' @return List of results for all projects and normalization methods
run_family_analysis_all <- function(ortholog_proteins, 
                                    norm_methods = c("raw", "full", "std_log", "log"),
                                    use_constant_healthy = TRUE) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("FAMILY TRAJECTORY ANALYSIS - ALL NORMALIZATION METHODS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  all_results <- list()
  
  for (norm_method in norm_methods) {
    norm_display <- get_norm_display(norm_method)
    cat(sprintf("\n>>> RUNNING WITH NORMALIZATION: %s <<<\n", norm_display))
    cat(paste(rep("-", 50), collapse = ""), "\n")
    
    for (project_id in PROJECTS) {
      cat(sprintf("\nProcessing: %s\n", project_id))
      
      tryCatch({
        results <- analyze_family_trajectories_norm(
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
      create_family_combined_summary(norm_results, norm_method)
    }
  }
  
  # Create cross-normalization comparison
  create_family_norm_comparison(all_results)
  
  return(all_results)
}

#' Create combined summary for one normalization method
#' @param all_results List of results for one normalization method
#' @param norm_method Normalization method
create_family_combined_summary <- function(all_results, norm_method) {
  
  norm_display <- get_norm_display(norm_method)
  norm_suffix <- get_norm_suffix(norm_method)
  output_dir <- get_norm_output_dir(FAMILY_OUTPUT_DIR, norm_method)
  
  summary_data <- data.frame()
  
  for (project_id in names(all_results)) {
    results <- all_results[[project_id]]
    
    summary_row <- data.frame(
      Project = project_id,
      Normalization = norm_display,
      Families = results$summary$n_families,
      Families_With_Orthologs = results$summary$n_families_with_orthologs,
      Stages_Analyzed = results$summary$n_stages_analyzed,
      Constant_Healthy = results$summary$use_constant_healthy,
      Total_Outliers_All = results$summary$n_total_outliers_all,
      Total_Outliers_All_Percent = results$summary$percent_total_outliers_all,
      Total_Outliers_Ortholog = results$summary$n_total_outliers_ortholog,
      Total_Outliers_Ortholog_Percent = results$summary$percent_total_outliers_ortholog,
      Velocity_Outliers_All = results$summary$n_velocity_outliers_all,
      Velocity_Outliers_Ortholog = results$summary$n_velocity_outliers_ortholog,
      stringsAsFactors = FALSE
    )
    
    # Add stage-specific outliers
    for (stage in STAGES) {
      stage_key <- gsub(" ", "_", stage)
      col_all <- paste0("stage_", stage_key, "_outliers_all")
      col_ortholog <- paste0("stage_", stage_key, "_outliers_ortholog")
      
      if (col_all %in% names(results$summary)) {
        summary_row[[col_all]] <- results$summary[[col_all]]
        summary_row[[paste0("stage_", stage_key, "_up_all")]] <- 
          results$summary[[paste0("stage_", stage_key, "_up_all")]]
        summary_row[[paste0("stage_", stage_key, "_down_all")]] <- 
          results$summary[[paste0("stage_", stage_key, "_down_all")]]
      }
    }
    
    summary_data <- rbind(summary_data, summary_row)
  }
  
  # Save
  saveRDS(summary_data, file.path(output_dir, paste0("project_summary", norm_suffix, ".rds")))
  write.csv(summary_data,
            file.path(output_dir, paste0("project_summary", norm_suffix, ".csv")),
            row.names = FALSE)
  
  cat(sprintf("\n✓ Combined family summary saved for [%s] to %s\n", 
              norm_display, output_dir))
  
  return(summary_data)
}

#' Create cross-normalization comparison
#' @param all_results Nested list: all_results[[project_id]][[norm_method]]
create_family_norm_comparison <- function(all_results) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CROSS-NORMALIZATION COMPARISON\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Create comparison directory
  comp_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "family_norm_comparison")
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
        Families = results$summary$n_families,
        Outliers_All = results$summary$n_total_outliers_all,
        Outliers_All_Pct = results$summary$percent_total_outliers_all,
        Outliers_Ortholog = results$summary$n_total_outliers_ortholog,
        Outliers_Ortholog_Pct = results$summary$percent_total_outliers_ortholog,
        Velocity_Outliers = results$summary$n_velocity_outliers_all,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Save comparison
  saveRDS(comp_data, file.path(comp_dir, "family_norm_comparison.rds"))
  write.csv(comp_data, file.path(comp_dir, "family_norm_comparison.csv"), row.names = FALSE)
  
  cat(sprintf("✓ Family normalization comparison saved to %s\n", comp_dir))
  
  return(comp_data)
}

# If run directly
if (sys.nframe() == 0) {
  ortholog_proteins <- load_ortholog_info(ORTHOLOG_FILE)

  cat("Number of orthologs: ", length(ortholog_proteins))
  
  # Run all three normalization methods
  family_results <- run_family_analysis_all(
    ortholog_proteins, 
    norm_methods = c("raw", "full", "std_log", "log"),
    use_constant_healthy = TRUE
  )
  
  cat(sprintf("\n✓ Family analysis complete. Analyzed %d projects with 4 normalization methods\n", 
              length(family_results)))
  warnings()
}