# utils.R
# Utility functions for Family & Gene Trajectory Analysis
source("tensoromics_functions.R")
source("config.R")

# ==================== NORMALIZATION FUNCTIONS ====================

#' Apply normalization to expression matrix
#' @param expr_matrix samples × genes expression matrix
#' @param norm_method Normalization method: "raw", "full", "std_log", or "log"
#' @return Normalized expression vector (named)
apply_normalization <- function(expr_matrix, norm_method = "raw", apply_mean = TRUE) {
  
  if (norm_method == "raw") {
    # No normalization - use raw data
    if(apply_mean){
      cat("Mean applied \n")
      gene_means <- colMeans(expr_matrix)
      return(gene_means)
    } else {
      return(expr_matrix)
    }
    
  } else if (norm_method == "full") {
    # Method 1: gene-wise scaling + quantile + log2
    cat("      Applying FULL normalization: gene-wise scaling + quantile + log2\n")

    gene_ids <- colnames(expr_matrix)
    expr_t <- t(expr_matrix)
    expr_scaled <- tox_normalize_by_std_dev(expr_t)
    expr_quantile <- tox_quantile_normalization(expr_scaled)
    gene_means <- rowMeans(expr_quantile)
    expr_log <- tox_log2_transformation(matrix(gene_means, ncol=1))
    rownames(expr_log) <- gene_ids
    return(t(expr_log))
    
  } else if (norm_method == "std_log") {
    # Method 2: gene-wise scaling + log2 (no quantile)
    cat("      Applying STD_LOG normalization: gene-wise scaling + log2\n")
    
    expr_t <- t(expr_matrix)
    org_names <- rownames(expr_t)
    expr_scaled <- tox_normalize_by_std_dev(expr_t)
    expr_means <- rowMeans(expr_scaled)
    expr_log <- tox_log2_transformation(matrix(expr_means, ncol = 1))

    rownames(expr_log) <- org_names
    return(t(expr_log))

  } else if(norm_method == "log"){
    # Method 3: log only
    expr_t <- t(expr_matrix)
    original_names <- rownames(expr_t)
    expr_mean <- rowMeans(expr_t)
    expr_log <- tox_log2_transformation(matrix(expr_mean, ncol=1))
    rownames(expr_log) <- original_names
    return(t(expr_log))
  } else {
    stop(paste("Unknown normalization method:", norm_method))
  }
}

#' Get file suffix for normalization method
#' @param norm_method Normalization method
#' @return String suffix for file names
get_norm_suffix <- function(norm_method) {
  switch(norm_method,
         "raw" = "_raw",
         "full" = "_fullnorm",
         "std_log" = "_stdlognorm",
         "log" = "_lognorm",
         stop(paste("Unknown normalization method:", norm_method)))
}

#' Get display name for normalization method
#' @param norm_method Normalization method
#' @return String display name
get_norm_display <- function(norm_method) {
  switch(norm_method,
         "raw" = "Raw Data",
         "full" = "Gene-wise Scaling + Quantile + Log2",
         "std_log" = "Gene-wise Scaling + Log2",
         "log" = "Log2",
         stop(paste("Unknown normalization method:", norm_method)))
}

# ==================== DATA LOADING ====================

#' Load stage-specific data with normalization
#' @param project_id e.g., "TCGA-BRCA"
#' @param stage e.g., "Stage I"
#' @param data_type "cancer" or "healthy"
#' @param use_constant_healthy If TRUE, load reference healthy file
#' @param norm_method Normalization method: "raw", "full", or "std_log"
#' @return List with: expression_vectors, gene_ids, gene_to_fam, n_families, norm_method
load_stage_data <- function(project_id, stage, data_type = "cancer", 
                            use_constant_healthy = FALSE, norm_method = "raw", apply_mean = TRUE, normalize = TRUE) {
  
  # Construct filename based on pattern
  if (data_type == "healthy") {
    if (use_constant_healthy) {
      filename <- paste0("healthy_", project_id, "_reference.rds")
      if (!file.exists(file.path(BASE_DATA_DIR, project_id, filename))) {
        filename <- paste0("healthy_", project_id, ".rds")
      }
    } else {
      filename <- paste0("healthy_", project_id, ".rds")
    }
  } else {
    filename <- paste0(project_id, "-", gsub(" ", "-", stage), ".rds")
  }
  
  file_path <- file.path(BASE_DATA_DIR, project_id, filename)
  
  if (!file.exists(file_path)) {
    cat(sprintf("    Warning: File not found: %s\n", file_path))
    return(NULL)
  }
  
  # Load raw data
  data_obj <- readRDS(file_path)
  
  # Validate structure
  required_fields <- c("expression_vectors", "gene_ids", "gene_to_fam", "n_families")
  missing_fields <- setdiff(required_fields, names(data_obj))
  
  if (length(missing_fields) > 0) {
    cat(sprintf("    Warning: Missing fields in %s: %s\n", 
                file_path, paste(missing_fields, collapse = ", ")))
    return(NULL)
  }
  
  if (!is.matrix(data_obj$expression_vectors)) {
    cat(sprintf("    Warning: expression_vectors is not a matrix in %s\n", file_path))
    return(NULL)
  }
  
  # Apply normalization
  if(normalize){
    cat(sprintf("    Applying normalization: %s\n", get_norm_display(norm_method)))
    expr_norm <- apply_normalization(data_obj$expression_vectors, norm_method, apply_mean)
    data_obj$expression_vectors <- expr_norm
  }
  
  # Add normalization info to object - expr_norm is now a named vector
  data_obj$norm_method <- norm_method
  data_obj$norm_display <- get_norm_display(norm_method)

  if(apply_mean){
    cat(sprintf("    Loaded %s %s data: %d genes [%s]\n", 
              ifelse(use_constant_healthy, "reference", stage), 
              data_type, 
              length(data_obj$expression_vectors),
              get_norm_display(norm_method)))
  } else{
    cat(sprintf("    Loaded %s %s data: %d samples %d genes [%s]\n", 
              ifelse(use_constant_healthy, "reference", stage), 
              data_type,
              nrow(data_obj$expression_vectors),
              ncol(data_obj$expression_vectors),
              get_norm_display(norm_method)))
  }
  
  
  return(data_obj)
}

#' Load ortholog information
#' @param ortholog_file Path to the ortholog file
#' @return Character vector of ortholog protein IDs
load_ortholog_info <- function(ortholog_file) {
  if (!file.exists(ortholog_file)) {
    stop(sprintf("Ortholog file not found: %s", ortholog_file))
  }
  
  ortholog_data <- read.table(ortholog_file, sep = "\t", 
                              header = FALSE, stringsAsFactors = FALSE)
  ortholog_proteins <- unique(ortholog_data[, 2])
  
  cat(sprintf("✓ Loaded %d unique ortholog proteins\n", length(ortholog_proteins)))
  return(ortholog_proteins)
}

# ==================== CORE STATISTICS ====================

#' Compute within-family gene expression variability and LOESS smoothing
#' @return List with smoothed_sds, family_gene_means, raw_sds, valid_families, etc.
compute_family_gene_expression_stats <- function(expr_matrix, gene_to_fam, n_families,
                                                ortholog_set = NULL,
                                                min_genes_per_family = 3,
                                                sd_floor = NULL) {
  
  gene_means <- expr_matrix  # Vektor der Gen-Mittelwerte
  n_genes <- length(gene_means)
  
  # Bestimme welche Gene verwendet werden sollen
  if (!is.null(ortholog_set)) {
    use_gene <- ortholog_set
  } else {
    use_gene <- rep(TRUE, n_genes)
  }
  
  valid_genes <- use_gene & !is.na(gene_means) & is.finite(gene_means)
  expression_valid <- gene_means[valid_genes]
  gene_to_fam_valid <- gene_to_fam[valid_genes]
  
  if (length(expression_valid) < 10) {
    warning(sprintf("Only %d valid genes for family scaling", length(expression_valid)))
  }
  
  scaling_result <- tox_compute_family_scaling(
    distances = expression_valid, 
    gene_to_fam = gene_to_fam_valid,
    n_families = n_families
  )
  # Prüfe auf Extremwerte
  if (max(scaling_result$dscale, na.rm=TRUE) < 0.01) {
    warning("dscale values are extremely small (< 0.01)!")
  }
  
  if (min(scaling_result$dscale, na.rm=TRUE) == 0) {
    warning("dscale contains zeros!")
  }
  
  # Wende sd_floor an falls gewünscht
  if (!is.null(sd_floor)) {
    cat(sprintf("  Applying sd_floor = %.6f\n", sd_floor))
    before <- sum(scaling_result$dscale < sd_floor, na.rm=TRUE)
    scaling_result$dscale[scaling_result$dscale < sd_floor] <- sd_floor
    cat(sprintf("    %d values floored\n", before))
  }
  
  # Ergebnis
  list(
    smoothed_sds = scaling_result$dscale,
    loess_x = scaling_result$loess_x,
    loess_y = scaling_result$loess_y,
    indices_used = scaling_result$indices_used,
    family_gene_sds = scaling_result$dscale,
    family_n_genes = tabulate(gene_to_fam, n_families),
    gene_means = gene_means,
    gene_to_fam = gene_to_fam,
    valid_families = which(scaling_result$indices_used >= min_genes_per_family)
  )
}

# ==================== FAMILY MEANS WITH NORMALIZATION ====================

#' Compute family means using tox_group_centroid with normalization
#' @param stage_obj Loaded stage data object (with normalized expression)
#' @param ortholog_proteins Vector of ortholog protein IDs
#' @return List with family means, gene SDs, and all needed statistics
compute_family_means_with_orthologs <- function(stage_obj, ortholog_proteins) {
  if (is.null(stage_obj)) return(NULL)
  
  expression_vectors <- stage_obj$expression_vectors  # named vector of gene means
  gene_to_fam <- stage_obj$gene_to_fam
  gene_ids <- stage_obj$gene_ids
  n_families <- stage_obj$n_families
  n_genes <- length(gene_ids)
  norm_method <- stage_obj$norm_method
  norm_display <- stage_obj$norm_display
  
  # Identify ortholog genes
  is_ortholog <- gene_ids %in% ortholog_proteins
  
  cat(sprintf("    Found %d ortholog genes out of %d total genes\n", 
              sum(is_ortholog), n_genes))
  
  # ===== PART 1: COMPUTE FAMILY MEANS =====
  # With gene means vector, we can compute family means directly
  family_means_all <- numeric(n_families)
  family_means_ortholog <- numeric(n_families)

  family_means_all <- tox_group_centroid(
    expression_vectors = matrix(expression_vectors, nrow = 1),
    gene_to_family = gene_to_fam,
    n_families = n_families,
    ortholog_set = rep(TRUE, n_genes),
    mode = 'all'
  )$centroid_matrix

  family_means_ortholog <- tox_group_centroid(
    expression_vectors = matrix(expression_vectors, nrow = 1),
    gene_to_family = gene_to_fam,
    n_families = n_families,
    ortholog_set = is_ortholog,
    mode = 'orthologs'
  )$centroid_matrix
  
  # ===== PART 2: COMPUTE WITHIN-FAMILY GENE EXPRESSION VARIABILITY =====
  
  # All genes
  gene_stats_all <- compute_family_gene_expression_stats(
    expr_matrix = expression_vectors,
    gene_to_fam = gene_to_fam,
    n_families = n_families,
    ortholog_set = NULL
  )
  
  # Orthologs only
  gene_stats_ortholog <- compute_family_gene_expression_stats(
    expr_matrix = expression_vectors,
    gene_to_fam = gene_to_fam,
    n_families = n_families,
    ortholog_set = is_ortholog
  )
  
  # Count orthologs per family
  ortholog_counts <- integer(n_families)
  for (f in 1:n_families) {
    ortholog_counts[f] <- sum(gene_to_fam == f & is_ortholog)
  }
  
  # Return everything with CONSISTENT naming
  list(
    # Normalization info
    norm_method = norm_method,
    norm_display = norm_display,
    
    # Family means
    family_means_all = family_means_all,
    family_means_ortholog = family_means_ortholog,
    
    # Within-family gene expression SDs - ALL use consistent naming
    family_gene_sds_all = gene_stats_all$family_gene_sds,           # = smoothed_sds
    family_gene_sds_ortholog = gene_stats_ortholog$family_gene_sds, # = smoothed_sds
    smoothed_sds_all = gene_stats_all$smoothed_sds,                 # explizit
    smoothed_sds_ortholog = gene_stats_ortholog$smoothed_sds,       # explizit
    
    family_gene_means_all = gene_stats_all$family_gene_means,
    family_gene_means_ortholog = gene_stats_ortholog$family_gene_means,
    family_n_genes_all = gene_stats_all$family_n_genes,
    family_n_genes_ortholog = gene_stats_ortholog$family_n_genes,
    
    # LOESS diagnostic data
    loess_x_all = gene_stats_all$loess_x,
    loess_y_all = gene_stats_all$loess_y,
    loess_x_ortholog = gene_stats_ortholog$loess_x,
    loess_y_ortholog = gene_stats_ortholog$loess_y,
    
    # Valid families
    valid_families_all = gene_stats_all$valid_families,
    valid_families_ortholog = gene_stats_ortholog$valid_families,
    
    # Raw gene data
    gene_means = expression_vectors,
    gene_to_fam = gene_to_fam,
    gene_ids = gene_ids,
    is_ortholog = is_ortholog,
    ortholog_counts = ortholog_counts,
    
    # Metadata
    n_families = n_families,
    n_genes = n_genes,
    n_orthologs = sum(is_ortholog)
  )
}

# ==================== FAMILY-LEVEL ANALYSIS ====================

#' Compute normalized shifts for families
#' JETZT MIT tox_identify_outliers für konsistente Ausreißererkennung
compute_family_normalized_shifts <- function(cancer_means,
                                             healthy_means,
                                             family_stats,
                                             n_families,
                                             percentile = 95.0) {
  
  # family_stats kommt von compute_family_gene_expression_stats
  # Verwende smoothed_sds (oder family_gene_sds - beides gleich)
  smoothed_sds <- family_stats$smoothed_sds
  
  # Compute family normalized shifts
  normalized_shifts <- numeric(n_families)
  for (f in 1:n_families) {
    if (!is.na(cancer_means[f]) && !is.na(healthy_means[f]) &&
        !is.na(smoothed_sds[f]) && smoothed_sds[f] > 0) {
      normalized_shifts[f] <- (cancer_means[f] - healthy_means[f]) / smoothed_sds[f]
    } else {
      normalized_shifts[f] <- NA
    }
  }
  
  cat(sprintf("    Out of %d families, %d have NA normalized shifts\n", 
              n_families, sum(is.na(normalized_shifts))))
  
  # Detect family outliers mit tox_identify_outliers (zwei-tailed)
  valid_family_idx <- !is.na(normalized_shifts)
  family_outliers <- rep(FALSE, n_families)
  family_outlier_direction <- rep(NA, n_families)
  lower_threshold <- NA
  upper_threshold <- NA
  
  if (sum(valid_family_idx) > 0) {
    valid_shifts <- normalized_shifts[valid_family_idx]
    
    # Zwei-tailed: getrennt für positive und negative
    pos_idx <- valid_shifts > 0
    neg_idx <- valid_shifts < 0
    
    # Up-regulated (positive shifts)
    if (sum(pos_idx) > 0) {
      pos_outliers <- tox_identify_outliers(valid_shifts[pos_idx], percentile)
      
      # Finde die Original-Indices
      orig_pos_indices <- which(valid_family_idx)[pos_idx]
      family_outliers[orig_pos_indices[pos_outliers$is_outlier]] <- TRUE
      family_outlier_direction[orig_pos_indices[pos_outliers$is_outlier]] <- "up"
      upper_threshold <- pos_outliers$threshold
    }
    
    # Down-regulated (negative shifts) - arbeite mit Absolutwerten
    if (sum(neg_idx) > 0) {
      neg_abs <- abs(valid_shifts[neg_idx])
      neg_outliers <- tox_identify_outliers(neg_abs, percentile)
      
      orig_neg_indices <- which(valid_family_idx)[neg_idx]
      family_outliers[orig_neg_indices[neg_outliers$is_outlier]] <- TRUE
      family_outlier_direction[orig_neg_indices[neg_outliers$is_outlier]] <- "down"
      lower_threshold <- -neg_outliers$threshold  # Negativ für untere Grenze
    }
  }
  
  list(
    normalized_shifts = normalized_shifts,
    family_outliers = family_outliers,
    family_outlier_direction = family_outlier_direction,
    smoothed_sds = smoothed_sds,
    family_means_healthy = healthy_means,
    family_means_cancer = cancer_means,
    loess_x = family_stats$loess_x,
    loess_y = family_stats$loess_y,
    valid_families_smooth = family_stats$valid_families,
    lower_threshold = lower_threshold,
    upper_threshold = upper_threshold,
    percentile = percentile,
    n_families = n_families
  )
}

# ==================== VELOCITY ANALYSIS ====================

#' Compute velocities and detect outliers with p-values
compute_velocity_outliers <- function(shifts_matrix, stages, percentile = 95.0) {
  n_items <- nrow(shifts_matrix)
  n_stages <- length(stages)
  
  velocities <- list()
  outliers <- list()
  thresholds <- list()
  p_values_list <- list()
  directions <- list()
  
  for (i in 1:(n_stages - 1)) {
    stage1 <- stages[i]
    stage2 <- stages[i + 1]
    transition <- paste(stage1, "->", stage2)
    
    vel <- shifts_matrix[, stage2] - shifts_matrix[, stage1]
    valid_idx <- !is.na(vel)
    
    if (sum(valid_idx) >= 10) {
      abs_vel <- abs(vel[valid_idx])
      
      # tox_identify_outliers für p-Werte!
      outlier_result <- tox_identify_outliers(
        rdi = abs_vel,
        percentile = percentile
      )
      
      is_outlier <- rep(FALSE, n_items)
      is_outlier[valid_idx] <- outlier_result$is_outlier
      
      p_values_full <- rep(NA, n_items)
      p_values_full[valid_idx] <- outlier_result$p_values
      
      direction <- rep(NA, n_items)
      direction[valid_idx] <- ifelse(vel[valid_idx] > 0, "up", "down")
      
      velocities[[transition]] <- vel
      outliers[[transition]] <- is_outlier
      thresholds[[transition]] <- outlier_result$threshold
      p_values_list[[transition]] <- p_values_full
      directions[[transition]] <- direction
    } else {
      # Zu wenige Werte
      velocities[[transition]] <- vel
      outliers[[transition]] <- rep(FALSE, n_items)
      thresholds[[transition]] <- NA
      p_values_list[[transition]] <- rep(NA, n_items)
      directions[[transition]] <- rep(NA, n_items)
    }
  }
  
  list(
    velocities = velocities,
    is_outlier = outliers,
    thresholds = thresholds,
    p_values = p_values_list,
    direction = directions,
    stages = stages
  )
}

# ==================== TOTAL DISTANCE ANALYSIS ====================

#' Compute scaled total distance outliers with p-values
#' @param total_distances Vector of total distances (already normalized once)
#' @param percentile Percentile threshold
#' @return List with scaled_distances, is_outlier, threshold, p_values
compute_scaled_total_outliers <- function(total_distances, percentile = 95.0) {
  
  n_items <- length(total_distances)
  valid_idx <- !is.na(total_distances)
  
  if (sum(valid_idx) < 10) {
    return(list(
      scaled_distances = total_distances,
      is_outlier = rep(FALSE, n_items),
      threshold = NA,
      p_values = rep(NA, n_items),
      percentile = percentile,
      n_valid = sum(valid_idx)
    ))
  }
  
  # tox_identify_outliers für p-Werte verwenden!
  outlier_result <- tox_identify_outliers(
    rdi = total_distances[valid_idx],
    percentile = percentile
  )
  
  is_outlier <- rep(FALSE, n_items)
  is_outlier[valid_idx] <- outlier_result$is_outlier
  
  p_values_full <- rep(NA, n_items)
  p_values_full[valid_idx] <- outlier_result$p_values
  
  list(
    scaled_distances = total_distances,
    is_outlier = is_outlier,
    threshold = outlier_result$threshold,
    p_values = p_values_full,
    percentile = percentile,
    n_valid = sum(valid_idx)
  )
}


# ==================== FILE HANDLING ====================

#' Get output directory for normalization method
#' @param base_dir Base output directory
#' @param norm_method Normalization method
#' @return Path to normalization-specific output directory
get_norm_output_dir <- function(base_dir, norm_method) {
  norm_suffix <- get_norm_suffix(norm_method)
  file.path(dirname(base_dir), paste0(basename(base_dir), norm_suffix))
}

#' Save intermediate results with normalization info
#' @param results Results object to save
#' @param project_id Project ID
#' @param type Type of results
#' @param output_dir Base output directory
#' @param norm_method Normalization method
save_intermediate_results <- function(results, project_id, type, output_dir, norm_method) {
  # Create normalization-specific directory
  norm_dir <- get_norm_output_dir(output_dir, norm_method)
  proj_dir <- file.path(norm_dir, project_id, "intermediate/")
  
  cat(sprintf("    Saving intermediate results to: %s\n", proj_dir))
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  
  norm_suffix <- get_norm_suffix(norm_method)
  filename <- paste0(project_id, "_", type, norm_suffix, ".rds")
  saveRDS(results, file.path(proj_dir, filename))
}

#' Save final results with normalization info
#' @param results Results object to save
#' @param project_id Project ID
#' @param filename_base Base filename (without extension)
#' @param output_dir Base output directory
#' @param norm_method Normalization method
save_final_results <- function(results, project_id, filename_base, output_dir, norm_method) {
  # Create normalization-specific directory
  norm_dir <- get_norm_output_dir(output_dir, norm_method)
  proj_dir <- file.path(norm_dir, project_id)
  
  cat(sprintf("    Saving final results to: %s\n", proj_dir))
  dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)
  
  norm_suffix <- get_norm_suffix(norm_method)
  
  # Save RDS
  rds_file <- file.path(proj_dir, paste0(project_id, "_", filename_base, norm_suffix, ".rds"))
  saveRDS(results, rds_file)
  
  # For data frames, also save CSV
  if (is.data.frame(results)) {
    csv_file <- file.path(proj_dir, paste0(project_id, "_", filename_base, norm_suffix, ".csv"))
    write.csv(results, csv_file, row.names = FALSE)
  }
}

#' Load gene results for a specific normalization method
#' @param project_id Project ID
#' @param norm_method Normalization method: "raw", "full", "log", or "std_log"
#' @return List with all gene results
load_gene_results_norm <- function(project_id, norm_method) {
  
  norm_display <- get_norm_display(norm_method)
  norm_suffix <- get_norm_suffix(norm_method)
  
  # Get normalization-specific directory
  norm_dir <- get_norm_output_dir(GENE_OUTPUT_DIR, norm_method)
  proj_dir <- file.path(norm_dir, project_id)
  
  if (!dir.exists(proj_dir)) {
    return(NULL)
  }
  
  # Check for required files
  required_files <- c(
    paste0(project_id, "_gene_results", norm_suffix, ".rds"),
    paste0(project_id, "_gene_shifts", norm_suffix, ".rds"),
    paste0(project_id, "_gene_summary", norm_suffix, ".rds")
  )
  
  missing_files <- required_files[!file.exists(file.path(proj_dir, required_files))]
  if (length(missing_files) > 0) {
    return(NULL)
  }
  
  results <- list(
    project_id = project_id,
    norm_method = norm_method,
    norm_display = norm_display,
    norm_suffix = norm_suffix,
    trajectory = readRDS(file.path(proj_dir, paste0(project_id, "_gene_results", norm_suffix, ".rds"))),
    stage_shifts = readRDS(file.path(proj_dir, paste0(project_id, "_gene_shifts", norm_suffix, ".rds"))),
    summary = readRDS(file.path(proj_dir, paste0(project_id, "_gene_summary", norm_suffix, ".rds")))
  )
  
  # Load optional files
  outliers_file <- file.path(proj_dir, paste0(project_id, "_gene_outliers", norm_suffix, ".rds"))
  velocities_file <- file.path(proj_dir, paste0(project_id, "_gene_velocities", norm_suffix, ".rds"))
  
  if (file.exists(outliers_file)) {
    results$outliers <- readRDS(outliers_file)
  }
  
  if (file.exists(velocities_file)) {
    results$velocities <- readRDS(velocities_file)
  }
  
  return(results)
}

# ==================== COLOR MANAGER ====================

create_gene_color_manager <- function(initial_map = NULL) {
  color_env <- new.env()
  color_env$color_map <- list()
  color_env$used_hues <- numeric(0)

  # Helper to ensure gene_ids is a character vector
  safe_as_character <- function(x) {
    if (is.null(x)) return(character(0))
    if (is.matrix(x) || is.data.frame(x)) {
      if (ncol(x) > 1) {
        warning("Input has multiple columns, using first column only")
        x <- x[, 1]
      } else {
        x <- x[, 1]
      }
    }
    x <- as.vector(x)
    result <- as.character(x)
    result <- result[!is.na(result) & result != ""]
    return(result)
  }

  infer_hue <- function(color) {
    rgb <- grDevices::col2rgb(color) / 255
    hsv <- grDevices::rgb2hsv(rgb[1], rgb[2], rgb[3])
    as.numeric(hsv[1]) * 360
  }

  # Calculate perceptual distance between two colors in CIE Lab
  perceptual_distance <- function(color1, color2) {
    rgb1 <- t(grDevices::col2rgb(color1) / 255)
    rgb2 <- t(grDevices::col2rgb(color2) / 255)

    lab1 <- grDevices::convertColor(rgb1, from = "sRGB", to = "Lab", scale.in = 1)
    lab2 <- grDevices::convertColor(rgb2, from = "sRGB", to = "Lab", scale.in = 1)

    sqrt(sum((lab1 - lab2)^2))
  }

  is_color_distinct <- function(candidate_color, existing_colors, min_distance = 20) {
    if (length(existing_colors) == 0) return(TRUE)
    all(vapply(existing_colors, function(existing) {
      perceptual_distance(candidate_color, existing) >= min_distance
    }, logical(1)))
  }

  generate_candidate_palette <- function() {
    hues <- seq(15, 375, length.out = 49)[1:48]
    schemes <- list(
      list(c = 85, l = 55),
      list(c = 70, l = 68),
      list(c = 95, l = 42),
      list(c = 55, l = 76)
    )

    candidates <- character(0)
    candidate_hues <- numeric(0)
    for (scheme in schemes) {
      cols <- grDevices::hcl(h = hues, c = scheme$c, l = scheme$l, fixup = TRUE)
      candidates <- c(candidates, cols)
      candidate_hues <- c(candidate_hues, hues)
    }

    keep <- !duplicated(candidates)
    list(colors = candidates[keep], hues = candidate_hues[keep])
  }

  choose_best_candidate <- function(existing_colors, used_colors = character(0)) {
    candidates <- generate_candidate_palette()
    blocked <- unique(c(existing_colors, used_colors))

    available_idx <- which(!(candidates$colors %in% blocked))
    if (length(available_idx) == 0) {
      stop("No candidate colors available")
    }

    if (length(existing_colors) == 0) {
      idx <- available_idx[1]
      return(list(color = candidates$colors[idx], hue = candidates$hues[idx]))
    }

    min_distances <- vapply(available_idx, function(i) {
      min(vapply(existing_colors, function(existing) {
        perceptual_distance(candidates$colors[i], existing)
      }, numeric(1)))
    }, numeric(1))

    best_local <- which.max(min_distances)
    idx <- available_idx[best_local]
    list(color = candidates$colors[idx], hue = candidates$hues[idx])
  }

  generate_distinct_colors <- function(n, existing_colors) {
    if (n <= 0) return(list(colors = character(0), hues = numeric(0)))

    new_colors <- character(n)
    new_hues <- numeric(n)
    current_colors <- existing_colors

    for (i in seq_len(n)) {
      result <- choose_best_candidate(current_colors, used_colors = new_colors[seq_len(max(0, i - 1))])
      new_colors[i] <- result$color
      new_hues[i] <- result$hue
      current_colors <- c(current_colors, result$color)
    }

    list(colors = new_colors, hues = new_hues)
  }

  load_map <- function(color_map) {
    if (is.null(color_map) || length(color_map) == 0) {
      color_env$color_map <- list()
      color_env$used_hues <- numeric(0)
      return(invisible(FALSE))
    }

    if (is.list(color_map)) {
      color_map <- unlist(color_map)
    }
    color_map <- color_map[!is.na(names(color_map)) & names(color_map) != "" & !is.na(color_map)]

    color_env$color_map <- as.list(color_map)
    color_env$used_hues <- vapply(unname(color_map), infer_hue, numeric(1))
    invisible(TRUE)
  }

  get_colors <- function(gene_ids) {
    gene_ids <- safe_as_character(gene_ids)

    if (length(gene_ids) == 0) {
      return(character(0))
    }

    existing_genes <- names(color_env$color_map)
    new_genes <- setdiff(unique(gene_ids), existing_genes)

    if (length(new_genes) > 0) {
      existing_colors <- unlist(color_env$color_map)
      result <- generate_distinct_colors(length(new_genes), existing_colors)
      new_colors <- result$colors
      names(new_colors) <- new_genes

      color_env$color_map <- c(color_env$color_map, as.list(new_colors))
      color_env$used_hues <- c(color_env$used_hues, result$hues)
    }

    result <- vapply(gene_ids, function(gid) {
      if (gid %in% names(color_env$color_map)) {
        color_env$color_map[[gid]]
      } else {
        "#CCCCCC"
      }
    }, FUN.VALUE = character(1))

    names(result) <- gene_ids
    result
  }

  get_color <- function(gene_id) {
    colors <- get_colors(gene_id)
    if (length(colors) > 0) colors[1] else "#CCCCCC"
  }

  get_full_map <- function() {
    unlist(color_env$color_map)
  }

  save_map <- function(file) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(get_full_map(), file)
    invisible(file)
  }

  reset_colors <- function() {
    color_env$color_map <- list()
    color_env$used_hues <- numeric(0)
    cat("Color map reset\n")
  }

  diagnose_colors <- function(threshold = 20) {
    if (length(color_env$color_map) == 0) {
      cat("No colors assigned yet\n")
      return(NULL)
    }

    colors <- unlist(color_env$color_map)
    n_colors <- length(colors)
    cat(sprintf("Color manager: %d colors assigned\n", n_colors))

    if (any(duplicated(colors))) {
      cat("WARNING: Exact duplicate colors found!\n")
      dup_colors <- unique(colors[duplicated(colors)])
      for (dup in dup_colors) {
        genes_with_dup <- names(colors[colors == dup])
        cat("  Color", dup, "assigned to:", paste(genes_with_dup, collapse = ", "), "\n")
      }
    } else {
      cat("All colors are unique\n")
    }

    if (n_colors > 1) {
      similar_pairs <- list()
      for (i in 1:(n_colors - 1)) {
        for (j in (i + 1):n_colors) {
          dist <- perceptual_distance(colors[i], colors[j])
          if (dist < threshold) {
            similar_pairs[[length(similar_pairs) + 1]] <- list(
              gene1 = names(colors)[i],
              gene2 = names(colors)[j],
              color1 = colors[i],
              color2 = colors[j],
              distance = dist
            )
          }
        }
      }

      if (length(similar_pairs) > 0) {
        cat(sprintf("\nWARNING: Found %d perceptually similar color pairs (Lab distance < %d):\n",
                    length(similar_pairs), threshold))
        for (pair in similar_pairs) {
          cat(sprintf("  %s and %s (distance: %.1f)\n    %s vs %s\n",
                      pair$gene1, pair$gene2, pair$distance,
                      pair$color1, pair$color2))
        }
      } else {
        cat(sprintf("\nAll color pairs are perceptually distinct (Lab distance >= %d)\n", threshold))
      }
    }

    invisible(TRUE)
  }

  fix_similar_colors <- function(threshold = 20) {
    if (length(color_env$color_map) < 2) {
      cat("Not enough colors to check\n")
      return(invisible(FALSE))
    }

    colors <- unlist(color_env$color_map)
    fixed_count <- 0

    for (i in 1:(length(colors) - 1)) {
      for (j in (i + 1):length(colors)) {
        dist <- perceptual_distance(colors[i], colors[j])
        if (dist < threshold) {
          gene_to_fix <- names(colors)[j]
          existing_colors <- colors[names(colors) != gene_to_fix]
          result <- generate_distinct_colors(1, existing_colors)
          color_env$color_map[[gene_to_fix]] <- result$colors[1]
          colors <- unlist(color_env$color_map)
          fixed_count <- fixed_count + 1
        }
      }
    }

    color_env$used_hues <- vapply(unname(unlist(color_env$color_map)), infer_hue, numeric(1))
    cat(sprintf("Fixed %d similar color assignments\n", fixed_count))
    invisible(TRUE)
  }

  load_map(initial_map)

  return(list(
    get_colors = get_colors,
    get_color = get_color,
    get_full_map = get_full_map,
    save_map = save_map,
    load_map = load_map,
    reset_colors = reset_colors,
    diagnose_colors = diagnose_colors,
    fix_similar_colors = fix_similar_colors
  ))
}
# ==================== INITIALIZATION ====================

cat("✓ Utility functions loaded with NORMALIZATION FRAMEWORK\n")
cat("  Available normalization methods:\n")
cat("    - raw: No normalization\n")
cat("    - full: Gene-wise scaling + quantile + log2\n")
cat("    - std_log: Gene-wise scaling + log2\n")
cat("    - log: log2 scaling\n")