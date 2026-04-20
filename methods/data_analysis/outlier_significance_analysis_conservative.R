##############
## OUTDATED ##
##############

#!/usr/bin/env Rscript
# outlier_significance_analysis.R
# Outlier Significance Analysis mit kNN-basierten Noise-p-Werten
# SPEICHERT NULL-DISTANZEN für spätere diagnostische Plots

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

# Für Plot 4: Speichere Perzentile für globale Visualisierung
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
      
      healthy_data_raw <- load_stage_data(
        project_id, 
        STAGES[1], 
        "healthy", 
        use_constant_healthy = TRUE, 
        norm_method = "raw", # load raw for noise computation
        apply_mean = FALSE
      )
      
      if (is.null(healthy_data_raw) || !is.matrix(healthy_data_raw$expression_vectors)) {
        cat(sprintf("  Healthy RAW Data nicht gefunden für %s\n", norm_method))
        next
      }
      
      all_results[[paste(project_id, norm_method, sep = "_")]] <- list(
        gene_results = gene_results,
        family_stats = family_stats,
        healthy_replicates = healthy_data_raw$expression_vectors,
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

# ==================== DISTANCE P-WERTE EXTRAHIEREN ====================

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

# ==================== NOISE P-WERT BERECHNUNG ====================
#' Berechnet Noise-p-Werte nach der kNN-basierten Monte-Carlo Methode
#' @param healthy_replicates Matrix: samples × genes mit Healthy-Expressionen (TPM!)
#' @param gene_results Data.Frame mit Gene-Ergebnissen (enthält observed differences)
#' @param family_stats List mit family_sds pro Gen
#' @param gene_to_fam Vector mit Family-Zuordnung
#' @param norm_method Normalisierungsmethode ("raw", "log", "std_log", "full")
#' @param k Anzahl der nächsten Nachbarn für kNN-Binning (default: 100)
#' @param B Anzahl der Monte-Carlo Iterationen (default: 1000)
#' @return List mit:
#'   - pvalues: Array mit Noise-p-Werten pro (Gen, Stage, Comparison)
#'   - all_null_distances: Vektor mit allen Null-Distanzen für Plot 3
#'   - percentiles: Data.Frame für Plot 4 (observed distance vs percentile)
# ==================== NOISE P-WERT BERECHNUNG ====================
compute_noise_pvalues <- function(healthy_replicates, gene_results, 
                                  family_stats, gene_to_fam,
                                  norm_method = "raw",
                                  k = 100, B = 1000) {
  
  set.seed(42)
  n_genes <- ncol(healthy_replicates)
  n_stages <- length(STAGES)
  n_comparisons <- length(COMP_TYPES)
  
  # Ergebnis-Array
  noise_pvalues <- array(NA, dim = c(n_genes, n_stages, n_comparisons),
                        dimnames = list(
                          colnames(healthy_replicates),
                          STAGES,
                          COMP_TYPES
                        ))
  
  # ===== 1. PRE-LOG SPACE =====
  cat(sprintf("    Berechne pre-log space für %s...\n", norm_method))
  
  if (norm_method %in% c("raw", "log")) {
    values_prelog <- healthy_replicates
    means_prelog <- colMeans(healthy_replicates, na.rm = TRUE)
    names(means_prelog) <- colnames(healthy_replicates)
    
  } else if (norm_method == "std_log") {
    expr_t <- t(healthy_replicates)
    scaled <- tox_normalize_by_std_dev(expr_t)
    values_prelog <- t(scaled)
    means_prelog <- colMeans(values_prelog, na.rm = TRUE)
    names(means_prelog) <- colnames(values_prelog)
    
  } else if (norm_method == "full") {
    expr_t <- t(healthy_replicates)
    scaled <- tox_normalize_by_std_dev(expr_t)
    quantile <- tox_quantile_normalization(scaled)
    values_prelog <- t(quantile)
    means_prelog <- colMeans(values_prelog, na.rm = TRUE)
    names(means_prelog) <- colnames(values_prelog)
  }
  
  # ===== 2. RESIDUEN =====
  cat("    Berechne Residuen im pre-log space...\n")
  residuals_matrix <- tox_compute_residuals(values_prelog, means_prelog)
  
  # ===== 3. OPTIMIERTE NACHBARSCHAFTSSUCHE =====
  cat(sprintf("    Baue Nachbarschaftsstruktur mit k=%d...\n", k))
  
  # Sortiere Gene nach Mittelwert
  gene_order <- order(means_prelog)
  sorted_means <- means_prelog[gene_order]
  sorted_indices <- which(!is.na(means_prelog))[gene_order]
  n_valid <- length(sorted_indices)
  
  # Position-Lookup
  pos_lookup <- integer(n_genes)
  for (i in seq_along(sorted_indices)) {
    pos_lookup[sorted_indices[i]] <- i
  }
  
  # Nachbarschaftssuche
  find_k_nearest <- function(g) {
    pos <- pos_lookup[g]
    if (is.na(pos)) return(integer(0))
    
    window <- min(k * 2, n_valid)
    left <- max(1, pos - window)
    right <- min(n_valid, pos + window)
    
    candidate_positions <- setdiff(left:right, pos)
    if (length(candidate_positions) == 0) return(integer(0))
    
    cand_means <- sorted_means[candidate_positions]
    target_mean <- sorted_means[pos]
    distances <- abs(cand_means - target_mean)
    
    if (length(distances) >= k) {
      selected <- order(distances)[1:k]
      return(sorted_indices[candidate_positions[selected]])
    } else {
      return(sorted_indices[candidate_positions])
    }
  }
  
  # ===== 4. SAMMLE RESIDUEN ALS LISTE =====
  cat("    Sammle Residuen als Liste...\n")
  
  neighbors_list <- vector("list", n_genes)
  for (g in 1:n_genes) {
    if (!is.na(means_prelog[g])) {
      neighbors_list[[g]] <- find_k_nearest(g)
    } else {
      neighbors_list[[g]] <- integer(0)
    }
  }
  
  resid_list <- vector("list", n_genes)
  resid_counts <- integer(n_genes)
  
  for (g in 1:n_genes) {    
    neighbors <- neighbors_list[[g]]
    if (length(neighbors) == 0) next
    
    all_resid <- c()
    for (neighbor in neighbors) {
      all_resid <- c(all_resid, residuals_matrix[, neighbor])
    }
    all_resid <- all_resid[!is.na(all_resid)]
    
    if (length(all_resid) >= 10) {
      resid_list[[g]] <- all_resid
      resid_counts[g] <- length(all_resid)
    }
  }
  
  # residuals_matrix kann jetzt gelöscht werden
  rm(residuals_matrix, values_prelog)
  gc()
  
  # ===== 5. PARTNER-GENE =====
  cat("    Bestimme mögliche Gen-Partner...\n")
  has_noise <- resid_counts >= 10
  genes_with_noise <- which(has_noise)
  
  other_genes <- vector("list", n_genes)
  for (g in genes_with_noise) {
    other_genes[[g]] <- genes_with_noise[genes_with_noise != g]
  }
  
  # Pre-berechnung für log-Daten
  if (norm_method != "raw") {
    log_means_all <- log2(means_prelog + 1)
  }
  
  # ===== 6. STICHPROBE FÜR NULL-DISTANZEN =====
  MAX_NULL_SAMPLE <- 100000  # 100.000 Werte reichen für Plots
  null_sample <- c()
  
  # ===== 7. BATCH-WEISE PARALLELE VERARBEITUNG =====
  cat("    Starte Batch-Verarbeitung...\n")
  
  # Kleinere Chunks
  chunk_size <- 600
  gene_chunks <- split(genes_with_noise, 
                       ceiling(seq_along(genes_with_noise) / chunk_size))
  
  cat(sprintf("    %d Chunks mit je ~%d Genen\n", length(gene_chunks), chunk_size))
  
  # Anzahl Worker pro Batch
  n_cores_per_batch <- min(8, length(gene_chunks))
  
  # Anzahl Batches
  batch_size <- 8  # Chunks pro Batch
  n_batches <- ceiling(length(gene_chunks) / batch_size)
  
  cat(sprintf("    %d Batches mit je %d Chunks\n", n_batches, batch_size))
  
  # Batch-Schleife
  for (batch_idx in 1:n_batches) {
    batch_start <- (batch_idx - 1) * batch_size + 1
    batch_end <- min(batch_idx * batch_size, length(gene_chunks))
    batch_chunks <- gene_chunks[batch_start:batch_end]
    
    # ===== PARALLELE BERECHNUNG FÜR DIESEN BATCH =====
    library(parallel)
    
    batch_results <- mclapply(seq_along(batch_chunks), function(chunk_idx_in_batch) {
      set.seed(42)
      chunk_genes <- batch_chunks[[chunk_idx_in_batch]]
      actual_chunk_idx <- batch_start + chunk_idx_in_batch - 1
            
      # Pre-allocate für Chunk-Ergebnisse
      max_null_per_chunk <- length(chunk_genes) * 500 * B
      chunk_null_dists <- numeric(max_null_per_chunk)
      chunk_null_pos <- 1
      
      # Array für p-Werte
      chunk_pvalues <- array(NA, dim = c(length(chunk_genes), n_stages, n_comparisons))
      
      for (local_idx in seq_along(chunk_genes)) {
        g <- chunk_genes[local_idx]
        
        n_resid <- resid_counts[g]
        resid_g <- resid_list[[g]]
        
        partners <- other_genes[[g]]
        if (length(partners) == 0) next
        
        n_partner_actual <- min(500, length(partners))
        partner_genes <- sample(partners, n_partner_actual)
        
        # Pre-allocate für dieses Gen
        max_dists_gen <- n_partner_actual * B
        gen_null_dists <- numeric(max_dists_gen)
        gen_pos <- 1
        
        # Pre-generiere Zufallsindizes für g
        idx_g_vec <- sample(1:n_resid, n_partner_actual * B, replace = TRUE)
        
        if (norm_method != "raw") {
          log_mean_g <- log_means_all[g]
        }
        
        for (partner_idx in 1:n_partner_actual) {
          g2 <- partner_genes[partner_idx]
          
          n_resid2 <- resid_counts[g2]
          resid_g2 <- resid_list[[g2]]
          
          idx_g2 <- sample(1:n_resid2, B, replace = TRUE)
          
          start_idx <- (partner_idx - 1) * B + 1
          end_idx <- partner_idx * B
          eps1 <- resid_g[idx_g_vec[start_idx:end_idx]]
          eps2 <- resid_g2[idx_g2]
          
          x1_prelog <- means_prelog[g] + eps1
          x2_prelog <- means_prelog[g2] + eps2
          x1_prelog[x1_prelog < 0] <- 0
          x2_prelog[x2_prelog < 0] <- 0
          
          if (norm_method == "raw") {
            null_dists <- abs(x1_prelog - x2_prelog)
          } else {
            log_mean_g2 <- log_means_all[g2]
            x1_log <- log2(x1_prelog + 1)
            x2_log <- log2(x2_prelog + 1)
            delta1 <- x1_log - log_mean_g
            delta2 <- x2_log - log_mean_g2
            null_dists <- abs(delta1 - delta2)
          }
          
          fam_id <- gene_to_fam[g]
          family_sd <- family_stats$gene_smoothed_sds_all[fam_id]
          
          if (!is.na(family_sd) && family_sd > 0) {
            scaled_null <- null_dists / family_sd
            n_new <- length(scaled_null)
            gen_null_dists[gen_pos:(gen_pos + n_new - 1)] <- scaled_null
            gen_pos <- gen_pos + n_new
          }
        }
        
        if (gen_pos > 1) {
          gen_null_dists <- gen_null_dists[1:(gen_pos - 1)]
        } else {
          next
        }
        
        # Füge zu Chunk-Null-Distanzen hinzu
        n_gen_dists <- length(gen_null_dists)
        chunk_null_dists[chunk_null_pos:(chunk_null_pos + n_gen_dists - 1)] <- gen_null_dists
        chunk_null_pos <- chunk_null_pos + n_gen_dists
        
        # p-Wert-Berechnung
        for (s in 1:n_stages) {
          stage <- STAGES[s]
          stage_col <- gsub(" ", "_", stage)
          
          for (comp_idx in 1:n_comparisons) {
            comp_type <- COMP_TYPES[comp_idx]
            obs_col <- switch(comp_type,
              "own_healthy" = paste0("shift_vs_own_healthy_", stage_col),
              "family_mean" = paste0("shift_vs_family_mean_", stage_col),
              "ortholog_mean" = paste0("shift_vs_ortholog_mean_", stage_col)
            )
            
            if (!obs_col %in% colnames(gene_results)) next
            
            observed <- gene_results[g, obs_col]
            if (is.na(observed)) next
            
            d_obs <- abs(observed)
            p_value <- (1 + sum(gen_null_dists >= d_obs)) / (length(gen_null_dists) + 1)
            
            chunk_pvalues[local_idx, s, comp_idx] <- p_value
          }
        }
        
        # Explizite GC alle 10 Gene
        if (local_idx %% 10 == 0) gc()
      }
      
      # Trimme Chunk-Null-Distanzen
      if (chunk_null_pos > 1) {
        chunk_null_dists <- chunk_null_dists[1:(chunk_null_pos - 1)]
      } else {
        chunk_null_dists <- numeric(0)
      }
      
      list(
        null_dists = chunk_null_dists,
        pvalues = chunk_pvalues,
        genes = chunk_genes
      )
      
    }, mc.cores = n_cores_per_batch, mc.preschedule = TRUE)
    
    # ===== 8. ERGEBNISSE DIESES BATCHES VERARBEITEN =====
    cat("    Verarbeite Batch-Ergebnisse...\n")
    
    # p-Werte in Haupt-Array einfügen
    for (chunk_res in batch_results) {
      if (is.null(chunk_res)) next
      
      for (local_idx in seq_along(chunk_res$genes)) {
        g <- chunk_res$genes[local_idx]
        for (s in 1:n_stages) {
          for (comp_idx in 1:n_comparisons) {
            p_val <- chunk_res$pvalues[local_idx, s, comp_idx]
            if (!is.na(p_val)) {
              noise_pvalues[g, s, COMP_TYPES[comp_idx]] <- p_val
            }
          }
        }
      }
    }
    
    # ===== 9. STICHPROBE DER NULL-DISTANZEN SAMMELN =====
    for (chunk_res in batch_results) {
      if (!is.null(chunk_res) && length(chunk_res$null_dists) > 0) {
        if (length(null_sample) < MAX_NULL_SAMPLE) {
          n_needed <- MAX_NULL_SAMPLE - length(null_sample)
          n_take <- min(length(chunk_res$null_dists), n_needed)
          null_sample <- c(null_sample, sample(chunk_res$null_dists, n_take))
        }
      }
    }
    
    # ===== 10. BATCH-DATEN FREIGEBEN =====
    rm(batch_results)
    gc()
  }
  
  # ===== 11. FINALE STICHPROBE =====
  cat("\n  Finale Stichprobe:\n")
  cat(sprintf("    %d Null-Distanzen gesammelt (von möglichen Milliarden)\n", 
              length(null_sample)))
  
  return(list(
    pvalues = noise_pvalues,
    all_null_distances = null_sample  # Nur Stichprobe für Plots
  ))
}

# ==================== GESAMTTABELLE ====================

create_all_genes_table <- function(all_results, k = 100, B = 1000) {
  
  cat("\n>>> Schätze Gesamtgröße für Pre-allocation...\n")
  
  total_estimated <- 0
  for (res_name in names(all_results)) {
    res <- all_results[[res_name]]
    n_genes <- nrow(res$gene_results)
    total_estimated <- total_estimated + (n_genes * length(STAGES) * length(COMP_TYPES))
  }
  
  buffer <- 1.1
  alloc_rows <- ceiling(total_estimated * buffer)
  
  cat(sprintf("  Geschätzte Gesamtzeilen (max): %d\n", total_estimated))
  cat(sprintf("  Pre-allocated Liste: %d Zeilen\n", alloc_rows))
  
  # ===== 1. PRE-ALLOCATE ALS LISTE VON LISTEN (viel schneller!) =====
  all_genes_list <- vector("list", alloc_rows)
  global_idx <- 1
  
  # ===== 2. NULL-DISTANZEN: Liste statt rbindlist =====
  null_dist_list <- vector("list", length(all_results))
  null_idx <- 1
  
  # ===== 3. PERZENTILE: Auch als Liste =====
  percentiles_list <- vector("list", length(all_results) * 10)  # Puffer
  pct_idx <- 1
  
  for (res_name in names(all_results)) {
    res <- all_results[[res_name]]
    cat(sprintf("\n  Verarbeite %s...\n", res_name))
    
    family_n_genes <- tabulate(res$gene_to_fam, nbins = res$n_families)
    
    cat("    Extrahiere Distance p-Werte...\n")
    distance_p <- extract_distance_pvalues(res$gene_results)
    
    cat(sprintf("    Berechne Noise p-Werte mit k=%d, B=%d...\n", k, B))
    noise_result <- compute_noise_pvalues(
      healthy_replicates = res$healthy_replicates,
      gene_results = res$gene_results,
      family_stats = res$family_stats,
      gene_to_fam = res$gene_to_fam,
      norm_method = res$norm_method,
      k = k,
      B = B
    )
    
    noise_p <- noise_result$pvalues
    percentiles <- noise_result$percentiles
    
    # ===== 4. NULL-DISTANZEN: In Liste sammeln  =====
    if (length(noise_result$all_null_distances) > 0) {
      null_dist_list[[null_idx]] <- noise_result$all_null_distances
      null_idx <- null_idx + 1
    }
    
    # ===== 5. PERZENTILE: Auch in Liste sammeln =====
    if (!is.null(percentiles) && nrow(percentiles) > 0) {
      # Füge Metadaten in einem Schritt hinzu
      percentiles$cancer_id <- res$cancer_id
      percentiles$cancer_type <- res$cancer_type
      percentiles$norm_method <- res$norm_method
      percentiles_list[[pct_idx]] <- percentiles
      pct_idx <- pct_idx + 1
    }
    
    n_genes <- nrow(res$gene_results)
    rows_added <- 0
    
    # ===== 6. OPTIMIERTE DATENEXTRAKTION =====
    # Sammle alle Daten in Matrizen statt Data.Frames
    gene_ids_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    cancer_types_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    cancer_ids_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    stages_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    normalizations_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    norm_methods_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    comparisons_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    directions_all <- character(n_genes * length(STAGES) * length(COMP_TYPES))
    signed_diffs_all <- numeric(n_genes * length(STAGES) * length(COMP_TYPES))
    dist_p_all <- numeric(n_genes * length(STAGES) * length(COMP_TYPES))
    noise_p_all <- numeric(n_genes * length(STAGES) * length(COMP_TYPES))
    family_sizes_all <- integer(n_genes * length(STAGES) * length(COMP_TYPES))
    family_means_all <- numeric(n_genes * length(STAGES) * length(COMP_TYPES))
    
    local_idx <- 1
    
    for (g in 1:n_genes) {      
      gene_id <- res$gene_results$gene_id[g]
      
      for (s in 1:length(STAGES)) {
        stage <- STAGES[s]
        stage_col <- gsub(" ", "_", stage)
        
        for (comp_type in COMP_TYPES) {
          
          dist_p_val <- distance_p[g, s, comp_type]
          noise_p_val <- noise_p[g, s, comp_type]
          
          if (is.na(dist_p_val) || is.na(noise_p_val)) next
          
          shift_col <- switch(comp_type,
            "own_healthy" = paste0("shift_vs_own_healthy_", stage_col),
            "family_mean" = paste0("shift_vs_family_mean_", stage_col),
            "ortholog_mean" = paste0("shift_vs_ortholog_mean_", stage_col)
          )
          
          if (!shift_col %in% colnames(res$gene_results)) next
          
          shift <- res$gene_results[g, shift_col]
          direction <- ifelse(is.na(shift), NA, ifelse(shift > 0, "up", "down"))
          
          fam_id <- res$gene_to_fam[g]
          family_size <- if (!is.na(fam_id) && fam_id > 0) family_n_genes[fam_id] else NA
          family_mean_val <- if (!is.na(fam_id) && fam_id > 0) 
            res$family_stats$gene_family_means_all[fam_id] else NA
          
          # In Matrizen speichern (viel schneller als Data.Frame!)
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
    }
    
    # ===== 7. EINEN Data.Frame pro Projekt erstellen =====
    if (rows_added > 0) {
      project_df <- data.frame(
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
      
      all_genes_list[[global_idx]] <- project_df
      global_idx <- global_idx + 1
    }
    
    cat(sprintf("    → %d Zeilen für %s\n", rows_added, res_name))
  }
  
  cat(sprintf("\n>>> Gesamt: %d Zeilen generiert\n", global_idx - 1))
  
  # ===== 8. KOMBINIERE ALLE PROJEKTE =====
  all_genes_list <- all_genes_list[1:(global_idx - 1)]
  result <- rbindlist(all_genes_list)  # data.table ist hier viel schneller!
  
  # ===== 9. NULL-DISTANZEN SPEICHERN =====
  null_dist_list <- null_dist_list[1:(null_idx - 1)]
  all_null_distances <- unlist(null_dist_list)
  
  null_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "diagnostic_plots", "null_distances")
  dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(all_null_distances, 
          file.path(null_dir, sprintf("global_null_distances_k%d_B%d.rds", k, B)))
  
  # ===== 10. PERZENTILE SPEICHERN =====
  percentiles_list <- percentiles_list[1:(pct_idx - 1)]
  if (length(percentiles_list) > 0) {
    all_percentiles <- rbindlist(percentiles_list, fill = TRUE)
    percentiles_file <- file.path(dirname(FAMILY_OUTPUT_DIR), "diagnostic_plots", 
                                   sprintf("percentiles_k%d_B%d.rds", k, B))
    saveRDS(all_percentiles, percentiles_file)
    cat(sprintf("\n  Perzentil-Daten gespeichert: %d Einträge\n", nrow(all_percentiles)))
  }
  
  cat(sprintf("\n FINALE TABELLE: %d Zeilen\n", nrow(result)))
  return(result)
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
analyze_outlier_significance <- function(k = 100, B = 1000, output_dir = NULL) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat(sprintf("OUTLIER SIGNIFICANCE ANALYSIS (k=%d, B=%d)\n", k, B))
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  all_results <- load_all_gene_results()
  if (length(all_results) == 0) stop("Keine Daten gefunden!")
  
  if (is.null(output_dir)) {
    output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  param_suffix <- sprintf("_k%d_B%d", k, B)
  
  # ===== 1. GESAMTTABELLE =====
  cat("\n>>> Erstelle Gesamttabelle...\n")
  all_genes_table <- create_all_genes_table(all_results, k = k, B = B)
  
  # ===== 2. SPEICHERE ROHDATEN (UNVERÄNDERTE DATEINAMEN) =====
  saveRDS(all_genes_table, file.path(output_dir, paste0("all_genes_pvalues", param_suffix, ".rds")))
  write.csv(all_genes_table, 
            file.path(output_dir, paste0("all_genes_pvalues", param_suffix, ".csv")), 
            row.names = FALSE)
  
  # ===== 3. BASELINE AUSGABE =====
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("UNKORRIGIERTE WERTE (BASELINE)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  cat("\n--- Summary distance_p_values ---\n")
  print(summary(all_genes_table$distance_p_value))
  cat("\n--- Summary noise_p_values ---\n")
  print(summary(all_genes_table$noise_p_value))
  
  # ===== 4. NEUE SELEKTION: NUR NOISE-P-WERT < 0.05, SORTIERT NACH DISTANCE =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("SELEKTION: noise_p_adj < 0.05 (sortiert nach distance_p_value)\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  ALPHA_NOISE <- 0.05
  
  significant_by_norm <- list()
  
  for (current_norm in unique(all_genes_table$norm_method)) {
    
    cat(sprintf("\n>>> Verarbeite Normalisierung: %s <<<\n", current_norm))
    
    norm_data <- all_genes_table %>% 
      filter(norm_method == current_norm)
    
    # BH-Korrektur für noise_p_value (pro cancer_id)
    norm_data <- norm_data %>%
      group_by(cancer_id) %>%
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
    
    # ===== NEUE SELEKTION: Nur noise_p_adj < 0.05 =====
    # (distance_p_value wird nicht gefiltert, nur für Ranking verwendet)
    significant <- norm_data %>%
      filter(noise_p_value_adj < ALPHA_NOISE) %>%
      arrange(distance_p_value)  # Sortiere nach distance p-Wert (kleinste zuerst = höchste Signifikanz)
    
    cat(sprintf("\n  → Signifikante Einträge (noise_p_adj < 0.05): %d\n", nrow(significant)))
    
    if (nrow(significant) > 0) {
      # Top 10 anzeigen
      cat("\n  Top 10 nach distance_p_value:\n")
      print(significant %>%
        select(gene_id, stage, comparison, distance_p_value, noise_p_value_adj) %>%
        head(10))
      
      # ===== SPEICHERN MIT GLEICHEN DATEINAMEN (ÜBERSCHREIBT VORHERIGE) =====
      # Wichtig: Gleicher Dateiname wie vorher, damit downstream Skripte funktionieren
      filename <- sprintf("significant_%s%s.csv", 
                         gsub("[^a-zA-Z0-9]", "_", current_norm),
                         param_suffix)
      write.csv(significant, file.path(output_dir, filename), row.names = FALSE)
      
      significant_by_norm[[current_norm]] <- significant
    }
  }
  
  # ===== 5. GESAMTTABELLE ÜBER ALLE NORMALISIERUNGEN =====
  cat("\n>>> Gesamttabelle über alle Normalisierungen <<<\n")
  
  if (length(significant_by_norm) > 0) {
    all_significant <- rbindlist(significant_by_norm, fill = TRUE)
    
    # Gesamt nach distance_p_value sortieren
    all_significant <- all_significant %>%
      arrange(distance_p_value)
    
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
    
    # ===== SPEICHERN MIT GLEICHEM DATEINAMEN =====
    # Wichtig: Gleicher Dateiname wie vorher (significant_all...)
    write.csv(all_significant, 
              file.path(output_dir, paste0("significant_all", param_suffix, ".csv")), 
              row.names = FALSE)
    saveRDS(all_significant, 
            file.path(output_dir, paste0("significant_all", param_suffix, ".rds")))
  } else {
    cat("\n  Keine signifikanten Einträge gefunden!\n")
    all_significant <- data.frame()
  }
  
  # ===== 6. JACCARD PLOTS (basierend auf noise-filterten Daten) =====
  # Wichtig: Gleiche Plot-Namen wie vorher
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
          # ===== PLOT 1: JACCARD FÜR KORRIGIERTE DATEN (GLEICHER NAME) =====
          p_adj <- plot_jaccard_heatmap(
            plot_data,
            title = sprintf("%s (min %d, k=%d, B=%d) - KORRIGIERT", 
                            comp_type, min_size, k, B),
            subtitle = sprintf("noise_p_adj < 0.05, n=%d", nrow(plot_data))
          )
          
          filename_adj <- sprintf("jaccard_ADJ_%s_min%d%s.png", 
                                 comp_type, min_size, param_suffix)
          ggsave(file.path(plots_dir, filename_adj), p_adj, 
                 width = 7, height = 6, dpi = 300)
          
          # ===== PLOT 2: VERGLEICH ROHDATEN VS KORRIGIERT (OPTIONAL) =====
          # Rohdaten für Vergleich (mit noise_p < 0.05, unkorrigiert)
          raw_data <- all_genes_table %>%
            filter(
              comparison == comp_type,
              family_size >= min_size,
              distance_p_value < 0.05,
              noise_p_value < 0.05
            )
          
          if (nrow(raw_data) > 0 && nrow(plot_data) > 0) {
            # Jaccard zwischen Rohdaten und korrigiert
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
            
            filename_comp <- sprintf("jaccard_COMP_%s_min%d%s.png", 
                                    comp_type, min_size, param_suffix)
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
    significant = all_significant,
    params = list(k = k, B = B)
  )))
}

# ==================== AUSFÜHRUNG ====================

if (sys.nframe() == 0) {
  set.seed(42)
  base_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_conservative_grouping")
  
  # for now
  k_values <- c(25, 50, 75, 100)
  B_values <- c(1000, 750, 500)
  
  for (k in k_values) {
    for (B in B_values) {
      cat("\n", paste(rep("#", 80), collapse = ""), "\n")
      cat(sprintf("### TEST: k=%d, B=%d ###\n", k, B))
      cat(paste(rep("#", 80), collapse = ""), "\n")
      
      param_dir <- file.path(base_dir, sprintf("k%d_B%d", k, B))
      dir.create(param_dir, recursive = TRUE, showWarnings = FALSE)
      
      results <- analyze_outlier_significance(k = k, B = B, output_dir = param_dir)
      
      saveRDS(list(
        k = k,
        B = B,
        timestamp = Sys.time(),
        n_significant = nrow(results$significant)
      ), file.path(param_dir, "metadata.rds"))
    }
  }
}
