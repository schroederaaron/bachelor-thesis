#!/usr/bin/env Rscript
# pathway_analysis_parallel.R
# Parallelisierte Pathway-Analyse mit 16 Kernen

library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(patchwork)
library(ggwordcloud)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(future)
library(future.apply)

source("config.R")
source("utils.R")

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

NORM_METHODS <- c("raw", "log", "std_log", "full")
COMP_TYPES <- c("own_healthy", "family_mean", "ortholog_mean")
DIRECTIONS <- c("all", "up", "down")

TOP_N_GENES <- 200
P_CUTOFF <- 0.05
Q_CUTOFF <- 0.2

# Parallelisierung
N_CORES <- 12
options(future.globals.maxSize = 50 * 1024^3)  # 50 GB

# ==================== CACHING ====================

# Globaler Cache für Pathway-Ergebnisse
pathway_cache <- list()

#' Führe GO-Analyse mit Caching durch
cached_enrichGO <- function(genes, cache_key) {
  cache_key_go <- paste0(cache_key, "_GO")
  
  if (!is.null(pathway_cache[[cache_key_go]])) {
    return(pathway_cache[[cache_key_go]])
  }
  
  result <- tryCatch({
    enrichGO(
      gene = genes,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = P_CUTOFF,
      qvalueCutoff = Q_CUTOFF,
      readable = TRUE
    )
  }, error = function(e) NULL)
  
  pathway_cache[[cache_key_go]] <<- result
  return(result)
}

#' Führe KEGG-Analyse mit Caching durch
cached_enrichKEGG <- function(genes, cache_key) {
  cache_key_kegg <- paste0(cache_key, "_KEGG")
  
  if (!is.null(pathway_cache[[cache_key_kegg]])) {
    return(pathway_cache[[cache_key_kegg]])
  }
  
  result <- tryCatch({
    enrichKEGG(
      gene = genes,
      organism = "hsa",
      pvalueCutoff = P_CUTOFF,
      qvalueCutoff = Q_CUTOFF
    )
  }, error = function(e) NULL)
  
  pathway_cache[[cache_key_kegg]] <<- result
  return(result)
}

# ==================== MAPPING ====================

#' Erstelle Mapping mit Lookup-Table
create_mapping_lookup <- function(mapping_df) {
  mapping_vec <- mapping_df$entrez_id
  names(mapping_vec) <- mapping_df$uniprot_id
  return(mapping_vec)
}

# ==================== PLOT-FUNKTIONEN ====================

plot_go_wordcloud <- function(ego, title, filename, max_words = 15) {
  
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  
  set.seed(42)  # für reproduzierbare Platzierung
  
  go_df <- as.data.frame(ego) %>%
    arrange(pvalue) %>%
    head(max_words) %>%
    mutate(
      log_p = -log10(pvalue),
      weight = log_p * Count,
      angle = 0  # Feste horizontale Ausrichtung
    )
  
  # Skaliere Gewichte für bessere Darstellung
  weight_range <- range(go_df$weight, na.rm = TRUE)
  if (diff(weight_range) < 1e-6) {
    go_df$size <- 8
  } else {
    go_df$size <- 4 + 12 * (go_df$weight - weight_range[1]) / diff(weight_range)
  }
  
  p <- ggplot(go_df, aes(label = Description, size = size, color = log_p)) +
    geom_text_wordcloud(
      shape = "circle",
      eccentricity = 1,
      rm_outside = TRUE,
      grid_size = 2,
      grid_margin = 2,
      area_corr_power = 1,
      placement = "williams"
    ) +
    scale_size_identity() +
    scale_color_gradient(low = "darkblue", high = "red") +
    theme_minimal() +
    labs(
      title = title,
      subtitle = sprintf("Top %d GO Biological Process Terms", max_words)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10)
    )
  
  ggsave(filename, p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

plot_go_barplot <- function(ego, title, filename, max_terms = 20) {
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  
  go_df <- as.data.frame(ego) %>%
    arrange(pvalue) %>%
    head(max_terms) %>%
    mutate(
      log_p = -log10(pvalue),
      short_desc = ifelse(nchar(Description) > 50, 
                          paste0(substr(Description, 1, 47), "..."), 
                          Description),
      short_desc = factor(short_desc, levels = rev(short_desc))
    )
  
  p <- ggplot(go_df, aes(x = short_desc, y = log_p, fill = Count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(title = title, x = "", y = "-log10(p-value)", fill = "Gen-Anzahl") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.text.y = element_text(size = 9))
  
  ggsave(filename, p, width = 12, height = max(6, max_terms * 0.4), dpi = 300)
  return(p)
}

plot_go_dotplot <- function(ego, title, filename, max_terms = 20) {
  if (is.null(ego) || nrow(ego) == 0) return(NULL)
  
  p <- dotplot(ego, showCategory = max_terms) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(filename, p, width = 12, height = max(6, max_terms * 0.4), dpi = 300)
  return(p)
}

plot_kegg_dotplot <- function(kk, title, filename, max_terms = 20) {
  if (is.null(kk) || nrow(kk) == 0) return(NULL)
  
  p <- dotplot(kk, showCategory = max_terms) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  ggsave(filename, p, width = 12, height = max(6, max_terms * 0.4), dpi = 300)
  return(p)
}

plot_top_genes_barplot <- function(gene_summary, title, filename, n_top = 30) {
  plot_data <- gene_summary %>%
    head(n_top) %>%
    mutate(
      gene_label = ifelse(nchar(gene_id) > 12, 
                          paste0(substr(gene_id, 1, 10), "..."),
                          gene_id),
      direction = factor(direction, levels = c("up", "down"))
    ) %>%
    arrange(desc(mean_distance))
  
  p <- ggplot(plot_data, aes(x = reorder(gene_label, mean_distance), 
                              y = mean_distance, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
    labs(title = title, x = "Gen", y = "Mittlere |signed difference|") +
    theme_minimal() +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
  
  ggsave(filename, p, width = 10, height = max(6, n_top * 0.3), dpi = 300)
  return(p)
}

# ==================== PARALLELE PATHWAY-ANALYSE ====================

#' Verarbeite eine einzelne Gruppe (für parallele Ausführung)
process_group <- function(group_info, mapping_lookup, base_dir) {

  cancer_id <- group_info$cancer_id
  norm_method <- group_info$norm_method
  comparison <- group_info$comparison
  direction <- group_info$direction
  gene_ids <- group_info$gene_ids
  cache_key <- group_info$cache_key
  
  # Erstelle Ausgabeverzeichnis
  out_dir <- file.path(base_dir, cancer_id, norm_method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Mappe zu Entrez-IDs
  entrez_ids <- mapping_lookup[gene_ids]
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  if (length(entrez_ids) < 5) {
    return(list(
      cache_key = cache_key,
      status = "skipped",
      reason = "too_few_mapped",
      n_mapped = length(entrez_ids)
    ))
  }
  
  results <- list(
    cache_key = cache_key,
    status = "success",
    n_genes = length(gene_ids),
    n_mapped = length(entrez_ids)
  )
  
  # GO-Analyse
  ego <- cached_enrichGO(entrez_ids, cache_key)
  if (!is.null(ego) && nrow(ego) > 0) {
    suffix <- paste(comparison, direction, sep = "_")
    
    tryCatch({
      p <- plot_go_wordcloud(ego, 
              title = sprintf("GO - %s", cache_key),
              filename = file.path(out_dir, sprintf("go_wordcloud_%s.png", suffix)))
      results$go_wordcloud <- TRUE
    }, error = function(e) {})
    
    tryCatch({
      p <- plot_go_barplot(ego, 
              title = sprintf("GO - %s", cache_key),
              filename = file.path(out_dir, sprintf("go_barplot_%s.png", suffix)))
      results$go_barplot <- TRUE
    }, error = function(e) {})
    
    tryCatch({
      p <- plot_go_dotplot(ego, 
              title = sprintf("GO - %s", cache_key),
              filename = file.path(out_dir, sprintf("go_dotplot_%s.png", suffix)))
      results$go_dotplot <- TRUE
    }, error = function(e) {})
  }
  
  # KEGG-Analyse
  kk <- cached_enrichKEGG(entrez_ids, cache_key)
  if (!is.null(kk) && nrow(kk) > 0) {
    suffix <- paste(comparison, direction, sep = "_")
    tryCatch({
      p <- plot_kegg_dotplot(kk, 
              title = sprintf("KEGG - %s", cache_key),
              filename = file.path(out_dir, sprintf("kegg_%s.png", suffix)))
      results$kegg <- TRUE
    }, error = function(e) {})
  }
  
  return(results)
}

# ==================== HAUPTFUNKTION ====================

run_pathway_analysis_parallel <- function(output_dir) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("PARALLELIZED PATHWAY ANALYSIS\n")
  cat(sprintf("Cores=%d\n", N_CORES))
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  start_time <- Sys.time()
  
  # ===== 1. AUSGABEVERZEICHNIS =====
  base_plots_dir <- file.path(output_dir, "pathway_plots")
  dir.create(base_plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===== 2. LADE MAPPING =====
  mapping_file <- file.path(output_dir, "uniprot_entrez_mapping.rds")
  if (!file.exists(mapping_file)) {
    stop("Mapping-Datei nicht gefunden. Bitte erstellen Sie sie zuerst.")
  }
  mapping <- readRDS(mapping_file)
  mapping_lookup <- create_mapping_lookup(mapping)
  cat(sprintf("  Mapping geladen: %d Einträge\n", nrow(mapping)))
  
  # ===== 3. LADE DATEN =====
  cat("\n  Lade TensorOmics Daten...\n")
  sig_dir <- file.path(dirname(FAMILY_OUTPUT_DIR), "outlier_significance_adaptive_knn_fortran")
  data_file <- file.path(sig_dir, paste0("significant_all.rds"))
  
  if (!file.exists(data_file)) {
    cat("significance file not found, using all_genes_pvalues\n")
    data_file <- file.path(sig_dir, paste0("all_genes_pvalues.rds"))
    all_data <- readRDS(data_file)
    significant <- all_data %>% filter(noise_p_value_adj < 0.01)
  } else {
    significant <- readRDS(data_file)
  }
  
  # ===== 4. EXTRAHIERE TOP-GENE =====
  cat("\n  Extrahiere Top-Gene...\n")
  
  gene_summaries <- significant %>%
    group_by(cancer_id, norm_method, comparison, gene_id, direction) %>%
    summarise(
      mean_distance = mean(abs(signed_difference), na.rm = TRUE),
      n_stages = n(),
      .groups = "drop"
    ) %>%
    arrange(cancer_id, norm_method, comparison, desc(mean_distance))
  
  # Top N pro Gruppe
  top_genes <- gene_summaries %>%
    group_by(cancer_id, norm_method, comparison, direction) %>%
    slice_head(n = TOP_N_GENES) %>%
    ungroup()
  
  top_genes_all <- gene_summaries %>%
    group_by(cancer_id, norm_method, comparison) %>%
    slice_head(n = TOP_N_GENES) %>%
    ungroup() %>%
    mutate(direction = "all")
  
  all_top_genes <- bind_rows(top_genes, top_genes_all)
  cat(sprintf("  → %d Gen-Kombinationen extrahiert\n", nrow(all_top_genes)))
  
  # ===== 5. GRUPPIERE FÜR PARALLELE VERARBEITUNG =====
  groups <- all_top_genes %>%
    group_by(cancer_id, norm_method, comparison, direction) %>%
        summarise(
            gene_ids = list(gene_id),
            cache_key = paste(unique(cancer_id), unique(norm_method), 
                            unique(comparison), unique(direction), sep = "_"),
            .groups = "drop"
        )
  
  cat(sprintf("\n  → %d Gruppen für Pathway-Analyse\n", nrow(groups)))
  
  # ===== 6. PARALLELE VERARBEITUNG =====
  cat("\n  Starte parallele Verarbeitung mit", N_CORES, "Kernen...\n")
  
  # Setup für parallele Verarbeitung
  plan(multisession, workers = N_CORES)

  cat("converting to list\n")
  
  # Konvertiere groups zu Liste für future_lapply
  groups_list <- split(groups, seq(nrow(groups)))
  groups_list <- lapply(groups_list, as.list)

  cat("done converting\n")
  
  batch_start <- Sys.time()
  
  results <- future_lapply(groups_list, function(g) {
    process_group(
      group_info = list(
        cancer_id = g$cancer_id,
        norm_method = g$norm_method,
        comparison = g$comparison,
        direction = g$direction,
        gene_ids = unlist(g$gene_ids),
        cache_key = g$cache_key
      ),
      mapping_lookup = mapping_lookup,
      base_dir = base_plots_dir
    )
  }, future.seed = TRUE)
  
  batch_time <- difftime(Sys.time(), batch_start, units = "mins")
  cat(sprintf("\n  ✓ Parallele Verarbeitung: %.1f Minuten\n", batch_time))
  
  # ===== 7. ZUSAMMENFASSUNG =====
  successful <- sum(sapply(results, function(r) r$status == "success"))
  skipped <- sum(sapply(results, function(r) r$status == "skipped"))
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("PATHWAY ANALYSIS COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("  Erfolgreich: %d/%d Gruppen\n", successful, nrow(groups)))
  cat(sprintf("  Übersprungen: %d (zu wenige gemappte Gene)\n", skipped))
  cat(sprintf("  Ergebnisse in: %s\n", base_plots_dir))
  
  total_time <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Gesamtzeit: %.1f Minuten\n", total_time))
  
  return(results)
}

# ==================== AUSFÜHRUNG ====================

if (sys.nframe() == 0) {
  output_dir <- file.path(dirname(FAMILY_OUTPUT_DIR))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  results <- run_pathway_analysis_parallel(output_dir)
}