# QC-Skript für erweiterte TCGA-Daten

library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(scales)

# ==================== KONFIGURATION ====================
data_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/material/TCGA_data"
output_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/results/TCGA_QC_Reports"

# ==================== FUNKTION FÜR ERWEITERTE QC ====================

run_enhanced_qc <- function(project_id, data_dir, output_dir) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("ENHANCED QC ANALYSIS FOR:", project_id, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Schritt 1: Lade erweiterte Daten
  cat("STEP 1: LOADING ENHANCED DATA\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  project_data_dir <- file.path(data_dir, project_id)
  rds_file <- file.path(project_data_dir, paste0(project_id, "_enhanced.rds"))
  
  if (!file.exists(rds_file)) {
    stop("Enhanced data not found for ", project_id, 
         "\nPlease run the data pipeline first.")
  }
  
  enhanced_data <- readRDS(rds_file)
  cat("✓ Loaded enhanced data from:", rds_file, "\n")
  cat("  - Samples:", ncol(enhanced_data), "\n")
  cat("  - Genes:", nrow(enhanced_data), "\n")
  
  # Schritt 2: Extrahiere Daten
  clinical_data <- as.data.frame(colData(enhanced_data))
  
  if ("unstranded" %in% names(assays(enhanced_data))) {
    count_data <- assay(enhanced_data, "unstranded")
  } else {
    count_data <- assay(enhanced_data, 1)
  }
  
  # Schritt 3: Erstelle erweiterten QC-Bericht
  cat("\nSTEP 2: CREATING ENHANCED QC REPORT\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # ==================== ERWEITERTE STATISTIKEN ====================
  
  # Basic statistics
  basic_stats <- list(
    project = project_id,
    total_samples = ncol(count_data),
    total_genes = nrow(count_data),
    tumor_samples = sum(clinical_data$sample_type %in% 
                         c("Primary Tumor", "Metastatic"), na.rm = TRUE),
    normal_samples = sum(clinical_data$sample_type == "Solid Tissue Normal", na.rm = TRUE)
  )
  
  # Expression statistics
  expr_stats <- data.frame(
    sample = colnames(count_data),
    barcode = colnames(count_data),
    raw_total_counts = colSums(count_data),
    raw_genes_expressed = colSums(count_data > 0),
    raw_percent_zeros = colSums(count_data == 0) / nrow(count_data) * 100,
    stringsAsFactors = FALSE
  ) %>%
    left_join(
      clinical_data %>% 
        select(barcode, sample_type, bcr_patient_barcode, 
               final_stage, simple_stage, stage_source),
      by = "barcode"
    )
  
  # ==================== Visualisation ====================
  
  colors <- brewer.pal(8, "Set2")
  
  # ----- PLOT 1: Stage Distribution (Enhanced vs Original) -----
  if ("final_stage" %in% names(clinical_data)) {
    
    # Bereite Daten für Vergleich vor
    if ("ajcc_pathologic_stage" %in% names(clinical_data)) {
      comparison_data <- clinical_data %>%
        filter(sample_type %in% c("Primary Tumor", "Metastatic")) %>%
        select(bcr_patient_barcode, 
               original = ajcc_pathologic_stage, 
               enhanced = final_stage) %>%
        distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
        pivot_longer(cols = -bcr_patient_barcode, 
                    names_to = "source", 
                    values_to = "stage") %>%
        filter(!is.na(stage))
      
      if (nrow(comparison_data) > 0) {
        p1a <- ggplot(comparison_data, aes(x = stage, fill = source)) +
          geom_bar(position = "dodge") +
          scale_fill_manual(values = c("original" = colors[1], 
                                      "enhanced" = colors[2])) +
          labs(title = "Stage Distribution: Original vs Enhanced",
               subtitle = paste("Tumors only | Enhanced adds", 
                               length(unique(comparison_data$stage[comparison_data$source == "enhanced"])) - 
                               length(unique(comparison_data$stage[comparison_data$source == "original"])),
                               "more patients"),
               x = "Stage", y = "Count") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      }
    }
    
    # Enhanced stage distribution
    enhanced_stage_dist <- clinical_data %>%
      filter(sample_type %in% c("Primary Tumor", "Metastatic")) %>%
      group_by(final_stage) %>%
      summarise(n = n(), .groups = 'drop') %>%
      mutate(percentage = round(n / sum(n) * 100, 1))
    
    p1b <- ggplot(enhanced_stage_dist %>% filter(!is.na(final_stage)), 
                  aes(x = reorder(final_stage, -n), y = n, fill = final_stage)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0(n, "\n(", percentage, "%)")), 
                vjust = -0.3, size = 3) +
      scale_fill_hue() +
      labs(title = "Enhanced Stage Distribution",
           subtitle = paste("Using final_stage column |", 
                           sum(enhanced_stage_dist$n), "tumors"),
           x = "Stage", y = "Count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
  }
  
  # ----- PLOT 2: Stage Completion by Source -----
  if ("stage_source" %in% names(clinical_data)) {
    source_summary <- clinical_data %>%
      filter(sample_type %in% c("Primary Tumor", "Metastatic")) %>%
      group_by(stage_source) %>%
      summarise(
        n = n(),
        .groups = 'drop'
      ) %>%
      mutate(
        percentage = n / sum(n) * 100,
        stage_source = factor(stage_source, 
                             levels = c("Expression_Data", "Enhanced_Clinical", "Missing"))
      )
    
    p2 <- ggplot(source_summary, aes(x = stage_source, y = n, fill = stage_source)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")), 
                vjust = -0.3, size = 4) +
      scale_fill_manual(values = c("Expression_Data" = colors[3], 
                                  "Enhanced_Clinical" = colors[4],
                                  "Missing" = colors[5])) +
      labs(title = "Stage Information Sources",
           subtitle = "Where did the stage information come from?",
           x = "Source", y = "Count") +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  # ----- PLOT 3: Enhanced Clinical Overview -----
  clinical_summary <- data.frame(
    Variable = c("Tumor Samples", "Normal Samples", 
                "Patients with Stage", "Patients Missing Stage"),
    Count = c(
      basic_stats$tumor_samples,
      basic_stats$normal_samples,
      sum(!is.na(clinical_data$final_stage[clinical_data$sample_type %in% 
                                           c("Primary Tumor", "Metastatic")])),
      sum(is.na(clinical_data$final_stage[clinical_data$sample_type %in% 
                                         c("Primary Tumor", "Metastatic")]))
    )
  )

  clinical_summary <- clinical_summary %>%
    mutate(
      Percentage = round(Count / sum(clinical_summary$Count[1:2]) * 100, 1),
      Variable = factor(Variable, levels = Variable)
    )
  
  p3 <- ggplot(clinical_summary, aes(x = Variable, y = Count, fill = Variable)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(Count, "\n(", Percentage, "%)")), 
              vjust = -0.3, size = 3.5) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Enhanced Clinical Overview",
         subtitle = paste("Project:", project_id),
         x = "", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  # ----- PLOT 4: Quality Metrics by Tissue Type -----
  p4 <- ggplot(expr_stats, aes(x = sample_type, y = raw_total_counts/1e6, fill = sample_type)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
    scale_fill_manual(values = c("Primary Tumor" = colors[2], 
                                "Solid Tissue Normal" = colors[1],
                                "Metastatic" = colors[3])) +
    labs(title = "Library Size by Sample Type",
         subtitle = "Enhanced clinical data",
         x = "Sample Type", y = "Total Reads (Millions)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # ----- PLOT 5: Genes Expressed by Sample Type -----
  p5 <- ggplot(expr_stats, aes(x = sample_type, y = raw_genes_expressed/1000, fill = sample_type)) +
    geom_violin(alpha = 0.6) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.3) +
    scale_fill_manual(values = c("Primary Tumor" = colors[2], 
                                "Solid Tissue Normal" = colors[1],
                                "Metastatic" = colors[3])) +
    labs(title = "Genes Expressed by Sample Type",
         subtitle = "Genes with counts > 0",
         x = "Sample Type", y = "Genes Expressed (Thousands)") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # ----- PLOT 6: Stage vs Quality Metrics -----
  if ("final_stage" %in% names(expr_stats)) {
    stage_metrics <- expr_stats %>%
      filter(sample_type %in% c("Primary Tumor", "Metastatic"),
             !is.na(final_stage))
    
    if (nrow(stage_metrics) > 5) {
      p6 <- ggplot(stage_metrics, 
                   aes(x = final_stage, y = raw_total_counts/1e6, fill = final_stage)) +
        geom_boxplot() +
        geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
        stat_summary(fun = median, geom = "point", 
                    shape = 23, size = 3, fill = "red") +
        labs(title = "Library Size by Stage",
             subtitle = "Only tumors with stage information",
             x = "Stage", y = "Total Reads (Millions)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")
    }
  }
  
  # ==================== COMBINE PLOTS ====================
  
  # Erstelle Layout basierend auf verfügbaren Plots
  plot_list <- list(p3, p4, p5)
  
  if (exists("p1b")) plot_list <- c(list(p1b), plot_list)
  if (exists("p2")) plot_list <- c(plot_list, list(p2))
  if (exists("p6")) plot_list <- c(plot_list, list(p6))
  
  # Kombiniere Plots
  combined_plot <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = paste("Enhanced QC Report -", project_id),
      subtitle = paste("Generated:", Sys.Date(), 
                      "| Samples:", basic_stats$total_samples,
                      "| Tumors:", basic_stats$tumor_samples,
                      "| Normals:", basic_stats$normal_samples),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
    )
  
  # ==================== SPEICHERN ====================
  
  cat("\nSTEP 3: SAVING RESULTS\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Erstelle Ausgabeordner
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Speichere Plot
  plot_file <- file.path(output_dir, paste0(project_id, "_enhanced_QC.png"))
  ggsave(plot_file, combined_plot, 
         width = 16, height = 20, dpi = 300, bg = "white")
  cat("✓ Saved plot to:", plot_file, "\n")
  
  # Speichere QC-Daten
  qc_data_file <- file.path(output_dir, paste0(project_id, "_enhanced_QC.rds"))
  saveRDS(list(
    enhanced_data = enhanced_data,
    clinical_data = clinical_data,
    expr_stats = expr_stats,
    basic_stats = basic_stats,
    plot = combined_plot
  ), file = qc_data_file)
  cat("✓ Saved QC data to:", qc_data_file, "\n")
  
  # ==================== ZUSAMMENFASSUNG ====================
  
  cat("\n", paste(rep("-", 80), collapse = ""), "\n")
  cat("ENHANCED QC SUMMARY FOR", project_id, "\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  cat("✓ Samples:", basic_stats$total_samples, "\n")
  cat("✓ Tumors:", basic_stats$tumor_samples, "\n")
  cat("✓ Normals:", basic_stats$normal_samples, "\n")
  cat("✓ Genes:", basic_stats$total_genes, "\n\n")
  
  if ("final_stage" %in% names(clinical_data)) {
    tumor_data <- clinical_data %>%
      filter(sample_type %in% c("Primary Tumor", "Metastatic"))
    
    stage_completion <- sum(!is.na(tumor_data$final_stage)) / nrow(tumor_data) * 100
    
    cat("STAGE INFORMATION (Enhanced):\n")
    cat("  - Tumors with stage:", sum(!is.na(tumor_data$final_stage)), 
        paste0("(", round(stage_completion, 1), "%)\n"))
    cat("  - Stage sources:\n")
    
    if ("stage_source" %in% names(tumor_data)) {
      sources <- table(tumor_data$stage_source, useNA = "always")
      for (src in names(sources)) {
        if (!is.na(src)) {
          cat("    * ", src, ": ", sources[src], "\n", sep = "")
        }
      }
    }
    
    cat("  - Stage distribution:\n")
    stage_table <- table(tumor_data$final_stage, useNA = "always")
    for (stage in names(stage_table)) {
      if (!is.na(stage) && stage_table[stage] > 0) {
        cat("    * ", stage, ": ", stage_table[stage], "\n", sep = "")
      }
    }
  }
  
  cat("\nQUALITY METRICS (median):\n")
  cat("  - Library size:", 
      round(median(expr_stats$raw_total_counts/1e6), 1), "M reads\n")
  cat("  - Genes expressed:", 
      round(median(expr_stats$raw_genes_expressed/1000), 1), "k genes\n")
  cat("  - Zero-count genes:", 
      round(median(expr_stats$raw_percent_zeros), 1), "%\n")
  
  cat("\nOUTPUT FILES:\n")
  cat("  -", paste0(project_id, "_enhanced_QC.png"), "\n")
  cat("  -", paste0(project_id, "_enhanced_QC.rds"), "\n")
  cat("  - Location:", output_dir, "\n")
  
  return(list(
    enhanced_data = enhanced_data,
    clinical_data = clinical_data,
    expr_stats = expr_stats,
    basic_stats = basic_stats,
    plot = combined_plot
  ))
}

# ==================== HAUPTPROGRAMM ====================
# Für mehrere Projekte:
projects <- c("TCGA-LUAD", "TCGA-BRCA", "TCGA-COAD", "TCGA-BLCA", "TCGA-HNSC", "TCGA-KIRC", "TCGA-LUSC", "TCGA-SKCM", "TCGA-STAD", "TCGA-THCA")
results <- lapply(projects, function(p) {
  run_enhanced_qc(p, data_dir, output_dir)
})

