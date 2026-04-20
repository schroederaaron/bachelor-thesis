#!/usr/bin/env Rscript
# Führt differentielle Expressionsanalyse mit edgeR und limma/voom durch
# für jedes Krebsprojekt (TCGA) und jedes Stadium (I-IV).
# Vergleich: Cancer vs. Healthy (constant reference)
# Verwendet prepare_stage_data.R Ausgabeformat

# Installiere und lade benötigte Pakete
required_packages <- c("BiocManager", "edgeR", "limma", "dplyr", "ggplot2", "VennDiagram", 
                       "rtracklayer", "GenomicFeatures")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("edgeR", "limma", "BiocManager", "rtracklayer", "GenomicFeatures")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

source("config.R")        # enthält BASE_DATA_DIR, FAMILY_OUTPUT_DIR, STAGES, PROJECTS, etc.
source("utils.R")         # für get_norm_suffix, get_norm_output_dir, etc.

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

# Ausgabeverzeichnis für differentielle Expression
DE_OUTPUT_DIR <- file.path(dirname(GENE_OUTPUT_DIR), "differential_expression")
dir.create(DE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Signifikanzschwellen
FDR_THRESHOLD <- 0.05
LOGFC_THRESHOLD <- 2

# ==================== DATEN LADEN ====================

#' Lade Cancer-Daten für ein Projekt und eine Stage (Counts)
#' @param project_id Projekt-ID (z.B. "TCGA-BRCA")
#' @param stage Stage (z.B. "Stage I")
#' @return Counts-Matrix (genes × samples)
load_cancer_counts <- function(project_id, stage) {
  
  filename <- paste0(project_id, "-", gsub(" ", "-", stage), "_counts.rds")
  file_path <- file.path(BASE_DATA_DIR, project_id, filename)
  
  if (!file.exists(file_path)) {
    warning(sprintf("Datei nicht gefunden: %s", file_path))
    return(NULL)
  }
  
  data_obj <- readRDS(file_path)
  
  # Extrahiere Counts (samples × genes -> genes × samples)
  if (is.matrix(data_obj$expression_vectors)) {
    counts <- t(data_obj$expression_vectors)
    colnames(counts) <- rownames(data_obj$expression_vectors)
    rownames(counts) <- colnames(data_obj$expression_vectors)
  } else {
    warning("expression_vectors ist keine Matrix")
    return(NULL)
  }
  
  # Extrahiere Gen-IDs (Uniprot)
  gene_ids <- data_obj$gene_ids
  
  # Stelle sicher, dass die Reihenfolge stimmt
  if (!identical(rownames(counts), gene_ids)) {
    warning("Gen-Reihenfolge stimmt nicht überein!")
    return(NULL)
  }
  
  return(counts)
}

#' Lade Healthy-Daten für ein Projekt (Counts)
#' @param project_id Projekt-ID (z.B. "TCGA-BRCA")
#' @return Counts-Matrix (genes × samples)
load_healthy_counts <- function(project_id) {
  
  filename <- paste0("healthy_", project_id, "_counts.rds")
  file_path <- file.path(BASE_DATA_DIR, project_id, filename)
  
  if (!file.exists(file_path)) {
    warning(sprintf("Healthy-Datei nicht gefunden: %s", file_path))
    return(NULL)
  }
  
  data_obj <- readRDS(file_path)
  
  # Extrahiere Counts (samples × genes -> genes × samples)
  if (is.matrix(data_obj$expression_vectors)) {
    counts <- t(data_obj$expression_vectors)
    colnames(counts) <- rownames(data_obj$expression_vectors)
    rownames(counts) <- colnames(data_obj$expression_vectors)
  } else {
    warning("expression_vectors ist keine Matrix")
    return(NULL)
  }
  
  return(counts)
}

# ==================== DIFFERENTIELLE EXPRESSIONSANALYSE ====================

#' Führe edgeR-Analyse durch (Quasi-Likelihood F-Test)
#' @param counts Counts-Matrix (genes × samples)
#' @param group Faktor mit Gruppenzugehörigkeit
#' @return Data.Frame mit Ergebnissen
run_edgeR <- function(counts, group) {
  
  # DGEList erstellen
  y <- DGEList(counts = counts, group = group)
  
  # Filtere niedrig exprimierte Gene
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  # TMM-Normalisierung
  y <- calcNormFactors(y)
  
  # Design-Matrix (nur zwei Gruppen: healthy vs cancer)
  design <- model.matrix(~ group)
  
  # Dispersion schätzen
  y <- estimateDisp(y, design)
  
  # Quasi-Likelihood F-Test
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = 2)  # cancer vs healthy
  
  # Ergebnisse extrahieren
  result <- topTags(qlf, n = Inf)
  result_df <- result$table
  result_df$gene_id <- rownames(result_df)
  
  return(result_df)
}

#' Führe limma-voom-Analyse durch
#' @param counts Counts-Matrix (genes × samples)
#' @param group Faktor mit Gruppenzugehörigkeit
#' @return Data.Frame mit Ergebnissen
run_limma_voom <- function(counts, group) {
  
  # DGEList erstellen
  dge <- DGEList(counts = counts, group = group)
  
  # Filtere niedrig exprimierte Gene
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # TMM-Normalisierung
  dge <- calcNormFactors(dge)
  
  # Design-Matrix
  design <- model.matrix(~ group)
  
  # voom-Transformation (ohne Genlängen-Offsets)
  v <- voom(dge, design, plot = FALSE)
  
  # Lineares Modell
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # Ergebnisse für den Kontrast cancer - healthy
  result <- topTable(fit, coef = 2, number = Inf)
  result$gene_id <- rownames(result)
  
  return(result)
}
# ==================== VISUALISIERUNGEN ====================

#' Erstelle Volcano Plot
plot_volcano <- function(results, pval_col = "PValue", adj_pval_col = "FDR", 
                         title = "", filename = NULL) {
  
  results$Significant <- ifelse(results[[adj_pval_col]] < FDR_THRESHOLD, "Yes", "No")
  
  p <- ggplot(results, aes(x = logFC, y = -log10(.data[[pval_col]]), color = Significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("grey", "red")) +
    labs(title = title, x = "log2 Fold Change", y = "-log10(p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  if (!is.null(filename)) {
    ggsave(filename, p, width = 8, height = 6, dpi = 300)
  }
  return(p)
}

#' Erstelle Venn-Diagramm für verschiedene Filterkriterien
create_venn <- function(results, prefix = "", filename = NULL) {
  
  # Filter nach logFC >= 2
  genes_logFC <- results$gene_id[abs(results$logFC) >= LOGFC_THRESHOLD]
  # Filter nach FDR <= 0.05
  genes_fdr05 <- results$gene_id[results$FDR <= 0.05]
  # Filter nach FDR <= 0.01
  genes_fdr01 <- results$gene_id[results$FDR <= 0.01]
  
  if (length(genes_logFC) == 0 || length(genes_fdr05) == 0 || length(genes_fdr01) == 0) {
    warning("Nicht genug Gene für Venn-Diagramm")
    return(NULL)
  }
  
  venn.diagram(
    x = list(
      "|logFC| >= 2" = genes_logFC,
      "FDR <= 0.05" = genes_fdr05,
      "FDR <= 0.01" = genes_fdr01
    ),
    filename = filename,
    col = "black",
    fill = c("red", "blue", "green"),
    alpha = 0.5,
    label.col = "white",
    cex = 1.5,
    cat.cex = 1.5,
    cat.pos = 0,
    main = prefix
  )
}

#' Extrahiere Top 20 up- und down-regulierte Gene
extract_top_genes <- function(results, fdr_threshold = 0.05, logfc_threshold = 2) {
  
  sig <- results[results$FDR < fdr_threshold, ]
  
  if (nrow(sig) == 0) {
    return(list(up = data.frame(), down = data.frame()))
  }
  
  up <- sig[sig$logFC >= logfc_threshold, ]
  down <- sig[sig$logFC <= -logfc_threshold, ]
  
  # Sortieren

# Für Plot 3: Anz
  up <- up[order(up$logFC, decreasing = TRUE), ]
  down <- down[order(down$logFC, decreasing = FALSE), ]
  
  top_up <- head(up, 20)
  top_down <- head(down, 20)
  
  return(list(up = top_up, down = top_down))
}

# Für die Top-Gene und Plots:
process_de_results <- function(res, method_name, out_dir) {
  if (is.null(res) || nrow(res) == 0) return()
  
  # Speichern
  write.csv(res, file.path(out_dir, paste0(method_name, "_results.csv")), row.names = FALSE)
  
  # Volcano Plot
  if (method_name == "edgeR") {
    plot_volcano(res, pval_col = "PValue", adj_pval_col = "FDR",
                 title = paste(basename(out_dir), "- edgeR"),
                 filename = file.path(out_dir, "volcano_edgeR.png"))
  } else {
    res$FDR <- res$adj.P.Val  # Für einheitliche Spaltennamen
    plot_volcano(res, pval_col = "P.Value", adj_pval_col = "adj.P.Val",
                 title = paste(basename(out_dir), "- voom"),
                 filename = file.path(out_dir, "volcano_voom.png"))
  }
  
  # Venn-Diagramm
  create_venn(res, prefix = paste(basename(out_dir), method_name),
              filename = file.path(out_dir, paste0("venn_", method_name, ".png")))
  
  # Top-Gene
  top_genes <- extract_top_genes(res)
  if (nrow(top_genes$up) > 0)
    write.csv(top_genes$up, file.path(out_dir, paste0("top20_up_", method_name, ".csv")), 
              row.names = FALSE)
  if (nrow(top_genes$down) > 0)
    write.csv(top_genes$down, file.path(out_dir, paste0("top20_down_", method_name, ".csv")), 
              row.names = FALSE)
}


# ==================== PROJEKTE VERARBEITEN ====================

for (cancer_name in names(CANCER_TYPES)) {
  project_id <- CANCER_TYPES[cancer_name]
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Verarbeite:", cancer_name, "(", project_id, ")\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Lade Healthy-Daten (konstant)
  healthy_counts <- load_healthy_counts(project_id)
  if (is.null(healthy_counts)) {
    cat("  Keine Healthy-Daten. Überspringe.\n")
    next
  }
  
  for (stage in STAGES) {
    cat("\n  --- Stage:", stage, "---\n")
    
    cancer_counts <- load_cancer_counts(project_id, stage)
    if (is.null(cancer_counts)) {
      cat("    Keine Cancer-Daten. Überspringe.\n")
      next
    }
    
    # Gemeinsame Gene
    common_genes <- intersect(rownames(healthy_counts), rownames(cancer_counts))
    if (length(common_genes) == 0) {
      cat("    Keine gemeinsamen Gene. Überspringe.\n")
      next
    }
    
    healthy_sub <- healthy_counts[common_genes, , drop = FALSE]
    cancer_sub <- cancer_counts[common_genes, , drop = FALSE]
    
    cat(sprintf("    Healthy: %d Samples, Cancer: %d Samples\n", 
                ncol(healthy_sub), ncol(cancer_sub)))
    
    # Kombinierte Counts
    counts <- cbind(healthy_sub, cancer_sub)
    group <- factor(c(rep("healthy", ncol(healthy_sub)),
                      rep("cancer", ncol(cancer_sub))))
    group <- relevel(group, ref = "healthy")
    
    # Ausgabeverzeichnis
    out_dir <- file.path(DE_OUTPUT_DIR, project_id, stage)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # edgeR
    cat("    Führe edgeR aus...\n")
    edgeR_res <- tryCatch(run_edgeR(counts, group), error = function(e) NULL)
    process_de_results(edgeR_res, "edgeR", out_dir)
    
    # limma/voom
    cat("    Führe limma/voom aus...\n")
    voom_res <- tryCatch(run_limma_voom(counts, group), error = function(e) NULL)
    process_de_results(voom_res, "voom", out_dir)
  }
}

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("DIFFERENTIELLE EXPRESSIONSANALYSE ABGESCHLOSSEN\n")
cat("Ergebnisse in:", DE_OUTPUT_DIR, "\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
