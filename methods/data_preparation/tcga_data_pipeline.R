# tcga_data_pipeline.R
# Hauptskript für Datendownload und -erweiterung

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(stringr)
library(tidyr)

# ==================== Config ====================
output_base_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/material/TCGA_data"

# get the full clinical data, as some data is not present in the basic download
get_complete_tcga_clinical <- function(project_id, download_dir = NULL) {
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("GETTING COMPLETE CLINICAL DATA FOR:", project_id, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # create dir
  if (!is.null(download_dir)) {
    dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ==================== A. BCR BIOTAB (Main data) ====================
  cat("1. Downloading BCR Biotab (main clinical data)...\n")
  
  biotab_query <- GDCquery(
    project = project_id,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR Biotab"
  )
  
  # Download
  if (!is.null(download_dir)) {
    biotab_dir <- file.path(download_dir, "biotab")
    GDCdownload(biotab_query, directory = biotab_dir)
    capture.output({
      biotab_data <- GDCprepare(biotab_query, directory = biotab_dir)
    })
  } else {
    GDCdownload(biotab_query)
    capture.output({
      biotab_data <- GDCprepare(biotab_query)
    })
  }
  
  patient_table_name <- grep("clinical_patient", names(biotab_data), 
                            value = TRUE, ignore.case = TRUE)[1]
  
  if (is.na(patient_table_name)) {
    stop("Could not find patient table in Biotab data")
  }
  
  cat("   Found patient table:", patient_table_name, "\n")
  patient_data <- biotab_data[[patient_table_name]]
  
  # ==================== B. BCR XML (STAGE INFORMATION) ====================
  cat("\n2. Downloading BCR XML (stage information)...\n")
  
  xml_query <- GDCquery(
    project = project_id,
    data.category = "Clinical",
    data.type = "Clinical Supplement",
    data.format = "BCR XML"
  )
  
  # Download
  if (!is.null(download_dir)) {
    xml_dir <- file.path(download_dir, "xml")
    GDCdownload(xml_query, directory = xml_dir, files.per.chunk = 100)
    capture.output({
      stage_info <- GDCprepare_clinic(xml_query, clinical.info = "stage_event", 
                                   directory = xml_dir)
    })
  } else {
    GDCdownload(xml_query)
    capture.output({
      stage_info <- GDCprepare_clinic(xml_query, clinical.info = "stage_event")
    })
  }
  
  # ==================== C. Other clinical tables from biotab ====================
  cat("\n3. Extracting additional tables from Biotab...\n")
  
  all_tables <- names(biotab_data)
  cat("   Available tables:", paste(all_tables, collapse = ", "), "\n")
  
  important_tables <- list()
  
  for (table_name in all_tables) {
    if (grepl("clinical_(patient|drug|radiation|follow_up|new_tumor|stage_event)", 
              table_name, ignore.case = TRUE)) {
      important_tables[[table_name]] <- biotab_data[[table_name]]
      cat("   -", table_name, ":", nrow(biotab_data[[table_name]]), "rows\n")
    }
  }
  
  # ==================== D. Merge data ====================
  cat("\n4. Merging all clinical data...\n")
  
  # Start mit Patientendaten
  complete_data <- patient_data
  
  cat("   a) Adding stage from XML...\n")
  if (!is.null(stage_info) && nrow(stage_info) > 0) {
    complete_data <- complete_data %>%
      left_join(stage_info %>% 
                 select(bcr_patient_barcode, 
                        pathologic_stage, clinical_stage, tnm_categories),
               by = "bcr_patient_barcode",
               suffix = c("", "_xml"))
  }
  
  # 2. Füge Stage aus Biotab hinzu (falls vorhanden)
  stage_biotab_name <- grep("clinical_stage_event", names(important_tables), 
                           value = TRUE)[1]
  
  if (!is.na(stage_biotab_name)) {
    cat("   b) Adding stage from Biotab...\n")
    stage_biotab <- important_tables[[stage_biotab_name]]
    
    # Wähle Stage-Spalten aus Biotab
    stage_cols_biotab <- grep("stage|Stage", names(stage_biotab), 
                             value = TRUE, ignore.case = TRUE)
    
    if (length(stage_cols_biotab) > 0) {
      stage_biotab_subset <- stage_biotab %>%
        select(bcr_patient_barcode, all_of(stage_cols_biotab))
      
      complete_data <- complete_data %>%
        left_join(stage_biotab_subset,
                 by = "bcr_patient_barcode",
                 suffix = c("", "_biotab"))
    }
  }
  
  # 3. Füge Follow-up hinzu
  followup_name <- grep("clinical_follow_up", names(important_tables), 
                       value = TRUE)[1]
  
  if (!is.na(followup_name)) {
    followup_data <- important_tables[[followup_name]]
    
    # Prüfe ob die Spalte existiert
    if ("days_to_last_followup" %in% names(followup_data)) {
      cat("   c) Adding follow-up data...\n")
      
      latest_followup <- followup_data %>%
        group_by(bcr_patient_barcode) %>%
        arrange(desc(days_to_last_followup)) %>%
        slice(1) %>%
        ungroup() %>%
        select(bcr_patient_barcode, 
               days_to_last_followup, 
               vital_status_followup = vital_status)
      
      complete_data <- complete_data %>%
        left_join(latest_followup, by = "bcr_patient_barcode")
    } else {
      cat("   c) Skipping follow-up (column 'days_to_last_followup' not found)\n")
    }
  }
  
  # ==================== E. STAGE BEREINIGUNG & KONSOLIDIERUNG ====================
  cat("\n5. Consolidating stage information...\n")
  
  complete_data <- complete_data %>%
    mutate(
      # PRIORITÄT: 1. XML pathologic_stage, 2. Biotab stage, 3. patient table stage
      final_pathologic_stage = case_when(
        !is.na(pathologic_stage) & pathologic_stage != "" ~ pathologic_stage,
        !is.na(ajcc_pathologic_tumor_stage) & ajcc_pathologic_tumor_stage != "" ~ 
          paste("Stage", ajcc_pathologic_tumor_stage),
        !is.na(clinical_stage) & clinical_stage != "" ~ clinical_stage,
        TRUE ~ NA_character_
      ),
      
      # Vereinfachte Stage (I, II, III, IV)
      simple_stage = case_when(
        # Stage I: Must be "Stage I" optionally followed by A, B, C but NOT II, III, IV
        grepl("^Stage I([ABC]?)$", final_pathologic_stage) ~ "Stage I",
        
        # Stage II: Must be "Stage II" optionally followed by A, B, C
        grepl("^Stage II([ABC]?)$", final_pathologic_stage) ~ "Stage II",
        
        # Stage III: Must be "Stage III" optionally followed by A, B, C
        grepl("^Stage III([ABC]?)$", final_pathologic_stage) ~ "Stage III",
        
        # Stage IV: Must be "Stage IV" optionally followed by A, B, C
        grepl("^Stage IV([ABC]?)$", final_pathologic_stage) ~ "Stage IV",
        
        TRUE ~ NA_character_
      ),
      
      # Markiere die Stage-Quelle
      stage_source = case_when(
        !is.na(pathologic_stage) & pathologic_stage != "" ~ "XML",
        !is.na(ajcc_pathologic_tumor_stage) & ajcc_pathologic_tumor_stage != "" ~ "Biotab",
        !is.na(clinical_stage) & clinical_stage != "" ~ "XML (clinical)",
        TRUE ~ "Missing"
      )
    )
  
  # ==================== F. QUALITÄTSBERICHT ====================
  cat("\n6. Quality report:\n")
  cat("   Total patients:", nrow(complete_data), "\n")
  
  # Stage statistics
  stage_count <- sum(!is.na(complete_data$final_pathologic_stage))
  cat("   Patients with stage info:", stage_count, 
      paste0("(", round(stage_count/nrow(complete_data)*100, 1), "%)\n"))
  
  cat("   Stage sources:\n")
  stage_source_table <- table(complete_data$stage_source)
  for (source in names(stage_source_table)) {
    cat("     -", source, ":", stage_source_table[source], "\n")
  }
  
  cat("\n   Final stage distribution:\n")
  final_stage_table <- table(complete_data$final_pathologic_stage, useNA = "always")
  for (stage in names(final_stage_table)) {
    if (!is.na(stage) && final_stage_table[stage] > 0) {
      cat("     -", stage, ":", final_stage_table[stage], "\n")
    }
  }
  
  # ==================== G. RETURN ====================
  return(list(
    complete_clinical = complete_data,
    biotab_tables = important_tables,
    xml_stage = stage_info,
    patient_table = patient_data,
    summary = list(
      total_patients = nrow(complete_data),
      with_stage = stage_count,
      stage_completeness = round(stage_count/nrow(complete_data)*100, 1),
      stage_sources = stage_source_table,
      stage_distribution = final_stage_table
    )
  ))
}

# ==================== Get expression data and add missing clinical information ====================

download_and_enhance_tcga_data <- function(project_id, output_dir) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("TCGA DATA PIPELINE FOR:", project_id, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  # Schritt 1: Expression Data herunterladen
  cat("STEP 1: DOWNLOADING EXPRESSION DATA\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  expr_query <- GDCquery(
    project = project_id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download mit tryCatch für den Fall, dass Daten bereits existieren
  tryCatch({
    GDCdownload(expr_query, directory = output_dir)
    cat("✓ Downloaded expression data\n")
  }, error = function(e) {
    cat("⚠ Data already exists or download error:", e$message, "\n")
    cat("  Attempting to continue with existing data...\n")
  })
  capture.output({
    expr_data <- GDCprepare(expr_query, directory = output_dir)
  })
  
  cat("✓ Expression data prepared\n")
  cat("  - Samples:", ncol(expr_data), "\n")
  cat("  - Genes:", nrow(expr_data), "\n")

  original_colnames <- colnames(expr_data)
  original_rownames <- rownames(expr_data)
  
  # Extrahiere Aliquot IDs
  expr_patient_ids <- unique(substr(colnames(expr_data), 1, 12))
  cat("  - Aliquots:", length(expr_patient_ids), "\n")
  
  # Schritt 2: Klinische Daten extrahieren und Stage-Informationen finden
  cat("\nSTEP 2: EXTRACTING CLINICAL DATA\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  clinical_data <- as.data.frame(colData(expr_data))
  
  # Finde Stage-Spalte in den klinischen Daten
  stage_cols <- grep("ajcc.*stage|pathologic.*stage", names(clinical_data), 
                     value = TRUE, ignore.case = TRUE)
  
  if (length(stage_cols) > 0) {
    primary_stage_col <- stage_cols[1]
    cat("  Found stage column:", primary_stage_col, "\n")
    
    # Analysiere Stage-Verteilung
    stage_dist <- table(clinical_data[[primary_stage_col]], useNA = "always")
    cat("  Initial stage distribution:\n")
    print(stage_dist)
  } else {
    cat("  ⚠ No stage column found in clinical data\n")
    # Setze einen Standard-Namen
    primary_stage_col <- "ajcc_pathologic_stage"
  }
  
  # Schritt 3: Erweiterte klinische Daten holen
  cat("\nSTEP 3: ENHANCING CLINICAL DATA\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  enhanced_clinical <- get_complete_tcga_clinical(project_id, output_dir)
  
  # Schritt 4: Intelligente Zusammenführung (ergänze nur NAs)
  cat("\nSTEP 4: INTELLIGENT DATA INTEGRATION\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Erweiterte klinische Daten filtern
  enhanced_data <- enhanced_clinical$complete_clinical %>%
    filter(bcr_patient_barcode %in% expr_patient_ids)
  
  cat("  Patients in expression data:", length(expr_patient_ids), "\n")
  cat("  Patients in enhanced clinical:", nrow(enhanced_data), "\n")
  
  # Erstelle Mapping für Stage-Ergänzung
  stage_mapping <- data.frame(
    bcr_patient_barcode = expr_patient_ids,
    stringsAsFactors = FALSE
  )
  
  # Join mit original stage data (wenn verfügbar)
  if (primary_stage_col %in% names(clinical_data)) {
    stage_mapping <- stage_mapping %>%
      left_join(
        clinical_data %>%
          mutate(bcr_patient_barcode = substr(barcode, 1, 12)) %>%
          select(bcr_patient_barcode, stage_original = !!sym(primary_stage_col)) %>%
          distinct(bcr_patient_barcode, .keep_all = TRUE),
        by = "bcr_patient_barcode"
      )
  } else {
    stage_mapping$stage_original <- NA_character_
  }
  
  # Join mit enhanced stage data
  stage_mapping <- stage_mapping %>%
    left_join(
      enhanced_data %>%
        select(bcr_patient_barcode, stage_enhanced = final_pathologic_stage),
      by = "bcr_patient_barcode"
    ) %>%
    mutate(
      # join new information when it is missing in the current table
      final_stage = case_when(
        !is.na(stage_original) & 
          stage_original != "" & 
          stage_original != "[Not Available]" &
          stage_original != "[Discrepancy]" ~ stage_original,
        
        !is.na(stage_enhanced) & 
          stage_enhanced != "" & 
          stage_enhanced != "[Not Available]" &
          stage_enhanced != "[Discrepancy]" ~ stage_enhanced,
        
        TRUE ~ NA_character_
      ),
      
      stage_source = case_when(
        !is.na(stage_original) & 
          stage_original != "" & 
          stage_original != "[Not Available]" &
          stage_original != "[Discrepancy]" ~ "Expression_Data",
        
        !is.na(stage_enhanced) & 
          stage_enhanced != "" & 
          stage_enhanced != "[Not Available]" &
          stage_enhanced != "[Discrepancy]" ~ "Enhanced_Clinical",
        
        TRUE ~ "Missing"
      ),
      
      simple_stage = case_when(
        # Stage I: Must be "Stage I" optionally followed by A, B, C but NOT II, III, IV
        grepl("^Stage I([ABC]?)$", final_stage) ~ "Stage I",
        
        # Stage II: Must be "Stage II" optionally followed by A, B, C
        grepl("^Stage II([ABC]?)$", final_stage) ~ "Stage II",
        
        # Stage III: Must be "Stage III" optionally followed by A, B, C
        grepl("^Stage III([ABC]?)$", final_stage) ~ "Stage III",
        
        # Stage IV: Must be "Stage IV" optionally followed by A, B, C
        grepl("^Stage IV([ABC]?)$", final_stage) ~ "Stage IV",
        
        TRUE ~ NA_character_
      ),
      
      # Flag für ergänzte Patienten
      was_supplemented = ifelse(
        (is.na(stage_original) | 
           stage_original == "" | 
           stage_original == "[Not Available]" |
           stage_original == "[Discrepancy]") &
          (!is.na(stage_enhanced) & 
             stage_enhanced != "" & 
             stage_enhanced != "[Not Available]" &
             stage_enhanced != "[Discrepancy]"),
        "Yes", "No"
      )
    )
  
  # Bericht über Ergänzung
  cat("\n  Stage supplementation results:\n")
  original_with_stage <- sum(!is.na(stage_mapping$stage_original) & 
                              stage_mapping$stage_original != "[Not Available]" &
                              stage_mapping$stage_original != "[Discrepancy]", na.rm = TRUE)
  
  final_with_stage <- sum(!is.na(stage_mapping$final_stage), na.rm = TRUE)
  supplemented <- sum(stage_mapping$was_supplemented == "Yes", na.rm = TRUE)
  
  cat("    Original: ", original_with_stage, 
      " patients with stage (", 
      round(original_with_stage/nrow(stage_mapping)*100, 1), "%)\n", sep = "")
  
  cat("    Final:    ", final_with_stage, 
      " patients with stage (", 
      round(final_with_stage/nrow(stage_mapping)*100, 1), "%)\n", sep = "")
  
  cat("    Supplemented: ", supplemented, 
      " additional patients (+", 
      round(supplemented/nrow(stage_mapping)*100, 1), "%)\n", sep = "")
  
  # Detaillierte Stage-Verteilung - KORRIGIERTE VERSION
  cat("\n  Final stage distribution:\n")
  
  # Erstelle eine sichere Version der Stage-Verteilung
  final_dist_df <- stage_mapping %>%
    mutate(final_stage_display = ifelse(is.na(final_stage), "Missing", final_stage)) %>%
    count(final_stage_display) %>%
    arrange(desc(n)) %>%
    mutate(percentage = round(n / sum(n) * 100, 1))
  
  for (i in 1:nrow(final_dist_df)) {
    cat("    - ", final_dist_df$final_stage_display[i], ": ", 
        final_dist_df$n[i], " (", final_dist_df$percentage[i], "%)\n", sep = "")
  }
  
  # Schritt 5: Erweiterte Daten zum SummarizedExperiment hinzufügen
  cat("\nSTEP 5: INTEGRATING ENHANCED DATA\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Erstelle erweitertes colData
  enhanced_col_data <- clinical_data %>%
    mutate(bcr_patient_barcode = substr(barcode, 1, 12)) %>%
    left_join(
      stage_mapping %>%
        select(bcr_patient_barcode, 
               final_stage, 
               simple_stage,
               stage_source,
               stage_original,
               stage_enhanced,
               was_supplemented),
      by = "bcr_patient_barcode"
    )
  
  # Ersetze das colData
  colData(expr_data) <- DataFrame(enhanced_col_data)

  colnames(expr_data) <- original_colnames
  rownames(expr_data) <- original_rownames
  
  # Schritt 6: Speichern der Ergebnisse
  cat("\nSTEP 6: SAVING RESULTS\n")
  cat(paste(rep("-", 60), collapse = ""), "\n")
  
  # Erstelle Projekt-spezifischen Ordner
  project_dir <- file.path(output_dir, project_id)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Dateinamen
  rds_file <- file.path(project_dir, paste0(project_id, "_enhanced.rds"))
  clinical_file <- file.path(project_dir, paste0(project_id, "_clinical_enhanced.csv"))
  summary_file <- file.path(project_dir, paste0(project_id, "_summary.txt"))
  mapping_file <- file.path(project_dir, paste0(project_id, "_stage_mapping.csv"))
  
  # Speichere erweitertes RDS
  saveRDS(expr_data, rds_file)
  cat("  ✓ Saved enhanced data to:", rds_file, "\n")
  
  cat("  Extracting essential clinical columns...\n")

  # Definiere die wichtigsten Spalten
  essential_columns <- c(
    "barcode", "patient", "sample", "tissue_type", "sample_type",
    "bcr_patient_barcode", "final_stage", "simple_stage", "stage_source",
    "stage_original", "stage_enhanced", "was_supplemented",
    "gender", "age_at_diagnosis", "vital_status", "days_to_death",
    "days_to_last_followup", "ajcc_pathologic_stage"
  )

  # Finde verfügbare Spalten
  available_columns <- intersect(essential_columns, names(enhanced_col_data))

  # Extrahiere nur einfache Spalten (keine Listen)
  simple_cols <- sapply(enhanced_col_data[, available_columns], 
                        function(x) !is.list(x) && !is(x, "DataFrame"))

  clinical_data_simple <- enhanced_col_data[, available_columns[simple_cols]]

  # Speichere als CSV
  write.csv(clinical_data_simple, clinical_file, row.names = FALSE)
  cat("  ✓ Saved clinical data to:", clinical_file, "\n")
  
  # Speichere Stage-Mapping als CSV
  write.csv(stage_mapping, mapping_file, row.names = FALSE)
  cat("  ✓ Saved stage mapping to:", mapping_file, "\n")
  
  # Erstelle Zusammenfassungsdatei
  sink(summary_file)
  
  cat("TCGA DATA SUMMARY -", project_id, "\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  cat("EXPRESSION DATA:\n")
  cat("  Samples:", ncol(expr_data), "\n")
  cat("  Genes:", nrow(expr_data), "\n")
  cat("  Patients:", length(expr_patient_ids), "\n\n")
  
  cat("STAGE INFORMATION:\n")
  cat("  Original stage column:", primary_stage_col, "\n")
  cat("  Patients with original stage:", original_with_stage, "\n")
  cat("  Patients with enhanced stage:", final_with_stage, "\n")
  cat("  Stage completion rate:", 
      round(final_with_stage / length(expr_patient_ids) * 100, 1), "%\n")
  cat("  Patients supplemented:", supplemented, "\n\n")
  
  cat("STAGE DISTRIBUTION (Enhanced):\n")
  for (i in 1:nrow(final_dist_df)) {
    cat("  - ", final_dist_df$final_stage_display[i], ": ", 
        final_dist_df$n[i], " (", final_dist_df$percentage[i], "%)\n", sep = "")
  }
  
  cat("\nSTAGE SOURCES:\n")
  source_table <- table(stage_mapping$stage_source, useNA = "always")
  for (source in names(source_table)) {
    if (!is.na(source) && source_table[source] > 0) {
      cat("  - ", source, ": ", source_table[source], "\n", sep = "")
    }
  }
  
  cat("\nENHANCED CLINICAL SUMMARY:\n")
  cat("  Total patients in enhanced clinical:", nrow(enhanced_data), "\n")
  cat("  Stage completeness in enhanced:", 
      enhanced_clinical$summary$stage_completeness, "%\n")
  
  sink()
  
  cat("  ✓ Saved summary to:", summary_file, "\n")
  
  # Schritt 7: Return
  return(list(
    enhanced_se = expr_data,
    stage_mapping = stage_mapping,
    enhanced_clinical = enhanced_clinical,
    enhanced_col_data = enhanced_col_data,
    summary = list(
      project = project_id,
      samples = ncol(expr_data),
      patients = length(expr_patient_ids),
      original_with_stage = original_with_stage,
      final_with_stage = final_with_stage,
      supplemented = supplemented,
      stage_completion = round(final_with_stage / length(expr_patient_ids) * 100, 1),
      stage_distribution = final_dist_df,
      files = list(
        rds = rds_file,
        clinical = clinical_file,
        mapping = mapping_file,
        summary = summary_file
      )
    )
  ))
}

analyze_multiple_projects <- function(project_list, output_dir) {
  results <- list()
  
  for (project in project_list) {
    cat("\n\n", paste(rep("#", 100), collapse = ""))
    cat("\nSTARTING PIPELINE FOR:", project, "\n")
    cat(paste(rep("#", 100), collapse = ""), "\n\n")
    
    result <- tryCatch({
      download_and_enhance_tcga_data(project, output_dir)
    }, error = function(e) {
      cat("Error processing", project, ":", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(result)) {
      results[[project]] <- result
    }
  }
  
  return(results)
}

tcga_projects <- c("TCGA-LUAD", "TCGA-BRCA", "TCGA-COAD", "TCGA-BLCA", "TCGA-HNSC", "TCGA-KIRC", "TCGA-LUSC", "TCGA-SKCM", "TCGA-STAD", "TCGA-THCA")

all_results <- analyze_multiple_projects(tcga_projects, output_base_dir)

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("PIPELINE COMPLETED!\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

if (!is.null(result$error)) {
  cat("⚠ Pipeline completed with warnings:\n")
  cat("  ", result$error, "\n")
  cat("  Loaded existing data instead.\n\n")
} else {
  cat("✓ Pipeline completed successfully!\n\n")
}

cat("Enhanced data saved to:\n")
if (!is.null(result$summary)) {
  cat("  RDS file:", result$summary$files$rds, "\n")
  cat("  Clinical CSV:", result$summary$files$clinical, "\n")
  cat("  Stage mapping:", result$summary$files$mapping, "\n")
  cat("  Summary:", result$summary$files$summary, "\n\n")
  
  if (result$summary$supplemented > 0) {
    cat("✓ Stage completion improved from", 
        result$summary$original_with_stage, "to", 
        result$summary$final_with_stage, "patients (",
        result$summary$supplemented, "added)\n")
  } else {
    cat("✓ Stage information complete: ", 
        result$summary$final_with_stage, "/", 
        result$summary$patients, " patients\n", sep = "")
  }
} else {
  cat("  RDS file:", file.path(output_base_dir, project_id, 
                              paste0(project_id, "_enhanced.rds")), "\n")
  cat("  (Using existing data)\n")
}

