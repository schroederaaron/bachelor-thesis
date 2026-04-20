# prepare_stage_data.R
# Script to prepare stage-specific data in TensorOmics format
# Creates RDS files for each cancer project and stage in results directory
# Also creates separate RDS files for healthy/normal tissue data
# Special handling for TCGA-HNSC: splits into two subprojects based on Primary site

library(SummarizedExperiment)
library(dplyr)

# ==================== CONFIGURATION ====================

# Define the cancer projects to process
projects <- c(
  "TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD",
  "TCGA-KIRC", "TCGA-LUAD", "TCGA-LUSC",
  "TCGA-STAD", "TCGA-THCA"
)

# Stages to process
stages <- c("Stage I", "Stage II", "Stage III", "Stage IV")

# Input directory containing filtered RDS files
input_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/results/TCGA_QC_Reports/"

# Output directory for prepared stage data
output_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/results/"

# Column name for tissue type
tissue_type_col <- "sample_type"

# ==================== HELPER FUNCTIONS ====================

#' Create family mapping arrays
#' @param families Dataframe with gene-family assignments
#' @param gene_ids Vector of all gene IDs in expression data
#' @return List with family mapping arrays
create_family_arrays <- function(families, gene_ids) {
  # Extract unique families
  family_ids <- unique(families$Familie)
  n_families <- length(family_ids)
  
  # Create mapping from family name to index
  family_index_map <- setNames(1:n_families, family_ids)
  
  # Initialize gene_to_fam with zeros (0 = no family assignment)
  gene_to_fam <- integer(length(gene_ids))
  names(gene_to_fam) <- gene_ids
  
  # Fill gene_to_fam array
  for (i in seq_along(gene_ids)) {
    gene <- gene_ids[i]
    family_row <- families[families$Protein == gene, ]
    if (nrow(family_row) > 0) {
      family_name <- family_row$Familie[1]
      gene_to_fam[i] <- family_index_map[[family_name]]
    }
  }
  
  return(list(
    family_ids = family_ids,
    family_index_map = family_index_map,
    gene_to_fam = gene_to_fam,
    n_families = n_families
  ))
}

#' Filter expression data for healthy/normal tissue samples
#' @param expr_matrix Expression matrix (genes x samples)
#' @param metadata Sample metadata
#' @param tissue_type_col Column name for tissue type
#' @param healthy_value Value indicating healthy/normal tissue
#' @return Filtered expression matrix with only healthy samples
filter_healthy_samples <- function(expr_matrix, metadata, 
                                  tissue_type_col = "sample_type",
                                  healthy_value = "Solid Tissue Normal") {
  
  if (is.null(metadata) || !(tissue_type_col %in% colnames(metadata))) {
    cat("  Warning: No tissue type metadata found, cannot filter healthy samples\n")
    return(NULL)
  }
  
  # Get sample tissue types
  sample_tissue_types <- metadata[[tissue_type_col]]
  names(sample_tissue_types) <- rownames(metadata)
  
  # Identify healthy samples
  healthy_samples <- names(sample_tissue_types)[
    sample_tissue_types == healthy_value & 
    !is.na(sample_tissue_types)
  ]
  
  if (length(healthy_samples) == 0) {
    cat(sprintf("  No healthy samples found (looking for '%s')\n", healthy_value))
    return(NULL)
  }
  
  # Filter expression matrix to healthy samples only
  healthy_expr <- expr_matrix[, colnames(expr_matrix) %in% healthy_samples, drop = FALSE]
  
  cat(sprintf("  Found %d healthy samples\n", ncol(healthy_expr)))
  
  return(healthy_expr)
}

#' Prepare data for a specific cancer project (stage or healthy)
#' @param project_id Cancer project ID (e.g., "TCGA-LUAD")
#' @param stage_name Stage name (e.g., "Stage I") or "healthy" for normal tissue
#' @param input_dir Input directory with filtered RDS files
#' @param stage_col Column name for stage information
#' @param tissue_type_col Column name for tissue type
#' @param data_type Type of data: "stage" or "healthy"
#' @param primary_site_filter Optional list with primary site filtering parameters
#' @return List with prepared data or NULL if no data
prepare_project_data <- function(project_id, stage_name = NULL, 
                                input_dir = input_dir, 
                                stage_col = "simple_stage",
                                tissue_type_col = "sample_type",
                                data_type = "stage",
                                primary_site_filter = NULL,
                                tpms = TRUE) {
  
  if (data_type == "stage") {
    cat(sprintf("Preparing %s - %s...\n", project_id, stage_name))
  } else if (data_type == "healthy") {
    cat(sprintf("Preparing %s - Healthy tissue...\n", project_id))
  }
  
  # Load filtered data with stage information
  if(tpms){
    data_file <- file.path(input_dir, paste0(project_id, "_filtered.rds"))
  } else{
    data_file <- file.path(input_dir, paste0(project_id, "_filtered_counts.rds"))
  }
  
  
  if (!file.exists(data_file)) {
    cat(sprintf("  Error: Data file not found: %s\n", data_file))
    return(NULL)
  }
  
  filtered_data <- readRDS(data_file)
  
  # Extract expression matrix (Uniprot IDs as rownames)
  expr_matrix <- filtered_data$dataframe
  # Extract metadata from the new format
  if ("metadata" %in% names(filtered_data) && 
      "sample_metadata" %in% names(filtered_data$metadata)) {
    # New format with embedded metadata
    metadata <- filtered_data$metadata$sample_metadata
    cat("  Using embedded sample metadata\n")
  } else {
    # Old format or no metadata
    metadata <- NULL
    cat("  Warning: No embedded metadata found\n")
  }
  
  cat("Filtering stage data\n") 
  # Filter samples based on data type
  if (data_type == "stage") {
    # Filter for specific cancer stage
    if (!is.null(metadata) && stage_col %in% colnames(metadata)) {
      
      # Get sample stages and tissue types from metadata
      sample_stages <- metadata[[stage_col]]
      names(sample_stages) <- rownames(metadata)
      
      sample_tissue_types <- metadata[[tissue_type_col]]
      names(sample_tissue_types) <- rownames(metadata)
      
      tumor_samples <- names(sample_tissue_types)[
        sample_tissue_types == "Primary Tumor" &
        !is.na(sample_tissue_types)
      ]
      
      # avoid healthy samples in stage data
      stage_samples <- intersect(
        names(sample_stages)[sample_stages == stage_name & !is.na(sample_stages)],
        tumor_samples
      )
      
      if (length(stage_samples) == 0) {
        cat(sprintf("  No tumor samples found for %s - %s\n", project_id, stage_name))
        return(NULL)
      }
      
      # Filter expression matrix
      expr_filtered <- expr_matrix[, colnames(expr_matrix) %in% stage_samples, drop = FALSE]
      sample_names <- stage_samples
      
    } else {
      cat("  Warning: No stage metadata, cannot filter reliably\n")
      return(NULL)
    }
    
  } else if (data_type == "healthy") {
    # Filter for healthy tissue samples (bleibt wie gehabt)
    expr_filtered <- filter_healthy_samples(expr_matrix, metadata, tissue_type_col)
    
    if (is.null(expr_filtered) || ncol(expr_filtered) == 0) {
      cat(sprintf("  No healthy tissue data for %s\n", project_id))
      return(NULL)
    }
    
    sample_names <- colnames(expr_filtered)
    stage_name <- "healthy"
  }
  
  if (ncol(expr_filtered) == 0) {
    cat(sprintf("  No expression data for %s - %s\n", project_id, stage_name))
    return(NULL)
  }
  
  # transpose to get tox format
  expr_tensoromics <- t(expr_filtered)
  
  cat(sprintf("  Found %d samples, %d genes\n", 
              nrow(expr_tensoromics), ncol(expr_tensoromics)))
  
  # Get families (already filtered to match expression genes)
  families <- filtered_data$families
  
  # Get gene IDs from expression matrix (Uniprot IDs)
  gene_ids <- colnames(expr_tensoromics)
  
  # Filter families to genes present in expression data
  families_filtered <- families[families$Protein %in% gene_ids, ]
  
  if (nrow(families_filtered) == 0) {
    cat(sprintf("  Warning: No family assignments for genes in %s - %s\n", project_id, stage_name))
    return(NULL)
  }
  
  # Create family arrays (from previous script)
  family_arrays <- create_family_arrays(families_filtered, gene_ids)
  
  # Count genes with family assignments
  genes_with_family <- sum(family_arrays$gene_to_fam > 0)
  cat(sprintf("  Genes with family assignments: %d/%d (%.1f%%)\n", 
              genes_with_family, length(gene_ids), 
              100 * genes_with_family / length(gene_ids)))
  
  # Prepare metadata for healthy samples
  if (data_type == "healthy" && !is.null(metadata)) {
    healthy_metadata <- metadata[sample_names, , drop = FALSE]
  } else {
    healthy_metadata <- NULL
  }
  
  # Prepare the final data structure
  prepared_data <- list(
    project_id = project_id,
    stage_name = stage_name,
    data_type = data_type,  # "stage" or "healthy"
    
    # TensorOmics format expression data
    expression_vectors = expr_tensoromics,  # samples x genes
    gene_ids = gene_ids,
    sample_names = rownames(expr_tensoromics),
    n_axes = nrow(expr_tensoromics),  # number of samples
    n_genes = ncol(expr_tensoromics),  # number of genes
    
    # Family information
    families = families_filtered,
    family_ids = family_arrays$family_ids,
    gene_to_fam = family_arrays$gene_to_fam,
    n_families = family_arrays$n_families,
    
    # Metadata (preserve original if available)
    metadata = list(
      stage_col = if(data_type == "stage") stage_col else NULL,
      tissue_type_col = tissue_type_col,
      n_samples = length(sample_names),
      sample_barcodes = sample_names,
      preparation_date = Sys.time(),
      data_type = if(data_type == "healthy") "Normal tissue (TPM)" else "Cancer stage (TPM)",
      format = "TensorOmics (samples × genes)",
      # Add original metadata if available
      original_metadata = if(!is.null(metadata) && data_type == "stage") {
        metadata[sample_names, ]
      } else if (data_type == "healthy") {
        healthy_metadata
      } else NULL,
      # For healthy data, add tissue type info
      tissue_types = if(data_type == "healthy" && !is.null(healthy_metadata)) {
        unique(healthy_metadata[[tissue_type_col]])
      } else NULL,
      # Add primary site info if available
      primary_site = if(!is.null(primary_site_filter)) {
        primary_site_filter$value
      } else NULL
    )
  )
  
  # Validate the data structure
  cat("  Validating data structure... ")
  
  # Check dimensions
  if (nrow(prepared_data$expression_vectors) != prepared_data$n_axes) {
    stop("Dimension mismatch: nrow(expression_vectors) != n_axes")
  }
  if (ncol(prepared_data$expression_vectors) != prepared_data$n_genes) {
    stop("Dimension mismatch: ncol(expression_vectors) != n_genes")
  }
  if (length(prepared_data$gene_ids) != prepared_data$n_genes) {
    stop("Dimension mismatch: length(gene_ids) != n_genes")
  }
  if (length(prepared_data$sample_names) != prepared_data$n_axes) {
    stop("Dimension mismatch: length(sample_names) != n_axes")
  }
  if (length(prepared_data$gene_to_fam) != prepared_data$n_genes) {
    stop("Dimension mismatch: length(gene_to_fam) != n_genes")
  }
  
  cat("OK\n")
  
  return(prepared_data)
}

#' Save prepared data to RDS file in results directory
#' @param prepared_data Prepared data list
#' @param output_dir Base output directory (results/)
#' @param suffix Optional suffix for filename
save_prepared_data <- function(prepared_data, output_dir, suffix = "", tpms = TRUE) {
  
  if (is.null(prepared_data)) return(NULL)
  
  # Create project directory
  project_dir <- file.path(output_dir, prepared_data$project_id)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create filename
  if(tpms){
    if (prepared_data$data_type == "stage") {
      filename <- paste0(prepared_data$project_id, "-", 
                        gsub(" ", "-", prepared_data$stage_name), 
                        suffix, ".rds")
    } else if (prepared_data$data_type == "healthy") {
      filename <- paste0("healthy_", prepared_data$project_id, 
                        suffix, ".rds")
    } else {
      filename <- paste0(prepared_data$project_id, "_", 
                        prepared_data$data_type, 
                        suffix, ".rds")
    }
  } else{
    if (prepared_data$data_type == "stage") {
      filename <- paste0(prepared_data$project_id, "-", 
                        gsub(" ", "-", prepared_data$stage_name), 
                        suffix, "_counts.rds")
    } else if (prepared_data$data_type == "healthy") {
      filename <- paste0("healthy_", prepared_data$project_id, 
                        suffix, "_counts.rds")
    } else {
      filename <- paste0(prepared_data$project_id, "_", 
                        prepared_data$data_type, 
                        suffix, "_counts.rds")
    }
  }
  
  
  # Save the data
  output_file <- file.path(project_dir, filename)
  saveRDS(prepared_data, output_file)
  
  cat(sprintf("  Saved to: %s\n", output_file))
  
  return(output_file)
}

#' Process healthy tissue data for a specific cancer project
#' @param project_id Cancer project ID
#' @param input_dir Input directory
#' @param output_dir Output directory
#' @param tissue_type_col Column name for tissue type
#' @param primary_site_filter Optional list with primary site filtering parameters
process_healthy_data <- function(project_id, input_dir = input_dir, 
                                output_dir = output_dir, 
                                tissue_type_col = tissue_type_col,
                                primary_site_filter = NULL,
                                tpms = TRUE) {
  
  cat(sprintf("\nProcessing healthy tissue for: %s\n", project_id))
  
  # Prepare healthy data
  healthy_data <- prepare_project_data(
    project_id = project_id,
    stage_name = "healthy",
    input_dir = input_dir,
    stage_col = "simple_stage",
    tissue_type_col = tissue_type_col,
    data_type = "healthy",
    primary_site_filter = primary_site_filter,
    tpms = tpms
  )
  
  if (!is.null(healthy_data)) {
    # Save to RDS file
    saved_file <- save_prepared_data(healthy_data, output_dir, "", tpms)
    
    return(list(
      data = healthy_data,
      file = saved_file,
      n_samples = healthy_data$n_axes,
      n_genes = healthy_data$n_genes,
      n_families = healthy_data$n_families,
      genes_with_family = sum(healthy_data$gene_to_fam > 0)
    ))
  }
  
  return(NULL)
}

#' Process a single cancer project for all stages and healthy tissue
#' @param project_id Cancer project ID
#' @param stages Vector of stages to process
#' @param input_dir Input directory
#' @param output_dir Output directory
#' @param stage_col Stage column name
#' @param tissue_type_col Tissue type column name
#' @param primary_site_filter Optional list with primary site filtering parameters
process_project <- function(project_id, stages = stages, input_dir = input_dir, 
                           output_dir = output_dir, stage_col = "simple_stage",
                           tissue_type_col = "sample_type",
                           primary_site_filter = NULL,
                           tpms = TRUE) {
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("PROCESSING: %s\n", project_id))
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  results <- list()
  
  # Process each stage
  for (stage in stages) {
    # Prepare data for this stage
    prepared_data <- prepare_project_data(
      project_id = project_id,
      stage_name = stage,
      input_dir = input_dir,
      stage_col = stage_col,
      tissue_type_col = tissue_type_col,
      data_type = "stage",
      primary_site_filter = primary_site_filter,
      tpms = tpms
    )
    
    if (!is.null(prepared_data)) {
      # Save to RDS file
      saved_file <- save_prepared_data(prepared_data, output_dir, "", tpms)
      
      # Store in results
      results[[stage]] <- list(
        data = prepared_data,
        file = saved_file,
        n_samples = prepared_data$n_axes,
        n_genes = prepared_data$n_genes,
        n_families = prepared_data$n_families,
        genes_with_family = sum(prepared_data$gene_to_fam > 0)
      )
    }
  }
  
  # Process healthy tissue data
  healthy_results <- process_healthy_data(
    project_id = project_id,
    input_dir = input_dir,
    output_dir = output_dir,
    tissue_type_col = tissue_type_col,
    primary_site_filter = primary_site_filter,
    tpms = tpms
  )
  
  if (!is.null(healthy_results)) {
    results$healthy <- healthy_results
  }
  
  # Print summary
  cat(sprintf("\nSummary for %s:\n", project_id))
  
  # Stage data summary
  for (stage in stages) {
    if (stage %in% names(results)) {
      cat(sprintf("  %s: %d samples, %d genes (%d with family), %d families\n", 
                  stage, 
                  results[[stage]]$n_samples,
                  results[[stage]]$n_genes,
                  results[[stage]]$genes_with_family,
                  results[[stage]]$n_families))
    } else {
      cat(sprintf("  %s: No data\n", stage))
    }
  }
  
  # Healthy data summary
  if ("healthy" %in% names(results)) {
    cat(sprintf("  Healthy: %d samples, %d genes (%d with family), %d families\n",
                results$healthy$n_samples,
                results$healthy$n_genes,
                results$healthy$genes_with_family,
                results$healthy$n_families))
  } else {
    cat(sprintf("  Healthy: No data\n"))
  }
  
  return(results)
}

#' Process all cancer projects (stages + healthy tissue)
#' @param projects Vector of project IDs
#' @param stages Vector of stages
#' @param input_dir Input directory
#' @param output_dir Output directory
process_all_projects <- function(projects = projects, stages = stages,
                                input_dir = input_dir, output_dir = output_dir,
                                tpms = TRUE) {
  
  cat("STARTING DATA PREPARATION\n")
  cat(sprintf("Input directory: %s\n", input_dir))
  cat(sprintf("Output directory: %s\n", output_dir))
  cat(sprintf("Projects to process: %d\n", length(projects)))
  cat(sprintf("Stages to process: %s\n", paste(stages, collapse = ", ")))
  cat(sprintf("Tissue type column: %s\n", tissue_type_col))
  cat("\n")
  
  all_results <- list()
  stages_local <- c("Stage I", "Stage II", "Stage III", "Stage IV") 
  for (project_id in projects) {
    project_results <- process_project(
      project_id = project_id,
      stages = stages_local,
      input_dir = input_dir,
      output_dir = output_dir,
      tpms = tpms
    )   
    all_results[[project_id]] = project_results
  }
  
  # Create summary statistics
  summary_df <- data.frame()
  healthy_summary_df <- data.frame()
  
  for (project_id in names(all_results)) {
    # Stage data
    for (stage in stages) {
      if (stage %in% names(all_results[[project_id]])) {
        result <- all_results[[project_id]][[stage]]
        
        summary_row <- data.frame(
          Project = project_id,
          Data_Type = "stage",
          Stage = stage,
          Samples = result$n_samples,
          Genes_Total = result$n_genes,
          Genes_With_Family = result$genes_with_family,
          Families = result$n_families,
          Percent_Genes_With_Family = round(100 * result$genes_with_family / result$n_genes, 1),
          File = basename(result$file),
          Path = dirname(result$file)
        )
        
        summary_df <- rbind(summary_df, summary_row)
      }
    }
    
    # Healthy data
    if ("healthy" %in% names(all_results[[project_id]])) {
      healthy_result <- all_results[[project_id]]$healthy
      
      healthy_row <- data.frame(
        Project = project_id,
        Data_Type = "healthy",
        Stage = "Normal",
        Samples = healthy_result$n_samples,
        Genes_Total = healthy_result$n_genes,
        Genes_With_Family = healthy_result$genes_with_family,
        Families = healthy_result$n_families,
        Percent_Genes_With_Family = round(100 * healthy_result$genes_with_family / healthy_result$n_genes, 1),
        File = basename(healthy_result$file),
        Path = dirname(healthy_result$file)
      )
      
      healthy_summary_df <- rbind(healthy_summary_df, healthy_row)
    }
  }
  
  # Save summaries
  summary_file <- file.path(output_dir, "data_preparation_summary.csv")
  write.csv(summary_df, summary_file, row.names = FALSE)
  cat(sprintf("\nStage data summary saved to: %s\n", summary_file))
  
  healthy_summary_file <- file.path(output_dir, "healthy_data_summary.csv")
  write.csv(healthy_summary_df, healthy_summary_file, row.names = FALSE)
  cat(sprintf("Healthy data summary saved to: %s\n", healthy_summary_file))
  
  # Combined summary
  combined_summary <- rbind(summary_df, healthy_summary_df)
  combined_summary_file <- file.path(output_dir, "combined_data_summary.csv")
  write.csv(combined_summary, combined_summary_file, row.names = FALSE)
  cat(sprintf("Combined summary saved to: %s\n", combined_summary_file))
  
  # Print overall summary
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("DATA PREPARATION COMPLETE\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  cat(sprintf("Total stage files created: %d\n", nrow(summary_df)))
  cat(sprintf("Total healthy files created: %d\n", nrow(healthy_summary_df)))
  cat(sprintf("Total samples (stage): %d\n", sum(summary_df$Samples)))
  cat(sprintf("Total samples (healthy): %d\n", sum(healthy_summary_df$Samples)))
  cat(sprintf("Total samples (all): %d\n", sum(combined_summary$Samples)))
  
  if (nrow(healthy_summary_df) > 0) {
    cat(sprintf("\nHealthy tissue statistics:\n"))
    cat(sprintf("  Projects with healthy data: %d/%d (%.1f%%)\n",
                nrow(healthy_summary_df), length(projects),
                100 * nrow(healthy_summary_df) / length(projects)))
    cat(sprintf("  Average healthy samples per project: %.1f\n", 
                mean(healthy_summary_df$Samples)))
    cat(sprintf("  Min healthy samples: %d\n", min(healthy_summary_df$Samples)))
    cat(sprintf("  Max healthy samples: %d\n", max(healthy_summary_df$Samples)))
  } else {
    cat(sprintf("\nWarning: No healthy tissue data found for any project!\n"))
    cat(sprintf("  Check tissue_type_col: '%s'\n", tissue_type_col))
    cat(sprintf("  Looking for value: 'Tissue Type Normal'\n"))
  }
  
  # Create project-wise summary
  project_summary <- combined_summary %>%
    group_by(Project, Data_Type) %>%
    summarise(
      Files = n(),
      Total_Samples = sum(Samples),
      Avg_Genes = mean(Genes_Total),
      Avg_Families = mean(Families)
    )
  
  if(tpms){
    project_summary_file <- file.path(output_dir, "project_summary.csv")
  } else{
    project_summary_file <- file.path(output_dir, "project_summary_counts.csv")
  }
  
  write.csv(project_summary, project_summary_file, row.names = FALSE)
  cat(sprintf("Project summary saved to: %s\n", project_summary_file))
  
  return(list(
    results = all_results,
    stage_summary = summary_df,
    healthy_summary = healthy_summary_df,
    combined_summary = combined_summary,
    project_summary = project_summary
  ))
}

# ==================== MAIN EXECUTION ====================

if (sys.nframe() == 0) {
  # This code runs when script is executed directly
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Process all projects (stages + healthy)
  results_tpm <- process_all_projects(
    projects = projects,
    stages = stages,
    input_dir = input_dir,
    output_dir = output_dir
  )

  results_counts <- process_all_projects(
    projects = projects,
    stages = stages,
    input_dir = input_dir,
    output_dir = output_dir,
    tpms = FALSE
  )
  
  # List all created files
  cat("\nAll created files:\n")
  
  # List RDS files in output directory
  rds_files <- list.files(output_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  
  if (length(rds_files) > 0) {
    for (file in rds_files) {
      file_info <- file.info(file)
      file_size <- file_info$size / 1024 / 1024  # Convert to MB
      cat(sprintf("  %s (%.2f MB)\n", file, file_size))
    }
    cat(sprintf("\nTotal files created: %d\n", length(rds_files)))
  } else {
    cat("  No RDS files found in output directory\n")
  }
}
