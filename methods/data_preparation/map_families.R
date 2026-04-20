library(biomaRt)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)
library(stringr)

material_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/material"
all_proteins_gene_families <- "all_proteins_gene_families.txt"
result_dir <- "/media/BioNAS2/bachelor_thesis_aaron_schroeder/results/TCGA_QC_Reports"
project_ids <- c("TCGA-LUAD", "TCGA-BRCA", "TCGA-COAD", "TCGA-BLCA", 
                 "TCGA-KIRC", "TCGA-LUSC", "TCGA-STAD", "TCGA-THCA")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

protein_mapping_with_stage <- function(data, project_id = NULL, get_tpm = TRUE) {
  
  # Extract expression data
  if (is(data$enhanced_data, "SummarizedExperiment")) {
    if(get_tpm){
      expr_data <- assay(data$enhanced_data, "tpm_unstrand")
    } else {
      expr_data <- assay(data$enhanced_data, "unstranded")
    }
    metadata <- as.data.frame(colData(data$enhanced_data))
  } else {
    stop("Expected SummarizedExperiment in enhanced_data")
  }
  
  cat("  Starting protein mapping...\n")
  cat(sprintf("    Initial: %d genes × %d samples\n", nrow(expr_data), ncol(expr_data)))
  
  # 1. Filter rows with non-zero expression
  expr_data <- expr_data[rowSums(expr_data != 0) > 0, ]
  cat(sprintf("    After zero-filter: %d genes\n", nrow(expr_data)))
  
  # 2. Remove versioning from Ensembl IDs
  clean_ids <- str_remove(rownames(expr_data), "\\..*$")
  rownames(expr_data) <- clean_ids
  
  # 3. Get Uniprot mappings
  cat("    Querying Ensembl for mappings...\n")
  map <- getBM(
    attributes = c("ensembl_gene_id", "uniprotswissprot", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = clean_ids,
    mart = ensembl,
    uniqueRows = TRUE
  )
  map <- map[map$uniprotswissprot != "", ]
  cat(sprintf("    Mapped to Uniprot: %d genes\n", nrow(map)))
  
  # 4. Read gene families
  cat("    Reading gene families...\n")
  gene_families <- read.delim(file.path(material_dir, all_proteins_gene_families), 
                              header = FALSE, stringsAsFactors = FALSE) %>%
    dplyr::rename(Familie = V1) %>%
    tidyr::separate_rows(V2, sep = ",") %>%
    dplyr::mutate(
      Protein = sapply(strsplit(V2, "\\|"), function(x) {
        if (length(x) >= 2 && x[2] != "") x[2] else x[1]
      })
    ) %>%
    dplyr::select(Familie, Protein) %>%
    dplyr::filter(Protein != "" & !is.na(Protein))
  
  cat(sprintf("    Gene families: %d entries\n", nrow(gene_families)))
  
  # 5. Find common proteins
  common_proteins <- intersect(map$uniprotswissprot, gene_families$Protein)
  cat(sprintf("    Common proteins: %d\n", length(common_proteins)))
  
  if (length(common_proteins) == 0) {
    warning("No common proteins between mappings and families!")
    return(NULL)
  }
  
  # 6. Filter and order
  map_filtered <- map[map$uniprotswissprot %in% common_proteins, ]
  # Remove duplicates (keep first if multiple Ensembl → same Uniprot)
  map_filtered <- map_filtered[!duplicated(map_filtered$uniprotswissprot), ]
  
  families_filtered <- gene_families[gene_families$Protein %in% common_proteins, ]
  
  # 7. Filter expression data
  expr_filtered <- expr_data[rownames(expr_data) %in% map_filtered$ensembl_gene_id, , drop = FALSE]
  expr_filtered <- expr_filtered[map_filtered$ensembl_gene_id, , drop = FALSE]
  
  # Convert rownames to Uniprot IDs
  rownames(expr_filtered) <- map_filtered$uniprotswissprot
  
  # 8. Order families to match expression
  families_ordered <- families_filtered[match(rownames(expr_filtered), families_filtered$Protein), ]
  
  # 9. Return comprehensive result
  result <- list(
    # Expression data (with Uniprot IDs as rownames)
    dataframe = expr_filtered,
    
    # Families
    families = families_ordered,
    
    # METADATA - THIS IS NEW AND CRITICAL
    metadata = list(
      # Keep essential sample metadata
      sample_metadata = metadata,
      
      # Stage distribution
      stage_counts = table(metadata$simple_stage),
      
      # Mapping info
      id_mapping = data.frame(
        ensembl_gene_id = map_filtered$ensembl_gene_id,
        uniprot_id = map_filtered$uniprotswissprot,
        gene_symbol = map_filtered$hgnc_symbol,
        stringsAsFactors = FALSE
      ),
      
      # Statistics
      statistics = list(
        n_genes_original = nrow(expr_data),
        n_genes_filtered = nrow(expr_filtered),
        n_samples = ncol(expr_filtered),
        n_families = length(unique(families_ordered$Familie)),
        genes_per_family = table(families_ordered$Familie)
      )
    )
  )
  
  cat("    ✓ Mapping complete\n")
  cat(sprintf("      Final: %d proteins × %d samples, %d families\n", 
              nrow(expr_filtered), ncol(expr_filtered), 
              length(unique(families_ordered$Familie))))
  
  return(result)
}

# Enhanced filtering function that preserves stage information

run_all_filters_with_stage <- function(project_ids, directory, get_tpm = TRUE) {
  
  if (is.null(material_dir) || is.null(all_proteins_gene_families)) {
    stop("material_dir and all_proteins_gene_families must be specified")
  }
  
  results <- list()
  
  for(project in project_ids) {
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("PROCESSING PROJECT:", project, "\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    
    tryCatch({
      # Load the full enhanced data
      file <- file.path(directory, paste0(project, "_enhanced_QC.rds"))
      
      if (!file.exists(file)) {
        warning("File not found: ", file)
        next
      }
      
      cat("Reading:", file, "\n")
      data <- readRDS(file)
      
      # Check for stage information
      if (is(data$enhanced_data, "SummarizedExperiment")) {
        metadata <- colData(data$enhanced_data)
        
        if (!"simple_stage" %in% colnames(metadata)) {
          warning(sprintf("No 'simple_stage' column in %s metadata!", project))
          next
        }
        
        # Show stage distribution
        stages <- table(metadata$simple_stage)
        cat("  Stage distribution:\n")
        for(stage in names(stages)) {
          cat(sprintf("    %s: %d samples\n", stage, stages[stage]))
        }
      }

      cat("  Running protein mapping with stage preservation...\n")
      mapping <- protein_mapping_with_stage(data, project, get_tpm = get_tpm)
      
      if (is.null(mapping)) {
        warning(sprintf("Mapping failed for %s", project))
        next
      }
      
      # Save the filtered data WITH metadata
      if(get_tpm){
        output_file <- file.path(directory, paste0(project, "_filtered.rds"))
        saveRDS(mapping, file = output_file)
        
        # Also save a CSV version of expression data
        csv_file <- file.path(directory, paste0(project, "_filtered_expression.csv"))
        write.csv(as.data.frame(mapping$dataframe), csv_file)
      } else{
        output_file <- file.path(directory, paste0(project, "_filtered_counts.rds"))
        saveRDS(mapping, file = output_file)
        
        # Also save a CSV version of expression data
        csv_file <- file.path(directory, paste0(project, "_filtered_expression_counts.csv"))
        write.csv(as.data.frame(mapping$dataframe), csv_file)
      }
      
      
      cat("  ✓ Saved to:", output_file, "\n")
      cat("  ✓ CSV saved to:", csv_file, "\n")
      
      # Store in results
      results[[project]] <- list(
        project = project,
        mapping = mapping,
        stats = mapping$metadata$statistics,
        files = list(rds = output_file, csv = csv_file)
      )
      
      # Print summary
      cat("  Summary:\n")
      cat(sprintf("    Proteins: %d\n", mapping$metadata$statistics$n_genes_filtered))
      cat(sprintf("    Samples: %d\n", mapping$metadata$statistics$n_samples))
      cat(sprintf("    Families: %d\n", mapping$metadata$statistics$n_families))
      
    }, error = function(e) {
      cat("  ERROR processing", project, ":", e$message, "\n")
      results[[project]] <- list(
        project = project,
        error = e$message,
        status = "failed"
      )
    })
  }
  
  # Print overall summary
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("FILTERING SUMMARY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  
  successful <- sum(sapply(results, function(x) !is.null(x$stats)))
  cat("Projects processed:", length(project_ids), "\n")
  cat("Successfully filtered:", successful, "\n")
  
  if (successful > 0) {
    cat("\nDetailed statistics:\n")
    for (proj in names(results)) {
      if (!is.null(results[[proj]]$stats)) {
        stats <- results[[proj]]$stats
        cat(sprintf("\n- %s:\n", proj))
        cat(sprintf("  Proteins: %d\n", stats$n_genes_filtered))
        cat(sprintf("  Samples: %d\n", stats$n_samples))
        cat(sprintf("  Families: %d\n", stats$n_families))
        if (!is.null(stats$genes_per_family) && length(stats$genes_per_family) > 0) {
          cat(sprintf("  Avg genes/family: %.1f\n", 
                      mean(stats$genes_per_family)))
        }
      }
    }
  }
  
  # Save comprehensive summary
  if(get_tpm){
    summary_file <- file.path(directory, "filtering_with_stage_summary.rds")
  } else{
    summary_file <- file.path(directory, "filtering_with_stage_counts_summary.rds")
  }
  
  saveRDS(results, summary_file)
  cat("\n✓ Overall summary saved to:", summary_file, "\n")
  
  return(invisible(results))
}

# Run the enhanced filtering
results_tpm <- run_all_filters_with_stage(project_ids, directory = result_dir, get_tpm = TRUE)
results_counts <- run_all_filters_with_stage(project_ids, directory = result_dir, get_tpm = FALSE)