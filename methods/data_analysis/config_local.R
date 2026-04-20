# config_local.R
# Configuration for Family & Gene Trajectory Analysis
# Asign paths based on local machine

CONFIG <- list(
  # Base directories
  base_data_dir = "~/Schreibtisch/Bachelor_thesis/bachelor-thesis-aaron/results/",
  ortholog_file = "~/Schreibtisch/Bachelor_thesis/bachelor-thesis-aaron/results/ortholog_pairs.csv",
  
  # Projects to analyze
  projects = c("TCGA-BLCA", "TCGA-BRCA", "TCGA-COAD", #"TCGA-HNSC",
               "TCGA-KIRC", "TCGA-LUAD", "TCGA-LUSC", 
               "TCGA-STAD", "TCGA-THCA"),
  
  # Cancer stages
  stages = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  
  # Analysis parameters
  outlier_percentile = 95.0,
  loess_span = 0.7,
  loess_degree = 2,
  loess_mode = 1,
  loess_n_iters = 3,
  
  
  # Output directories
  family_output_dir = "~/Schreibtisch/Bachelor_thesis/bachelor-thesis-aaron/results/family_analysis",
  gene_output_dir = "~/Schreibtisch/Bachelor_thesis/bachelor-thesis-aaron/results/gene_analysis",
  de_output_dir = "~/Schreibtisch/Bachelor_thesis/bachelor-thesis-aaron/results/differential_expression",
  
  # TensorOmics function file
  tensoromics_functions_file = "tensoromics_functions.R",
  min_genes_per_orth_fam = 3,
  min_genes_per_all_family = 3
)

# Export configuration to global environment
assign("BASE_DATA_DIR", CONFIG$base_data_dir, envir = .GlobalEnv)
assign("ORTHOLOG_FILE", CONFIG$ortholog_file, envir = .GlobalEnv)
assign("PROJECTS", CONFIG$projects, envir = .GlobalEnv)
assign("STAGES", CONFIG$stages, envir = .GlobalEnv)
assign("OUTLIER_PERCENTILE", CONFIG$OUTLIER_PERCENTILE, envir = .GlobalEnv)
assign("FAMILY_OUTPUT_DIR", CONFIG$family_output_dir, envir = .GlobalEnv)
assign("GENE_OUTPUT_DIR", CONFIG$gene_output_dir, envir = .GlobalEnv)
assign("MIN_GENES_PER_ORTH_FAM", CONFIG$min_genes_per_orth_fam, envir = .GlobalEnv)
assign("MIN_GENES_PER_ALL_FAM", CONFIG$min_genes_per_all_family, envir = .GlobalEnv)
assign("DE_OUTPUT_DIR", CONFIG$de_output_dir, envir = .GlobalEnv)

cat("✓ Configuration loaded\n")