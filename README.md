# Bachelor-thesis-aaron

## Folder Structure
```
bachelor-thesis-aaron/
|- material/
    |- orthologs.tsv # Ortholog mapping
    # all used TCGA data is stored in this folder on the BioServer
|- methods/
    |- data_analysis/
        |- config.R # contains basic global configurations
        |- family_trajectory_analysis.R # Runs the family based trajectory analysis
        |- gene_trajectory_analysis.R # Runs the gene based trajectory analysis -> depends on results from `family_trajectory_analysis.R`
        |- plot_trajectory_results.R # Uses results saved by trajectory analysis scripts to generate plots
        |- run_pipeline.R # Runs the full analysis pipeline
        |- utils.R # Contains utility functions to load/save intermediate results
    |- data_preparation/
        |- map_families.R # Maps all ensembl IDs to uniprot IDs; Applies basic filtering
        |- parse_orthologs.R # Reads orthologs file and creates protein -> family mapping for each gene ID
        |- prepare_stage_data.R # Reads intermediate results from `tcga_qc_enhanced.R` and creates dataframes for each cancer and each stage
        |- tcga_data_pipeline.R # Downloads and extracts the full TCGA data for a given set of projects
        |- tcga_qc_enhanced.R # Creates quality control reports
    |- tox/
        |- error_handling.R # Contains error handling utilies for tensor omics
        |- tensoromics_functions.cpp # contains all tox rcpp bindings
        |- tensoromics_functions.R # Contains a selection of tensoromics functions used in the analysis
|- results/
    |- cancer_enrichment_three_sets/ # contains combined and per-cancer tables resulting from fishers exact test.
    |- TCGA_QC_Reports/
    # contains quality control plots
    |- 
|- thesis/ 
    # will contain the thesis later
|- README.md # This file
```