#!/bin/bash

rsync -avz \
    --prune-empty-dirs \
    --include='*/' \
    --include='*/TCGA-*/plots**' \
    --include='*/TCGA-*/jaccard**' \
    --include='differential_expression/*/*/*.png' \
    --include='pvalue_distributions/ecdf_plots/*.png' \
    --include='pathway_plots/*/*/*.png' \
    --include='comparison_de/*.png' \
    --include='diagnostic_plots_adaptive_knn/*.png' \
    --include='significant_gene_plots/*/*/*.png' \
    --include='significant_gene_plots/*.png' \
    --include='outlier_significance_adaptive_knn_fortran/*/*.png' \
    --include='multicancer_heatmaps*/*.png' \
    --include='jaccard_plots/*.png' \
    --include='cancer_enrichment_three_sets/**' \
    --exclude='*' \
    schroder@143.93.91.124:/media/BioNAS2/bachelor_thesis_aaron_schroeder/results/ \
    ../results/
