# Tensor Omics in Cancer Transcriptomics

This repository contains the code and selected results for the bachelor thesis:

"Family-based comparison of gene expression reveals novel candidates in cancer transcriptomics"

The project evaluates the Tensor Omics framework, a geometry-based approach for analyzing gene expression, and compares it to classical differential expression methods (edgeR and limma) using TCGA datasets.

---

## Repository structure

- `material/`
  Input data required for the analysis (e.g., ortholog mappings and gene family assignments).
  Large TCGA datasets are stored externally (see Data Availability).

- `methods/`
  All analysis code:
  - `data_preparation/` – data download, preprocessing, and filtering
  - `data_analysis/` – trajectory analysis and pipeline execution

- `results/`
  Generated figures and tables, including:
  - enrichment analysis results
  - quality control reports
  - trajectory analysis results

- `thesis/`
  Final thesis document (PDF)

---

## Data availability

Large intermediate files (e.g., `.rds` objects) and raw TCGA datasets are not included in this repository due to size constraints.

- TCGA data can be obtained from:
  https://portal.gdc.cancer.gov/

- Additional intermediate files are available via internal storage or upon request (see Data Availability section in the thesis).

---

## Requirements

The analysis was performed using:

- R (≥ 4.3)
- Bioconductor packages:
  - `edgeR`
  - `limma`
  - `TCGAbiolinks`
  - `clusterProfiler`
- Additional R packages listed in the thesis
- Fortran compiler (`gfortran`) for Tensor Omics components

---

## Reproducibility

To reproduce the analysis:

- Get necessary input files from the materials folder
- Run scripts in the following order:
1. tcga_data_pipeline.R
2. tcga_qc_enhanced.R
3. map_families.R
4. parse_orthologs.R
5. prepare_stage_data.R
6. family_trajectory_analysis.R
7. gene_trajectory_analysis.R
8. outlier_significance_analysis_fortran.R
9. deg_comparison.R
10. compare_tox_with_deg.R
11. pathway_analysis.R
12. cancer_relation_analysis.R

Due to the size of intermediate data, full reproducibility requires access to the datasets described in the thesis.

---

## Citation

If you use this work, please cite:

Aaron Schroeder, Asis Hallab, Vivian Bass Vega (2026)
*Family-based comparison of gene expression reveals novel candidates in cancer transcriptomics*
Bachelor Thesis, conducted at Bingen University of Applied Sciences

---

## License

All original code in this repository is provided under the MIT License.
Figures and results may be reused with proper attribution.
