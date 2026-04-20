suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(httr2)
})

# ============================================================
# CONFIGURATION
# ============================================================

source("config.R")

# ---- cancer annotation files ----
CGC_FILE    <- file.path(MATERIAL_DIR, "cancerGenes/Cosmic_Genes_v103_GRCh38.tsv.gz")
NCG_FILE    <- file.path(MATERIAL_DIR, "cancerGenes/NCG_cancerdrivers_annotation_supporting_evidence.tsv")
ONCOKB_FILE <- file.path(MATERIAL_DIR, "cancerGenes/cancerGeneList.tsv")

# ---- input data ----
DE_BASE_DIR  <- DE_OUTPUT_DIR
TOX_RDS_FILE <- paste0(BASE_DATA_DIR, "outlier_significance_adaptive_knn_fortran/significant_all.rds")

# ---- output ----
OUTPUT_BASE_DIR <- paste0(BASE_DATA_DIR, "cancer_enrichment_three_sets")

# ---- cancers to analyze ----
CANCER_IDS <- c(
  "TCGA-BLCA",
  "TCGA-BRCA",
  "TCGA-KIRC",
  "TCGA-LUAD",
  "TCGA-LUSC",
  "TCGA-STAD",
  "TCGA-THCA",
  "TCGA-COAD"
)

# ---- stage folders for DGE ----
STAGE_DIRS <- c("Stage I", "Stage II", "Stage III", "Stage IV")

# ---- significance thresholds ----
EDGER_FDR_THRESHOLD <- 0.01
LIMMA_FDR_THRESHOLD <- 0.01
TOX_ADJ_P_THRESHOLD <- 0.01

# ---- UniProt mapping organism ----
ORGANISM_ID <- 9606
UNIPROT_PAGE_SIZE <- 500

# ============================================================
# PROVIDED HELPER
# ============================================================

generateContingencyTable <- function(row.1, row.2, col.1, col.2, row.category, col.category,
                                     categories = c("T", "F"), validate = TRUE) {
  cont.tbl <- matrix(
    c(length(intersect(row.1, col.1)),
      length(intersect(row.2, col.1)),
      length(intersect(row.1, col.2)),
      length(intersect(row.2, col.2))),
    nrow = 2, ncol = 2,
    dimnames = setNames(list(categories, categories), c(row.category, col.category))
  )

  if (validate) {
    val.res <- c(
      sum(cont.tbl[1, ]) == length(row.1),
      sum(cont.tbl[2, ]) == length(row.2),
      sum(cont.tbl[, 1]) == length(col.1),
      sum(cont.tbl[, 2]) == length(col.2),
      sum(cont.tbl) == length(unique(c(row.1, row.2, col.1, col.2)))
    )
    if (!all(val.res)) {
      stop("Validation of contingency table failed. Error-Code(s): ",
           paste(which(!val.res), collapse = ","))
    }
  }

  cont.tbl
}

# ============================================================
# HELPERS
# ============================================================

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

write_matrix_tsv <- function(mat, path) {
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df <- cbind(row_name = rownames(df), df, row.names = NULL)
  write_tsv(df, path)
}

chunk_vec <- function(x, size = 200) {
  split(x, ceiling(seq_along(x) / size))
}

fetch_uniprot_gene_map <- function(accessions, page_size = UNIPROT_PAGE_SIZE) {
  accessions <- unique(accessions)

  if (length(accessions) == 0) {
    return(data.frame(
      uniprot = character(),
      gene_symbol = character(),
      stringsAsFactors = FALSE
    ))
  }

  # 1) submit mapping job
  submit_resp <- request("https://rest.uniprot.org/idmapping/run") |>
    req_method("POST") |>
    req_body_form(
      ids = paste(accessions, collapse = ","),
      from = "UniProtKB_AC-ID",
      to = "UniProtKB"
    ) |>
    req_user_agent("R-cancer-enrichment-analysis") |>
    req_perform()

  submit_json <- resp_body_json(submit_resp)

  if (is.null(submit_json$jobId)) {
    stop("UniProt ID mapping submission failed: no jobId returned.")
  }

  job_id <- submit_json$jobId

  # 2) poll status
  status_url <- paste0("https://rest.uniprot.org/idmapping/status/", job_id)

  repeat {
    Sys.sleep(1)

    status_resp <- request(status_url) |>
      req_user_agent("R-cancer-enrichment-analysis") |>
      req_perform()

    status_json <- resp_body_json(status_resp)

    if (!is.null(status_json$results) || !is.null(status_json$failedIds)) {
      break
    }

    if (!is.null(status_json$jobStatus) &&
        status_json$jobStatus %in% c("RUNNING", "NEW")) {
      next
    }

    if (!is.null(status_json$jobStatus) &&
        status_json$jobStatus == "FAILED") {
      stop("UniProt ID mapping job failed.")
    }
  }

  # 3) fetch paginated TSV results
  next_url <- paste0(
    "https://rest.uniprot.org/idmapping/uniprotkb/results/",
    job_id,
    "?format=tsv&fields=accession,gene_primary&size=",
    page_size
  )

  all_pages <- list()
  first_colnames <- NULL

  repeat {
    resp <- request(next_url) |>
      req_user_agent("R-cancer-enrichment-analysis") |>
      req_perform()

    txt <- resp_body_string(resp)

    if (is.null(first_colnames)) {
      page_df <- read_tsv(I(txt), show_col_types = FALSE)
      first_colnames <- names(page_df)
    } else {
      page_df <- read_tsv(I(txt), show_col_types = FALSE, col_names = FALSE)
      names(page_df) <- first_colnames
    }

    all_pages[[length(all_pages) + 1]] <- page_df

    link_header <- resp_headers(resp)[["link"]]

    if (is.null(link_header)) {
      break
    }

    m <- regmatches(
      link_header,
      regexec('<([^>]+)>; rel="next"', link_header)
    )[[1]]

    if (length(m) >= 2) {
      next_url <- m[2]
    } else {
      break
    }
  }

  df <- bind_rows(all_pages)

  submitted <- data.frame(
    uniprot = accessions,
    stringsAsFactors = FALSE
  )

  mapped <- df %>%
    transmute(
      uniprot = as.character(From),
      gene_symbol = toupper(trimws(as.character(`Gene Names (primary)`)))
    ) %>%
    filter(!is.na(uniprot), uniprot != "") %>%
    distinct()

  out <- submitted %>%
    left_join(mapped, by = "uniprot")

  as.data.frame(out, stringsAsFactors = FALSE)
}

to_genes <- function(ids, id_map) {
  unique(id_map$gene_symbol[id_map$uniprot %in% ids & !is.na(id_map$gene_symbol)]) |>
    sort()
}

# ============================================================
# EXPLICIT DATABASE READERS
# ============================================================

read_cgc_genes <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  out <- df %>%
    transmute(
      gene_symbol = toupper(trimws(as.character(GENE_SYMBOL))),
      in_cancer_census = tolower(trimws(as.character(IN_CANCER_CENSUS)))
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "", in_cancer_census == "y") %>%
    distinct(gene_symbol) %>%
    mutate(is_CGC = TRUE)

  as.data.frame(out, stringsAsFactors = FALSE)
}

read_ncg_genes <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  out <- df %>%
    transmute(gene_symbol = toupper(trimws(as.character(symbol)))) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    distinct(gene_symbol) %>%
    mutate(is_NCG = TRUE)

  as.data.frame(out, stringsAsFactors = FALSE)
}

read_oncokb_genes <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  out <- df %>%
    transmute(
      gene_symbol = toupper(trimws(as.character(`Hugo Symbol`))),
      oncokb_annotated = trimws(as.character(`OncoKB Annotated`))
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    filter(is.na(oncokb_annotated) | oncokb_annotated != "") %>%
    distinct(gene_symbol) %>%
    mutate(is_OncoKB = TRUE)

  as.data.frame(out, stringsAsFactors = FALSE)
}

# ============================================================
# DATA EXTRACTION
# ============================================================

get_tox_groups_for_cancer <- function(tox_rds_file, cancer_id, adj_p_threshold = 0.05) {
  tox <- readRDS(tox_rds_file)

  tox_sub <- tox %>%
    filter(
      cancer_id == !!cancer_id,
      !is.na(gene_id),
      !is.na(noise_p_value_adj),
      noise_p_value_adj <= adj_p_threshold
    ) %>%
    mutate(
      norm_method = as.character(norm_method),
      comparison = as.character(comparison),
      tox_group = paste(norm_method, comparison, "TOX", sep = "_")
    )

  all_tox <- tox_sub %>%
    pull(gene_id) %>%
    unique() %>%
    sort()

  subgroup_tbl <- tox_sub %>%
    group_by(tox_group) %>%
    summarise(gene_ids = list(sort(unique(gene_id))), .groups = "drop")

  subgroup_list <- setNames(subgroup_tbl$gene_ids, subgroup_tbl$tox_group)

  c(list(ALL_TOX = all_tox), subgroup_list)
}

get_tox_groups_all_cancers <- function(tox_rds_file, cancer_ids, adj_p_threshold = 0.05) {
  tox <- readRDS(tox_rds_file)

  tox_sub <- tox %>%
    filter(
      cancer_id %in% cancer_ids,
      !is.na(gene_id),
      !is.na(noise_p_value_adj),
      noise_p_value_adj <= adj_p_threshold
    ) %>%
    mutate(
      norm_method = as.character(norm_method),
      comparison = as.character(comparison),
      tox_group = paste(norm_method, comparison, "TOX", sep = "_")
    )

  all_tox <- tox_sub %>%
    pull(gene_id) %>%
    unique() %>%
    sort()

  subgroup_tbl <- tox_sub %>%
    group_by(tox_group) %>%
    summarise(gene_ids = list(sort(unique(gene_id))), .groups = "drop")

  subgroup_list <- setNames(subgroup_tbl$gene_ids, subgroup_tbl$tox_group)

  c(list(ALL_TOX = all_tox), subgroup_list)
}

read_edger_file <- function(path, fdr_threshold = 0.01) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(gene_id), !is.na(FDR), FDR <= fdr_threshold) %>%
    pull(gene_id) %>%
    unique() %>%
    sort()
}

read_limma_file <- function(path, fdr_threshold = 0.01) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(gene_id), !is.na(adj.P.Val), adj.P.Val <= fdr_threshold) %>%
    pull(gene_id) %>%
    unique() %>%
    sort()
}

get_de_genes_for_cancer <- function(de_base_dir,
                                    cancer_id,
                                    stage_dirs = STAGE_DIRS,
                                    edger_fdr_threshold = 0.01,
                                    limma_fdr_threshold = 0.01) {
  project_dir <- file.path(de_base_dir, cancer_id)

  edger_all <- character()
  limma_all <- character()

  for (stage_name in stage_dirs) {
    stage_dir <- file.path(project_dir, stage_name)

    edger_path <- file.path(stage_dir, "edgeR_results.csv")
    limma_path <- file.path(stage_dir, "voom_results.csv")

    if (file.exists(edger_path)) {
      edger_all <- c(edger_all, read_edger_file(edger_path, edger_fdr_threshold))
    }

    if (file.exists(limma_path)) {
      limma_all <- c(limma_all, read_limma_file(limma_path, limma_fdr_threshold))
    }
  }

  list(
    edger = unique(edger_all) |> sort(),
    limma = unique(limma_all) |> sort()
  )
}

# ============================================================
# STATISTICS
# ============================================================

run_single_fisher <- function(group_name, label_name, detected_genes, cancer_genes, background_genes) {
  detected_genes <- sort(unique(intersect(detected_genes, background_genes)))
  not_detected   <- sort(setdiff(background_genes, detected_genes))

  cancer_genes_bg  <- sort(unique(intersect(cancer_genes, background_genes)))
  not_cancer_genes <- sort(setdiff(background_genes, cancer_genes_bg))

  tbl <- generateContingencyTable(
    row.1 = cancer_genes_bg,
    row.2 = not_cancer_genes,
    col.1 = detected_genes,
    col.2 = not_detected,
    row.category = "CancerRelated",
    col.category = "Detected",
    categories = c("T", "F"),
    validate = TRUE
  )

  ft <- fisher.test(tbl, alternative = "greater")

  a <- length(intersect(cancer_genes_bg, detected_genes))
  fp <- length(setdiff(detected_genes, cancer_genes_bg))

  list(
    result = data.frame(
      group = group_name,
      cancer_label = label_name,
      detected_n = length(detected_genes),
      cancer_detected_n = a,
      noncancer_detected_n = fp,
      precision = ifelse(length(detected_genes) > 0, a / length(detected_genes), NA_real_),
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    ),
    table = tbl
  )
}

# ============================================================
# CORE ANALYSIS FOR ONE GENE COLLECTION
# ============================================================

run_analysis_for_gene_sets <- function(tox_group_uniprot,
                                       edger_uniprot,
                                       limma_uniprot,
                                       result_dir,
                                       cgc_file = CGC_FILE,
                                       ncg_file = NCG_FILE,
                                       oncokb_file = ONCOKB_FILE) {
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  if (!("ALL_TOX" %in% names(tox_group_uniprot))) {
    stop("tox_group_uniprot must contain an 'ALL_TOX' entry.")
  }

  all_tox_uniprot <- tox_group_uniprot[["ALL_TOX"]]
  DGE_uniprot <- union(edger_uniprot, limma_uniprot)

  only_DGE_uniprot <- setdiff(DGE_uniprot, all_tox_uniprot)
  tox_and_DGE_uniprot <- intersect(all_tox_uniprot, DGE_uniprot)
  only_tox_uniprot <- setdiff(all_tox_uniprot, DGE_uniprot)

  background_uniprot <- sort(unique(c(all_tox_uniprot, edger_uniprot, limma_uniprot)))

  input_sizes <- data.frame(
    all_tox_n = length(all_tox_uniprot),
    edger_n = length(edger_uniprot),
    limma_n = length(limma_uniprot),
    all_DGE_n = length(DGE_uniprot),
    only_DGE_n = length(only_DGE_uniprot),
    tox_and_DGE_n = length(tox_and_DGE_uniprot),
    only_tox_n = length(only_tox_uniprot),
    background_n = length(background_uniprot),
    stringsAsFactors = FALSE
  )
  write_tsv(input_sizes, file.path(result_dir, "input_sizes.tsv"))

  writeLines(background_uniprot, file.path(result_dir, "background_uniprot_input_ids.txt"))

  id_map <- fetch_uniprot_gene_map(background_uniprot, page_size = UNIPROT_PAGE_SIZE)
  write_tsv(id_map, file.path(result_dir, "uniprot_to_gene_map.tsv"))

  mapped_ids <- id_map$uniprot[!is.na(id_map$gene_symbol)]
  unmapped_ids <- setdiff(background_uniprot, mapped_ids)
  writeLines(unmapped_ids, file.path(result_dir, "unmapped_uniprot_ids.txt"))

  background_symbols <- to_genes(background_uniprot, id_map)
  edger_symbols <- to_genes(edger_uniprot, id_map)
  limma_symbols <- to_genes(limma_uniprot, id_map)
  DGE_symbols <- to_genes(DGE_uniprot, id_map)

  tox_group_symbols <- lapply(tox_group_uniprot, function(x) to_genes(x, id_map))

  all_tox_symbols <- tox_group_symbols[["ALL_TOX"]]
  only_DGE_symbols <- setdiff(DGE_symbols, all_tox_symbols)
  tox_and_DGE_symbols <- intersect(all_tox_symbols, DGE_symbols)
  only_tox_symbols <- setdiff(all_tox_symbols, DGE_symbols)

  cgc_tbl    <- read_cgc_genes(cgc_file)
  ncg_tbl    <- read_ncg_genes(ncg_file)
  oncokb_tbl <- read_oncokb_genes(oncokb_file)

  annot_tbl <- data.frame(gene_symbol = sort(unique(background_symbols)), stringsAsFactors = FALSE) %>%
    left_join(cgc_tbl, by = "gene_symbol") %>%
    left_join(ncg_tbl, by = "gene_symbol") %>%
    left_join(oncokb_tbl, by = "gene_symbol") %>%
    mutate(
      is_CGC = ifelse(is.na(is_CGC), FALSE, is_CGC),
      is_NCG = ifelse(is.na(is_NCG), FALSE, is_NCG),
      is_OncoKB = ifelse(is.na(is_OncoKB), FALSE, is_OncoKB),
      is_MERGED = is_CGC | is_NCG | is_OncoKB
    )

  write_tsv(annot_tbl, file.path(result_dir, "background_gene_annotations.tsv"))

  gene_order <- sort(unique(background_symbols))

  detection_tbl <- data.frame(
    gene_symbol = gene_order,
    detected_ALL_edgeR = gene_order %in% edger_symbols,
    detected_ALL_limma = gene_order %in% limma_symbols,
    detected_ALL_DGE = gene_order %in% DGE_symbols,
    detected_ALL_TOX = gene_order %in% all_tox_symbols,
    detected_ONLY_DGE = gene_order %in% only_DGE_symbols,
    detected_TOX_AND_DGE = gene_order %in% tox_and_DGE_symbols,
    detected_ONLY_TOX = gene_order %in% only_tox_symbols,
    stringsAsFactors = FALSE
  )

  # Add TOX subgroup detection columns
  tox_subgroup_names <- setdiff(names(tox_group_symbols), "ALL_TOX")
  if (length(tox_subgroup_names) > 0) {
    for (grp in tox_subgroup_names) {
      col_name <- paste0("detected_", make.names(grp))
      detection_tbl[[col_name]] <- detection_tbl$gene_symbol %in% tox_group_symbols[[grp]]
    }
  }

  detection_tbl <- detection_tbl %>%
    left_join(annot_tbl, by = "gene_symbol")

  write_tsv(detection_tbl, file.path(result_dir, "gene_detection_and_annotations.tsv"))

  cgc_genes    <- annot_tbl %>% filter(is_CGC) %>% pull(gene_symbol)
  ncg_genes    <- annot_tbl %>% filter(is_NCG) %>% pull(gene_symbol)
  oncokb_genes <- annot_tbl %>% filter(is_OncoKB) %>% pull(gene_symbol)
  merged_genes <- annot_tbl %>% filter(is_MERGED) %>% pull(gene_symbol)

  cancer_sets <- list(
    CGC = cgc_genes,
    NCG = ncg_genes,
    OncoKB = oncokb_genes,
    MERGED = merged_genes
  )

  only_tox_hits_tbl <- data.frame(
    uniprot = sort(unique(only_tox_uniprot)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(id_map, by = c("uniprot" = "uniprot")) %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
    dplyr::left_join(annot_tbl, by = "gene_symbol") %>%
    dplyr::mutate(
      is_CGC = ifelse(is.na(is_CGC), FALSE, is_CGC),
      is_NCG = ifelse(is.na(is_NCG), FALSE, is_NCG),
      is_OncoKB = ifelse(is.na(is_OncoKB), FALSE, is_OncoKB),
      is_MERGED = ifelse(is.na(is_MERGED), FALSE, is_MERGED)
    ) %>%
    dplyr::distinct()

  readr::write_tsv(
    only_tox_hits_tbl,
    file.path(result_dir, "only_tox_database_hits.tsv")
  )

  group_sets <- c(
    list(
      ALL_edgeR = edger_symbols,
      ALL_limma = limma_symbols,
      ALL_DGE = DGE_symbols
    ),
    tox_group_symbols,
    list(
      ONLY_DGE = only_DGE_symbols,
      TOX_AND_DGE = tox_and_DGE_symbols,
      ONLY_TOX = only_tox_symbols
    )
  )

  result_rows <- list()

  for (group_name in names(group_sets)) {
    for (label_name in names(cancer_sets)) {
      out <- run_single_fisher(
        group_name = group_name,
        label_name = label_name,
        detected_genes = group_sets[[group_name]],
        cancer_genes = cancer_sets[[label_name]],
        background_genes = background_symbols
      )

      key <- paste(group_name, label_name, sep = "_")
      result_rows[[key]] <- out$result

      write_matrix_tsv(
        out$table,
        file.path(result_dir, paste0("contingency_", safe_name(group_name), "_", label_name, ".tsv"))
      )
    }
  }

  fisher_results <- bind_rows(result_rows) %>%
    mutate(
      p_adjust_BH = p.adjust(p_value, method = "BH"),
      p_adjust_BY = p.adjust(p_value, method = "BY"),
      precision_percent = round(100 * precision, 2),
      tool = case_when(
        grepl("_TOX$", group) ~ "TOX",
        group == "ALL_edgeR" ~ "edgeR",
        group == "ALL_limma" ~ "limma",
        group == "ALL_DGE" ~ "DGE",
        group == "ONLY_DGE" ~ "DGE_only",
        group == "TOX_AND_DGE" ~ "shared",
        group == "ONLY_TOX" ~ "TOX_only",
        TRUE ~ "other"
      ),
      norm_method = case_when(
        grepl("^std_log_", group) ~ "std_log",
        grepl("^full_", group) ~ "full",
        grepl("^log_", group) ~ "log",
        TRUE ~ NA_character_
      ),
      comparison_method = case_when(
        grepl("_family_mean_TOX$", group) ~ "family_mean",
        grepl("_ortholog_mean_TOX$", group) ~ "ortholog_mean",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(group, p_value)

  write_tsv(fisher_results, file.path(result_dir, "fisher_results.tsv"))

  invisible(list(
    fisher_results = fisher_results,
    annotations = annot_tbl,
    detection = detection_tbl
  ))
}

# ============================================================
# PER-CANCER ANALYSIS
# ============================================================

run_cancer_analysis <- function(cancer_id,
                                cgc_file = CGC_FILE,
                                ncg_file = NCG_FILE,
                                oncokb_file = ONCOKB_FILE,
                                de_base_dir = DE_BASE_DIR,
                                tox_rds_file = TOX_RDS_FILE,
                                output_base_dir = OUTPUT_BASE_DIR) {
  result_dir <- file.path(output_base_dir, safe_name(cancer_id))
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  tox_group_uniprot <- get_tox_groups_for_cancer(
    tox_rds_file = tox_rds_file,
    cancer_id = cancer_id,
    adj_p_threshold = TOX_ADJ_P_THRESHOLD
  )

  de_lists <- get_de_genes_for_cancer(
    de_base_dir = de_base_dir,
    cancer_id = cancer_id,
    stage_dirs = STAGE_DIRS,
    edger_fdr_threshold = EDGER_FDR_THRESHOLD,
    limma_fdr_threshold = LIMMA_FDR_THRESHOLD
  )

  edger_uniprot <- de_lists$edger
  limma_uniprot <- de_lists$limma

  res <- run_analysis_for_gene_sets(
    tox_group_uniprot = tox_group_uniprot,
    edger_uniprot = edger_uniprot,
    limma_uniprot = limma_uniprot,
    result_dir = result_dir,
    cgc_file = cgc_file,
    ncg_file = ncg_file,
    oncokb_file = oncokb_file
  )

  fisher_results <- res$fisher_results
  fisher_results$cancer_id <- cancer_id

  write_tsv(
    fisher_results %>% select(cancer_id, everything()),
    file.path(result_dir, "fisher_results_with_cancer_id.tsv")
  )

  invisible(fisher_results %>% select(cancer_id, everything()))
}

# ============================================================
# COMBINED ANALYSIS ACROSS ALL CANCERS
# ============================================================

run_combined_analysis <- function(cancer_ids = CANCER_IDS,
                                  cgc_file = CGC_FILE,
                                  ncg_file = NCG_FILE,
                                  oncokb_file = ONCOKB_FILE,
                                  de_base_dir = DE_BASE_DIR,
                                  tox_rds_file = TOX_RDS_FILE,
                                  output_base_dir = OUTPUT_BASE_DIR) {
  result_dir <- file.path(output_base_dir, "ALL_CANCERS_COMBINED")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  tox_group_uniprot <- get_tox_groups_all_cancers(
    tox_rds_file = tox_rds_file,
    cancer_ids = cancer_ids,
    adj_p_threshold = TOX_ADJ_P_THRESHOLD
  )

  edger_all <- character()
  limma_all <- character()

  for (cid in cancer_ids) {
    de_tmp <- get_de_genes_for_cancer(
      de_base_dir = de_base_dir,
      cancer_id = cid,
      stage_dirs = STAGE_DIRS,
      edger_fdr_threshold = EDGER_FDR_THRESHOLD,
      limma_fdr_threshold = LIMMA_FDR_THRESHOLD
    )

    edger_all <- c(edger_all, de_tmp$edger)
    limma_all <- c(limma_all, de_tmp$limma)
  }

  edger_all <- unique(edger_all) |> sort()
  limma_all <- unique(limma_all) |> sort()

  res <- run_analysis_for_gene_sets(
    tox_group_uniprot = tox_group_uniprot,
    edger_uniprot = edger_all,
    limma_uniprot = limma_all,
    result_dir = result_dir,
    cgc_file = cgc_file,
    ncg_file = ncg_file,
    oncokb_file = oncokb_file
  )

  fisher_results <- res$fisher_results
  fisher_results$cancer_id <- "ALL_CANCERS_COMBINED"

  write_tsv(
    fisher_results %>% select(cancer_id, everything()),
    file.path(result_dir, "fisher_results_with_cancer_id.tsv")
  )

  invisible(fisher_results %>% select(cancer_id, everything()))
}

create_diagnostics_tables <- function(results_df, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Ensure expected columns exist
  required_cols <- c(
    "cancer_id", "group", "cancer_label", "detected_n",
    "cancer_detected_n", "noncancer_detected_n",
    "precision", "odds_ratio", "p_value", "p_adjust_BH"
  )

  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in results_df: ",
         paste(missing_cols, collapse = ", "))
  }

  summary_by_group <- results_df %>%
    group_by(group) %>%
    summarise(
      n_tests = n(),
      n_significant_BH = sum(p_adjust_BH < 0.05, na.rm = TRUE),
      significance_percent = 100 * n_significant_BH / n_tests,

      mean_odds_ratio = mean(odds_ratio, na.rm = TRUE),
      median_odds_ratio = median(odds_ratio, na.rm = TRUE),

      mean_precision = mean(precision, na.rm = TRUE),
      median_precision = median(precision, na.rm = TRUE),
      mean_precision_percent = 100 * mean_precision,
      median_precision_percent = 100 * median_precision,

      mean_detected_n = mean(detected_n, na.rm = TRUE),
      median_detected_n = median(detected_n, na.rm = TRUE),

      mean_cancer_detected_n = mean(cancer_detected_n, na.rm = TRUE),
      median_cancer_detected_n = median(cancer_detected_n, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(desc(significance_percent), desc(mean_odds_ratio))

  write_tsv(summary_by_group,
            file.path(output_dir, "diagnostics_summary_by_group.tsv"))

  summary_by_group_label <- results_df %>%
    group_by(group, cancer_label) %>%
    summarise(
      n_tests = n(),
      n_significant_BH = sum(p_adjust_BH < 0.05, na.rm = TRUE),
      significance_percent = 100 * n_significant_BH / n_tests,

      mean_odds_ratio = mean(odds_ratio, na.rm = TRUE),
      median_odds_ratio = median(odds_ratio, na.rm = TRUE),

      mean_precision = mean(precision, na.rm = TRUE),
      mean_precision_percent = 100 * mean_precision,

      mean_detected_n = mean(detected_n, na.rm = TRUE),
      mean_cancer_detected_n = mean(cancer_detected_n, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(cancer_label, desc(significance_percent), desc(mean_odds_ratio))

  write_tsv(summary_by_group_label,
            file.path(output_dir, "diagnostics_summary_by_group_and_label.tsv"))

  summary_by_group_cancer <- results_df %>%
    group_by(cancer_id, group) %>%
    summarise(
      n_tests = n(),
      n_significant_BH = sum(p_adjust_BH < 0.05, na.rm = TRUE),
      significance_percent = 100 * n_significant_BH / n_tests,

      mean_odds_ratio = mean(odds_ratio, na.rm = TRUE),
      median_odds_ratio = median(odds_ratio, na.rm = TRUE),

      mean_precision = mean(precision, na.rm = TRUE),
      mean_precision_percent = 100 * mean_precision,

      mean_detected_n = mean(detected_n, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(cancer_id, desc(significance_percent), desc(mean_odds_ratio))

  write_tsv(summary_by_group_cancer,
            file.path(output_dir, "diagnostics_summary_by_group_and_cancer.tsv"))

  compact_overview <- results_df %>%
    group_by(group) %>%
    summarise(
      median_OR = round(median(odds_ratio, na.rm = TRUE), 3),
      sig_percent_BH = round(100 * mean(p_adjust_BH < 0.05, na.rm = TRUE), 1),
      mean_precision_percent = round(100 * mean(precision, na.rm = TRUE), 2),
      mean_detected_genes = round(mean(detected_n, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(sig_percent_BH), desc(median_OR))

  write_tsv(compact_overview,
            file.path(output_dir, "diagnostics_compact_overview.tsv"))

  invisible(list(
    summary_by_group = summary_by_group,
    summary_by_group_label = summary_by_group_label,
    summary_by_group_cancer = summary_by_group_cancer,
    compact_overview = compact_overview
  ))
}

create_hits_table <- function(results_df, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  top_groups <- c(
    "std_log_family_mean_TOX",
    "std_log_ortholog_mean_TOX",
    "full_ortholog_mean_TOX",
    "full_family_mean_TOX",
    "ALL_DGE",
    "ONLY_TOX",
    "ONLY_DGE"
  )

  hits_table <- results_df %>%
    filter(
      cancer_label == "MERGED",
      group %in% top_groups
    ) %>%
    select(
      cancer_id,
      group,
      detected_n,
      cancer_detected_n,
      noncancer_detected_n,
      precision_percent,
      odds_ratio,
      p_adjust_BH
    ) %>%
    arrange(cancer_id, match(group, top_groups))

  write_tsv(
    hits_table,
    file.path(output_dir, "hits_table_merged_selected_groups.tsv")
  )

  invisible(hits_table)
}

main <- function(cancer_ids = CANCER_IDS) {
  dir.create(OUTPUT_BASE_DIR, recursive = TRUE, showWarnings = FALSE)

  all_results <- list()

  for (cid in cancer_ids) {
    message("Running per-cancer analysis for ", cid)
    all_results[[cid]] <- run_cancer_analysis(cancer_id = cid)
  }

  message("Running combined analysis across all cancers")
  all_results[["ALL_CANCERS_COMBINED"]] <- run_combined_analysis(cancer_ids = cancer_ids)

  combined_results <- bind_rows(all_results) %>%
    select(cancer_id, everything())

  write_tsv(
    combined_results,
    file.path(OUTPUT_BASE_DIR, "all_fisher_results_combined.tsv")
  )

  diagnostics <- create_diagnostics_tables(
    results_df = combined_results,
    output_dir = file.path(OUTPUT_BASE_DIR, "diagnostics")
  )
  hits_table <- create_hits_table(
    results_df = combined_results,
    output_dir = file.path(OUTPUT_BASE_DIR, "diagnostics")
  )

  print(diagnostics$compact_overview)

  message("Done. Results written to: ", OUTPUT_BASE_DIR)
}

main()