## ------------------------------------------------------------
## Orthogroup / Family mapping for Homo sapiens–only analyses
## ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)

## -------------------------
## 1. Read orthogroup table
## -------------------------

ortho <- read_tsv(
  "/media/BioNAS2/bachelor_thesis_aaron_schroeder/material/orthologs.tsv",
  col_types = cols()
)

## Expected columns:
## Orthogroup | Protein1 | Species1 | Protein2 | Species2


## ----------------------------------------
## 2. Convert to long format (proteins)
## ----------------------------------------

ortho_long <- ortho %>%
  transmute(
    Orthogroup,
    Protein_ID_1 = sapply(strsplit(Protein1, "\\|"), `[`, 2),
    Species_1    = Species1,
    Protein_ID_2 = sapply(strsplit(Protein2, "\\|"), `[`, 2),
    Species_2    = Species2
  ) %>%
  pivot_longer(
    cols = starts_with("Protein_ID"),
    names_to = "source",
    values_to = "Protein_ID"
  ) %>%
  mutate(
    Species = ifelse(source == "Protein_ID_1", Species_1, Species_2)
  ) %>%
  select(Orthogroup, Protein_ID, Species) %>%
  distinct()

## Result:
## Orthogroup | Protein_ID | Species


## -------------------------------------------------
## 3. Identify orthogroups with non-human orthologs
## -------------------------------------------------

orthogroup_flags <- ortho_long %>%
  group_by(Orthogroup) %>%
  summarise(
    has_non_hsap = any(Species != "HUMAN"),
    .groups = "drop"
  )

## ---------------------------------------------------
## 4. Final mapping table (HUMAN-focused, annotated)
## ---------------------------------------------------

family_mapping <- ortho_long %>%
  left_join(orthogroup_flags, by = "Orthogroup")

## ---------------------------------------------------
## 5. Optional: restrict to Homo sapiens only
## ---------------------------------------------------

family_mapping_hsap <- family_mapping %>%
  filter(Species == "HUMAN")

## ---------------------------------------------------
## 6. Save mapping for downstream analyses
## ---------------------------------------------------

write_tsv(
  family_mapping_hsap,
  "/media/BioNAS2/bachelor_thesis_aaron_schroeder/results/ortholog_proteins.csv"
)
