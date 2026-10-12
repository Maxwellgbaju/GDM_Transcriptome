# ============================================================
# STANDALONE SCRIPT: Male GSVA Pathway Analysis (limma)
# Mirrors the existing Female_Only_Stats_Full.csv methodology exactly,
# applied to Male samples — addresses Reviewer 1 Comment #5
# ("Fig 6 GSVA shown only for female — do male GSVA too")
# Run via: Rscript Male_GSVA_limma_analysis.R
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma", update = FALSE, ask = FALSE)
for (pkg in c("dplyr", "readr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

suppressMessages({ library(limma); library(dplyr); library(readr); library(tibble) })

# ============================================================
# ---- 1. CONFIG ----
# ============================================================

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"

path_gsva     <- file.path(base_dir, "ALL_COUNTS", "GSVA_Pathway_Scores.csv")
path_metadata <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")

out_dir <- file.path(base_dir, "Male_GSVA_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_id_col_meta <- "run"
sex_col_meta       <- "sex"
condition_col_meta <- "group"
level_gdm    <- "gestational diabetes"
level_lean   <- "lean"
level_male   <- "male"

# ============================================================
# ---- 2. LOAD DATA ----
# ============================================================

gsva_raw <- read_csv(path_gsva, show_col_types = FALSE) %>% as.data.frame()
pathway_col1_name <- names(gsva_raw)[1]  # "...1"
rownames(gsva_raw) <- gsva_raw[[pathway_col1_name]]
gsva_mat <- as.matrix(gsva_raw[, -1])  # rows = pathways, cols = samples (SRR IDs)

metadata <- read_csv(path_metadata, show_col_types = FALSE) %>% as.data.frame()
metadata[[sample_id_col_meta]] <- trimws(metadata[[sample_id_col_meta]])
metadata[[sex_col_meta]]       <- tolower(trimws(metadata[[sex_col_meta]]))
metadata[[condition_col_meta]] <- tolower(trimws(metadata[[condition_col_meta]]))
rownames(metadata) <- metadata[[sample_id_col_meta]]

# ============================================================
# ---- 3. SUBSET TO MALE, LEAN+GDM SAMPLES ONLY ----
# ============================================================

male_samples <- metadata %>%
  filter(.data[[sex_col_meta]] == level_male,
         .data[[condition_col_meta]] %in% c(level_lean, level_gdm)) %>%
  pull(all_of(sample_id_col_meta))

cat("=== Male samples selected ===\n")
cat("n Lean:", sum(metadata[male_samples, condition_col_meta] == level_lean), "\n")
cat("n GDM: ", sum(metadata[male_samples, condition_col_meta] == level_gdm), "\n\n")

stopifnot(all(male_samples %in% colnames(gsva_mat)))
gsva_male <- gsva_mat[, male_samples]

meta_male <- metadata[male_samples, , drop = FALSE]
meta_male$Condition <- factor(meta_male[[condition_col_meta]], levels = c(level_lean, level_gdm))

# ============================================================
# ---- 4. LIMMA ANALYSIS (mirrors the female methodology exactly) ----
# ============================================================

design <- model.matrix(~ Condition, data = meta_male)
fit <- lmFit(gsva_male, design)
fit <- eBayes(fit)

male_results <- topTable(fit, coef = "Conditiongestational diabetes",
                          number = Inf, sort.by = "P") %>%
  rownames_to_column("Pathway_ID") %>%
  mutate(Pathway = gsub("HALLMARK_", "", gsub("_", " ", Pathway_ID)))

cat("=== Male GSVA results (GDM vs Lean), top 10 by p-value ===\n")
print(head(male_results, 10))
cat("\n")

cat("Significant pathways (nominal p<0.05):", sum(male_results$P.Value < 0.05), "\n")
cat("Significant pathways (adj.P.Val<0.05):", sum(male_results$adj.P.Val < 0.05), "\n\n")

# Explicitly flag the E2F1/cell-cycle related pathways for direct comparison to female
cell_cycle_pathways <- male_results %>%
  filter(grepl("E2F_TARGETS|G2M_CHECKPOINT|MITOTIC_SPINDLE|MYC_TARGETS", Pathway_ID))

cat("=== Cell-cycle/E2F1-related pathways in MALE GDM-vs-Lean ===\n")
print(cell_cycle_pathways %>% select(Pathway_ID, logFC, t, P.Value, adj.P.Val))
cat("\n")

write_csv(male_results, file.path(out_dir, "Male_Only_Stats_Full.csv"))
cat("Full results saved to:", file.path(out_dir, "Male_Only_Stats_Full.csv"), "\n")
