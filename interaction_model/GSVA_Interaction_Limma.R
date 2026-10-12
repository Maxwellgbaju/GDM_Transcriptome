if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma", update = FALSE, ask = FALSE)
for (pkg in c("dplyr", "readr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(limma); library(dplyr); library(readr); library(tibble) })

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
out_dir  <- file.path(base_dir, "Interaction_Model_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

path_gsva     <- file.path(base_dir, "ALL_COUNTS", "GSVA_Pathway_Scores.csv")
path_metadata <- file.path(base_dir, "METADATA_For_Sample_Sex_corrected.csv")

sample_id_col_meta <- "run"
sex_col_meta       <- "sex"
condition_col_meta <- "group"
level_lean   <- "lean"
level_over   <- "overweight"
level_gdm    <- "gestational diabetes"
level_male   <- "male"
level_female <- "female"

gsva_raw <- read_csv(path_gsva, show_col_types = FALSE) %>% as.data.frame()
pathway_col1_name <- names(gsva_raw)[1]
rownames(gsva_raw) <- gsva_raw[[pathway_col1_name]]
gsva_mat <- as.matrix(gsva_raw[, -1])

metadata <- read_csv(path_metadata, show_col_types = FALSE) %>% as.data.frame()
metadata[[sample_id_col_meta]] <- trimws(metadata[[sample_id_col_meta]])
rownames(metadata) <- metadata[[sample_id_col_meta]]
metadata <- metadata[colnames(gsva_mat), , drop = FALSE]
metadata[[sex_col_meta]]       <- tolower(trimws(metadata[[sex_col_meta]]))
metadata[[condition_col_meta]] <- tolower(trimws(metadata[[condition_col_meta]]))
metadata[[condition_col_meta]] <- factor(metadata[[condition_col_meta]], levels = c(level_lean, level_over, level_gdm))
metadata[[sex_col_meta]]       <- factor(metadata[[sex_col_meta]], levels = c(level_male, level_female))

design <- model.matrix(as.formula(paste("~", sex_col_meta, "*", condition_col_meta)), data = metadata)
cat("Design matrix columns:\n"); print(colnames(design))

fit <- lmFit(gsva_mat, design)
fit <- eBayes(fit)

interaction_cols <- grep(":", colnames(design), value = TRUE)
cat("\nInteraction coefficients found:\n"); print(interaction_cols)

for (coef_name in interaction_cols) {
  cat("\n\n=== GSVA interaction coefficient:", coef_name, "===\n")
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>% rownames_to_column("Pathway_ID")
  cat("Significant pathways (adj.P.Val<0.05):", sum(tt$adj.P.Val < 0.05), "\n")
  cat("\nE2F_TARGETS and G2M_CHECKPOINT rows:\n")
  print(tt %>% filter(grepl("E2F_TARGETS|G2M_CHECKPOINT", Pathway_ID)) %>%
          select(Pathway_ID, logFC, t, P.Value, adj.P.Val))
  cat("\nTop 10 pathways overall:\n")
  print(head(tt %>% select(Pathway_ID, logFC, t, P.Value, adj.P.Val), 10))
  safe_name <- gsub("[^A-Za-z0-9]", "_", coef_name)
  write_csv(tt, file.path(out_dir, paste0("GSVA_Interaction_", safe_name, ".csv")))
}
cat("\nAll outputs saved to:", out_dir, "\n")
