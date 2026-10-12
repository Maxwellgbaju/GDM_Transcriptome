# ============================================================
# Formal sex x condition interaction test (Reviewer 3, Major #1)
# Design: ~ sex + condition + sex:condition
# Run via: Rscript Interaction_Model_DESeq2.R
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
for (pkg in c("dplyr", "readr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(DESeq2); library(dplyr); library(readr); library(tibble) })

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
out_dir  <- file.path(base_dir, "Interaction_Model_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

path_counts     <- file.path(base_dir, "ALL_COUNTS", "Merged_Counts_Master_Final.csv")
path_metadata   <- file.path(base_dir, "METADATA_For_Sample_Sex_corrected.csv")
path_female_res <- file.path(base_dir, "DESeq2_All_60", "Female", "Female_gdm_vs_lean_Full_Results.csv")

sample_id_col_meta <- "run"
sex_col_meta       <- "sex"
condition_col_meta <- "group"
level_lean   <- "lean"
level_over   <- "overweight"
level_gdm    <- "gestational diabetes"
level_male   <- "male"
level_female <- "female"

# ---- Load counts ----
counts_raw <- read_csv(path_counts, show_col_types = FALSE) %>% as.data.frame()
names(counts_raw)[1] <- "ensembl_id"
counts_raw$ensembl_id <- sub("\\..*$", "", counts_raw$ensembl_id)
rownames(counts_raw) <- counts_raw$ensembl_id
counts_raw$ensembl_id <- NULL

# ---- Load metadata ----
metadata <- read_csv(path_metadata, show_col_types = FALSE) %>% as.data.frame()
metadata[[sample_id_col_meta]] <- trimws(metadata[[sample_id_col_meta]])
rownames(metadata) <- metadata[[sample_id_col_meta]]
metadata <- metadata[colnames(counts_raw), , drop = FALSE]
metadata[[sex_col_meta]]       <- tolower(trimws(metadata[[sex_col_meta]]))
metadata[[condition_col_meta]] <- tolower(trimws(metadata[[condition_col_meta]]))
metadata[[condition_col_meta]] <- factor(metadata[[condition_col_meta]], levels = c(level_lean, level_over, level_gdm))
metadata[[sex_col_meta]]       <- factor(metadata[[sex_col_meta]], levels = c(level_male, level_female))

cat("Sample counts by group:\n")
print(table(metadata[[sex_col_meta]], metadata[[condition_col_meta]]))

# ---- Symbol lookup from existing female results ----
female_res <- read_csv(path_female_res, show_col_types = FALSE)
names(female_res)[1] <- "ensembl_id"
female_res$ensembl_id <- sub("\\..*$", "", female_res$ensembl_id)
symbol_lookup <- female_res %>% select(ensembl_id, symbol) %>% distinct()

# ---- Fit interaction model ----
form <- as.formula(paste("~", sex_col_meta, "+", condition_col_meta, "+", sex_col_meta, ":", condition_col_meta))
dds <- DESeqDataSetFromMatrix(countData = round(counts_raw), colData = metadata, design = form)
dds <- DESeq(dds)

cat("\n=== Available result coefficients ===\n")
print(resultsNames(dds))

interaction_terms <- resultsNames(dds)[grepl(paste0("^", sex_col_meta), resultsNames(dds)) &
                                        grepl(condition_col_meta, resultsNames(dds)) &
                                        grepl("sex", resultsNames(dds)) & grepl("group", resultsNames(dds))]
cat("\nDetected interaction coefficients:\n")
print(interaction_terms)

for (term in interaction_terms) {
  cat("\n\n=== Interaction term:", term, "===\n")
  res <- results(dds, name = term) %>% as.data.frame() %>% rownames_to_column("ensembl_id") %>%
    left_join(symbol_lookup, by = "ensembl_id") %>% arrange(pvalue)

  n_nominal <- sum(res$pvalue < 0.05, na.rm = TRUE)
  n_fdr     <- sum(res$padj   < 0.05, na.rm = TRUE)
  cat("Nominal p<0.05:", n_nominal, " | FDR padj<0.05:", n_fdr, " | genes tested:", sum(!is.na(res$pvalue)), "\n")

  cat("Top 15 genes by interaction p-value:\n")
  print(head(res %>% select(symbol, log2FoldChange, stat, pvalue, padj), 15))

  safe_name <- gsub("[^A-Za-z0-9]", "_", term)
  write_csv(res, file.path(out_dir, paste0("Interaction_", safe_name, "_FullResults.csv")))
}

cat("\nAll outputs saved to:", out_dir, "\n")
