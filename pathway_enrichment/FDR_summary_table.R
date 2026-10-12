# ============================================================
# STANDALONE SCRIPT: Nominal vs. FDR-adjusted DEG Summary Table
# Addresses Reviewer 3, Major Comment #1
# Pulls counts directly from existing DESeq2 output files — no
# new statistical analysis is performed.
# Run via: Rscript FDR_summary_table.R
# ============================================================

for (pkg in c("dplyr", "readr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(dplyr); library(readr); library(tibble) })

# ============================================================
# ---- 1. CONFIG — list every contrast file to summarize ----
# ============================================================

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
out_dir  <- file.path(base_dir, "FDR_Summary_Table_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

contrasts <- list(
  list(label = "Female: GDM vs Lean",        file = "DESeq2_All_60/Female/Female_gdm_vs_lean_Full_Results.csv"),
  list(label = "Female: GDM vs Overweight",  file = "DESeq2_All_60/Female/Female_gdm_vs_overweight_Full_Results.csv"),
  list(label = "Female: Overweight vs Lean", file = "DESeq2_All_60/Female/Female_overweight_vs_lean_Full_Results.csv"),
  list(label = "Male: GDM vs Lean",          file = "DESeq2_All_60/Male/Male_gdm_vs_lean_Full_Results.csv"),
  list(label = "Male: GDM vs Overweight",    file = "DESeq2_All_60/Male/Male_gdm_vs_overweight_Full_Results.csv"),
  list(label = "Male: Overweight vs Lean",   file = "DESeq2_All_60/Male/Male_overweight_vs_lean_Full_Results.csv"),
  list(label = "Unstratified (all 60): GDM vs Lean", file = "DESeq2_All_60/Unified_60_Sample_Results.csv")
)

# ============================================================
# ---- 2. LOOP THROUGH EACH FILE, COUNT AT EACH THRESHOLD ----
# ============================================================

summary_rows <- list()

for (i in seq_along(contrasts)) {
  entry <- contrasts[[i]]
  full_path <- file.path(base_dir, entry$file)

  if (!file.exists(full_path)) {
    cat("WARNING: file not found, skipping:", full_path, "\n")
    next
  }

  df <- read_csv(full_path, show_col_types = FALSE) %>% as.data.frame()

  # Handle the unnamed first column (Ensembl ID) consistently with earlier scripts
  if (names(df)[1] == "" | grepl("^\\.\\.\\.", names(df)[1])) {
    names(df)[1] <- "ensembl_id"
  }

  # Drop genes never detected (baseMean = 0) before counting, consistent
  # with filtering used in prior analyses this session
  df <- df %>% filter(baseMean > 0)

  n_tested        <- nrow(df)
  n_nominal_p05   <- sum(df$pvalue < 0.05, na.rm = TRUE)
  n_fdr_p05       <- sum(df$padj   < 0.05, na.rm = TRUE)
  n_fdr_p10       <- sum(df$padj   < 0.10, na.rm = TRUE)
  pct_survive_fdr <- ifelse(n_nominal_p05 > 0, round(100 * n_fdr_p05 / n_nominal_p05, 1), NA)

  summary_rows[[length(summary_rows) + 1]] <- tibble(
    Contrast                  = entry$label,
    Genes_Tested               = n_tested,
    Nominal_p_lt_0.05          = n_nominal_p05,
    FDR_padj_lt_0.05           = n_fdr_p05,
    FDR_padj_lt_0.10           = n_fdr_p10,
    Pct_Nominal_Surviving_FDR  = pct_survive_fdr
  )

  cat("Processed:", entry$label, "-",
      n_nominal_p05, "nominal-significant,",
      n_fdr_p05, "FDR-significant (padj<0.05)\n")
}

summary_table <- bind_rows(summary_rows)

cat("\n=== FDR Summary Table ===\n")
print(summary_table)

write_csv(summary_table, file.path(out_dir, "FDR_Summary_Table.csv"))
cat("\nSaved to:", file.path(out_dir, "FDR_Summary_Table.csv"), "\n")
