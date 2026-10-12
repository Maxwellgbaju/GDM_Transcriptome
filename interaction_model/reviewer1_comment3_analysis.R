# ============================================================
# Reviewer 1 Comment #3 — Sex-Dimorphism Follow-up Analyses
# (1) Top-12 heatmap with forced group clustering
# (2) Male vs Female DEG correlation
# (3) Sex-covariate DESeq2 model (~ sex + condition)
# Run on Linux HPC via: Rscript reviewer1_comment3_analysis.R
# ============================================================

# ---- 0. PACKAGES ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap", repos = "https://cloud.r-project.org")
for (pkg in c("dplyr", "ggplot2", "tibble", "readr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

suppressMessages({
  library(DESeq2); library(pheatmap); library(dplyr)
  library(ggplot2); library(tibble); library(readr)
})

# ============================================================
# ---- 1. CONFIG (confirmed values as of your last two checks) ----
# ============================================================

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"

path_counts     <- file.path(base_dir, "ALL_COUNTS", "Merged_Counts_Master_Final.csv")
path_metadata   <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")
path_male_res   <- file.path(base_dir, "DESeq2_All_60", "Male",   "Male_gdm_vs_lean_Full_Results.csv")
path_female_res <- file.path(base_dir, "DESeq2_All_60", "Female", "Female_gdm_vs_lean_Full_Results.csv")

out_dir <- file.path(base_dir, "Reviewer1_Comment3_Outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_id_col_meta <- "run"      # confirmed: matches SRR... headers in counts matrix
sex_col_meta       <- "sex"      # confirmed
condition_col_meta <- "group"    # confirmed: lowercase "group"

level_gdm    <- "gestational diabetes"  # confirmed from your metadata paste
level_lean   <- "lean"                  # NOT YET CONFIRMED — script will print actual values and stop if this is wrong
level_male   <- "male"                  # confirmed
level_female <- "female"                # confirmed

# ============================================================
# ---- 2. LOAD DATA ----
# ============================================================

male_res   <- read_csv(path_male_res, show_col_types = FALSE)
female_res <- read_csv(path_female_res, show_col_types = FALSE)

names(male_res)[1]   <- "ensembl_id"
names(female_res)[1] <- "ensembl_id"
male_res$ensembl_id   <- sub("\\..*$", "", male_res$ensembl_id)
female_res$ensembl_id <- sub("\\..*$", "", female_res$ensembl_id)

male_res   <- male_res   %>% filter(baseMean > 0)
female_res <- female_res %>% filter(baseMean > 0)

counts_raw <- read_csv(path_counts, show_col_types = FALSE) %>% as.data.frame()
names(counts_raw)[1] <- "ensembl_id"   # renames "Geneid" to a consistent name
counts_raw$ensembl_id <- sub("\\..*$", "", counts_raw$ensembl_id)
rownames(counts_raw) <- counts_raw$ensembl_id
counts_raw$ensembl_id <- NULL

metadata <- read_csv(path_metadata, show_col_types = FALSE) %>% as.data.frame()
metadata[[sample_id_col_meta]] <- trimws(metadata[[sample_id_col_meta]])  # strip stray whitespace seen in some rows
rownames(metadata) <- metadata[[sample_id_col_meta]]

stopifnot(all(colnames(counts_raw) %in% rownames(metadata)))
metadata <- metadata[colnames(counts_raw), , drop = FALSE]
metadata[[sex_col_meta]]       <- tolower(trimws(metadata[[sex_col_meta]]))
metadata[[condition_col_meta]] <- tolower(trimws(metadata[[condition_col_meta]]))

actual_conditions <- unique(metadata[[condition_col_meta]])
actual_sexes       <- unique(metadata[[sex_col_meta]])

cat("=== Sanity check ===\n")
cat("Sex levels found:      ", paste(actual_sexes, collapse = ", "), "\n")
cat("Condition levels found:", paste(actual_conditions, collapse = ", "), "\n\n")

# SAFETY CHECK — stop here rather than silently running with wrong labels
if (!(level_gdm %in% actual_conditions) || !(level_lean %in% actual_conditions)) {
  stop(paste0(
    "\n\nSTOPPING: level_gdm ('", level_gdm, "') or level_lean ('", level_lean,
    "') not found in the actual condition values above.\n",
    "Update level_lean (and level_gdm if needed) in the CONFIG section to match exactly, then rerun.\n"
  ))
}
cat("Level check passed — proceeding with analysis.\n\n")

# ============================================================
# ---- 3. PART A — Male vs Female DEG correlation ----
# ============================================================

merged_mf <- inner_join(
  male_res   %>% select(ensembl_id, symbol, log2FoldChange_male = log2FoldChange, pvalue_male = pvalue),
  female_res %>% select(ensembl_id, log2FoldChange_female = log2FoldChange, pvalue_female = pvalue),
  by = "ensembl_id"
)

cor_pearson  <- cor(merged_mf$log2FoldChange_male, merged_mf$log2FoldChange_female, method = "pearson", use = "complete.obs")
cor_spearman <- cor(merged_mf$log2FoldChange_male, merged_mf$log2FoldChange_female, method = "spearman", use = "complete.obs")

cat("=== Male vs Female log2FC correlation (GDM vs Lean) ===\n")
cat("Pearson r  =", round(cor_pearson, 3), "\n")
cat("Spearman r =", round(cor_spearman, 3), "\n")
cat("n genes compared =", nrow(merged_mf), "\n\n")

p_corr <- ggplot(merged_mf, aes(x = log2FoldChange_male, y = log2FoldChange_female)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Male vs Female DEG log2FC (GDM vs Lean)",
       subtitle = paste0("Pearson r = ", round(cor_pearson, 3), " | Spearman r = ", round(cor_spearman, 3)),
       x = "log2FC (Male, GDM vs Lean)", y = "log2FC (Female, GDM vs Lean)") +
  theme_bw(base_size = 13)

ggsave(file.path(out_dir, "Male_vs_Female_DEG_correlation.png"), p_corr, width = 7, height = 6, dpi = 300)
write_csv(merged_mf, file.path(out_dir, "Male_Female_merged_log2FC.csv"))

# ============================================================
# ---- 4. PART B — Top-12 heatmap, forced group clustering ----
# ============================================================

keep_samples <- rownames(metadata)[metadata[[condition_col_meta]] %in% c(level_lean, level_gdm)]
counts_sub   <- counts_raw[, keep_samples]
meta_sub     <- metadata[keep_samples, , drop = FALSE]

dds_norm <- DESeqDataSetFromMatrix(countData = round(counts_sub), colData = meta_sub, design = ~1)
dds_norm <- estimateSizeFactors(dds_norm)
vst_mat  <- assay(vst(dds_norm, blind = TRUE))

top12_df <- female_res %>% filter(!is.na(padj)) %>% arrange(padj) %>% slice_head(n = 12)
top12_ids <- top12_df$ensembl_id
symbol_lookup <- setNames(top12_df$symbol, top12_df$ensembl_id)

heat_mat <- vst_mat[rownames(vst_mat) %in% top12_ids, , drop = FALSE]
rownames(heat_mat) <- ifelse(is.na(symbol_lookup[rownames(heat_mat)]),
                              rownames(heat_mat), symbol_lookup[rownames(heat_mat)])
heat_mat <- t(scale(t(heat_mat)))

meta_sub$group_label <- paste0(
  ifelse(meta_sub[[sex_col_meta]] == level_female, "Female", "Male"), "-",
  ifelse(meta_sub[[condition_col_meta]] == level_gdm, "GDM", "Lean")
)

group_order <- c("Female-Lean", "Female-GDM", "Male-Lean", "Male-GDM")
ordered_samples <- unlist(lapply(group_order, function(g) rownames(meta_sub)[meta_sub$group_label == g]))
heat_mat_ordered <- heat_mat[, ordered_samples]

annotation_col <- meta_sub[ordered_samples, "group_label", drop = FALSE]
colnames(annotation_col) <- "Group"

png(file.path(out_dir, "Top12_Heatmap_ForcedGrouping.png"), width = 2400, height = 1600, res = 300)
pheatmap(heat_mat_ordered, cluster_cols = FALSE, cluster_rows = TRUE,
         annotation_col = annotation_col, show_colnames = FALSE,
         main = "Top 12 Female GDM-vs-Lean DEGs — Forced Group Clustering", fontsize_row = 10)
dev.off()

# ============================================================
# ---- 5. PART C — Sex-covariate DESeq2 model (~ sex + condition) ----
# ============================================================

meta_sub[[condition_col_meta]] <- factor(meta_sub[[condition_col_meta]], levels = c(level_lean, level_gdm))
meta_sub[[sex_col_meta]]       <- factor(meta_sub[[sex_col_meta]], levels = c(level_male, level_female))

dds_cov <- DESeqDataSetFromMatrix(countData = round(counts_sub), colData = meta_sub,
                                   design = as.formula(paste("~", sex_col_meta, "+", condition_col_meta)))
dds_cov <- DESeq(dds_cov)
res_cov <- results(dds_cov, contrast = c(condition_col_meta, level_gdm, level_lean)) %>%
  as.data.frame() %>% rownames_to_column("ensembl_id") %>%
  left_join(female_res %>% select(ensembl_id, symbol), by = "ensembl_id")

write_csv(res_cov, file.path(out_dir, "SexCovariate_GDMvsLean_Results.csv"))

cat("\n=== Sex-covariate model (~ sex + condition) summary ===\n")
cat("Significant genes (padj < 0.05):", sum(res_cov$padj < 0.05, na.rm = TRUE), "\n")
cat("Significant genes (nominal p < 0.05):", sum(res_cov$pvalue < 0.05, na.rm = TRUE), "\n\n")

compare_cov <- res_cov %>% select(ensembl_id, log2FoldChange_cov = log2FoldChange) %>%
  inner_join(female_res %>% select(ensembl_id, log2FoldChange_female = log2FoldChange), by = "ensembl_id") %>%
  inner_join(male_res   %>% select(ensembl_id, log2FoldChange_male   = log2FoldChange), by = "ensembl_id")

cor_cov_female <- cor(compare_cov$log2FoldChange_cov, compare_cov$log2FoldChange_female, use = "complete.obs")
cor_cov_male   <- cor(compare_cov$log2FoldChange_cov, compare_cov$log2FoldChange_male,   use = "complete.obs")

cat("Correlation: covariate-model log2FC vs FEMALE-only log2FC =", round(cor_cov_female, 3), "\n")
cat("Correlation: covariate-model log2FC vs MALE-only log2FC   =", round(cor_cov_male, 3), "\n\n")

genes_of_interest <- c("E2F1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC6", "CDC25A")
res_cov_goi <- res_cov %>% filter(symbol %in% genes_of_interest)
cat("=== E2F1/cell-cycle gene status in sex-covariate model ===\n")
print(res_cov_goi)

write_csv(compare_cov, file.path(out_dir, "Covariate_vs_SexStratified_log2FC.csv"))
cat("\nAll outputs saved to:", out_dir, "\n")
