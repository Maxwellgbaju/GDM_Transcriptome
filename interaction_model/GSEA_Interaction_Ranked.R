# ============================================================
# GSEA using genes ranked by the sex x GDM interaction statistic
# Run AFTER Interaction_Model_DESeq2.R
# Run via: Rscript GSEA_Interaction_Ranked.R
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("fgsea", quietly = TRUE))
  BiocManager::install("fgsea", update = FALSE, ask = FALSE)
for (pkg in c("dplyr", "readr", "msigdbr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(fgsea); library(dplyr); library(readr); library(msigdbr) })

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
out_dir  <- file.path(base_dir, "Interaction_Model_Output")

interaction_file <- list.files(out_dir, pattern = "gestational", full.names = TRUE)[1]
cat("Using interaction results file:", interaction_file, "\n")

res <- read_csv(interaction_file, show_col_types = FALSE) %>%
  filter(!is.na(stat), !is.na(symbol), symbol != "NA") %>%
  group_by(symbol) %>% slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>% ungroup()

ranks <- setNames(res$stat, res$symbol)
ranks <- sort(ranks, decreasing = TRUE)

hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol) %>% distinct()
pathway_list <- split(hallmark$gene_symbol, hallmark$gs_name)

set.seed(42)
gsea_res <- fgsea(pathways = pathway_list, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
gsea_res <- gsea_res %>% arrange(padj)

cat("\n=== Top 15 pathways by interaction-ranked GSEA ===\n")
print(gsea_res %>% select(pathway, NES, pval, padj, size) %>% head(15))

cat("\nSignificant pathways (padj<0.05):", sum(gsea_res$padj < 0.05, na.rm = TRUE), "\n")

cat("\n=== E2F_TARGETS and G2M_CHECKPOINT specifically ===\n")
print(gsea_res %>% filter(grepl("E2F_TARGETS|G2M_CHECKPOINT", pathway)) %>%
        select(pathway, NES, pval, padj, size))

fwrite_safe <- gsea_res
fwrite_safe$leadingEdge <- sapply(fwrite_safe$leadingEdge, paste, collapse = ";")
write_csv(fwrite_safe, file.path(out_dir, "GSEA_Interaction_Ranked_Results.csv"))

cat("\nSaved to:", file.path(out_dir, "GSEA_Interaction_Ranked_Results.csv"), "\n")
