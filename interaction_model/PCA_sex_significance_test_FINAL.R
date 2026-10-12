# ============================================================
# CORRECTED: Reproduces Fig. 2 exactly (DESeq2::plotPCA, top 500 genes,
# design ~ Sex + Condition, blind = FALSE), then tests PC1 by sex.
# Run via: Rscript PCA_sex_significance_test_v2.R
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
for (pkg in c("dplyr", "vegan")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(DESeq2); library(dplyr); library(vegan) })

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
out_dir <- file.path(base_dir, "PCA_Sex_Test_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load exactly as the original script does ----
counts_f <- read.csv(file.path(base_dir, "counts_FEMALES.csv"), check.names = FALSE)
counts_m <- read.csv(file.path(base_dir, "counts_MALES.csv"), check.names = FALSE)

count_matrix <- full_join(counts_f, counts_m, by = "GeneID")
count_matrix[is.na(count_matrix)] <- 0
rownames(count_matrix) <- count_matrix$GeneID
count_matrix <- count_matrix[, -1]

sample_ids <- colnames(count_matrix)
sex_labels <- ifelse(sample_ids %in% colnames(counts_f), "Female", "Male")
srr_numbers <- as.numeric(gsub("SRR2407", "", sample_ids))
condition_labels <- case_when(
  srr_numbers >= 1876 & srr_numbers <= 1895 ~ "Lean",
  srr_numbers >= 1896 & srr_numbers <= 1915 ~ "Overweight",
  srr_numbers >= 1916 & srr_numbers <= 1935 ~ "GDM"
)

metadata <- data.frame(
  row.names = sample_ids,
  Sex = sex_labels,
  Condition = factor(condition_labels, levels = c("Lean", "Overweight", "GDM"))
)

# ---- 2. Run DESeq2 exactly as the original script does ----
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix), colData = metadata, design = ~ Sex + Condition)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# ---- 3. Reproduce plotPCA with default ntop=500 (matches Fig. 2 exactly) ----
pcaData <- plotPCA(vsd, intgroup = c("Sex", "Condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
cat("PC1 variance explained:", percentVar[1], "% (should match Fig. 2's 18%)\n")
cat("PC2 variance explained:", percentVar[2], "% (should match Fig. 2's 13%)\n\n")

write.csv(pcaData, file.path(out_dir, "PC_scores_matching_Fig2.csv"), row.names = TRUE)

# ---- 4. NOW run the significance tests on the CORRECT PC1/PC2 ----
female_pc1 <- pcaData %>% filter(Sex == "Female") %>% pull(PC1)
male_pc1   <- pcaData %>% filter(Sex == "Male")   %>% pull(PC1)

t_result <- t.test(female_pc1, male_pc1)
w_result <- wilcox.test(female_pc1, male_pc1)

cat("=== PC1 by sex (matching Fig. 2's actual PCA) ===\n")
cat("Female PC1: mean =", round(mean(female_pc1), 2), ", n =", length(female_pc1), "\n")
cat("Male PC1:   mean =", round(mean(male_pc1), 2), ", n =", length(male_pc1), "\n")
cat("Welch t-test p-value:      ", format.pval(t_result$p.value, digits = 4), "\n")
cat("Wilcoxon rank-sum p-value: ", format.pval(w_result$p.value, digits = 4), "\n\n")

pc_matrix <- as.matrix(pcaData[, c("PC1", "PC2")])
dist_matrix <- dist(pc_matrix, method = "euclidean")
set.seed(42)
permanova_result <- adonis2(dist_matrix ~ pcaData$Sex, permutations = 999)

cat("=== PERMANOVA: PC1+PC2 by sex (999 permutations) ===\n")
print(permanova_result)

pcaData$predicted_side <- ifelse(pcaData$PC1 < 0, "negative", "positive")
pcaData$expected_side  <- ifelse(pcaData$Sex == "Female", "negative", "positive")
misclassified <- pcaData %>% filter(predicted_side != expected_side)
cat("\nMisclassified:", nrow(misclassified), "out of", nrow(pcaData), "\n")
print(misclassified %>% select(name, Sex, Condition, PC1, PC2))

cat("\nAll outputs saved to:", out_dir, "\n")
