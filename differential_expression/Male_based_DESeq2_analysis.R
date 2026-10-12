# --- Integrated DESeq2 Analysis & Visualization ---
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)

# --- 1. CONFIGURATION ---
# Switching to MALES
SEX_TO_ANALYZE <- "MALES"

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
count_file <- file.path(base_dir, "counts_MALES.csv")
meta_file <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")

# Ensure output directory matches your requested folder
output_dir <- "/shared/rc/malariarna/ECFC/gene_count/MALE_Children_Results"
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

# --- 2. DATA LOADING & PREP ---
count_data <- read.csv(count_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_file)

# Standardize and filter for Lean vs GDM only 
# (This ensures "overweight" or "obese" are excluded for a clean G vs L baseline)
metadata <- metadata %>%
  rename(SampleID = run, Condition = group) %>%
  mutate(Condition = case_when(
    Condition == "gestational diabetes" ~ "GDM",
    Condition == "lean" ~ "Lean",
    TRUE ~ "Other"
  )) %>%
  filter(Condition %in% c("GDM", "Lean"))

# Align counts with filtered metadata
common_samples <- intersect(colnames(count_data), metadata$SampleID)
count_data <- count_data[, common_samples]
metadata <- metadata[match(colnames(count_data), metadata$SampleID), ]

# Factorize and set reference to "Lean"
metadata$Condition <- factor(metadata$Condition, levels = c("Lean", "GDM"))

# --- 3. RUN DESEQ2 ---
cat("Running DESeq2 for Male Cohort...\n")
dds <- DESeqDataSetFromMatrix(countData = round(count_data), colData = metadata, design = ~ Condition)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("Condition", "GDM", "Lean"))

# --- 4. GENERATE FIGURES ---

# A. PCA Plot (Quality Control)
vst_data <- vst(dds, blind=FALSE)
pca_plot <- plotPCA(vst_data, intgroup="Condition") +
  theme_minimal() +
  ggtitle(paste("PCA Plot - Male GDM vs Lean"))
ggsave(file.path(output_dir, "PCA_Plot_Male.pdf"), pca_plot)

# B. Volcano Plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue))) + # Using pvalue for visual density
  geom_point(aes(color=padj < 0.05), alpha=0.5) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(title=paste("Male GDM vs Lean"), x="Log2 Fold Change", y="-log10 P-value")
ggsave(file.path(output_dir, "Volcano_Plot_Male.pdf"), volcano)

# C. Heatmap of Top 30 Genes
top_genes <- head(order(res$padj), 30)
mat <- assay(vst_data)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=data.frame(Condition=metadata$Condition, row.names=metadata$SampleID),
         main=paste("Top 30 DEGs - Male GDM vs Lean"),
         filename=file.path(output_dir, "Heatmap_Top30_Male.pdf"))

# Save CSV results for GSEA
write.csv(as.data.frame(res[order(res$padj),]),
          file = file.path(output_dir, "DGE_Results_Male_G_vs_L.csv"))

cat("\nFinished! Results saved in:", output_dir, "\n")
