# --- Integrated DESeq2 Analysis & Visualization ---
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel) # For better label placement in Volcano plot

# --- 1. CONFIGURATION ---
# SET THIS MANUALLY: "MALES" or "FEMALES"
SEX_TO_ANALYZE <- "FEMALES" 

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
count_file <- file.path(base_dir, paste0("counts_", SEX_TO_ANALYZE, ".csv"))
meta_file <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")

# Custom output directory as requested
output_dir <- paste0(ifelse(SEX_TO_ANALYZE == "MALES", "MALE", "FEMALE"), "_Children_Results")
if (!dir.exists(output_dir)) { dir.create(output_dir) }

# --- 2. DATA LOADING & PREP ---
count_data <- read.csv(count_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_file)

# [cite_start]Standardize and filter for Lean vs GDM only [cite: 1, 4]
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

# Factorize and set reference
metadata$Condition <- factor(metadata$Condition, levels = c("Lean", "GDM"))

# --- 3. RUN DESEQ2 ---
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
  ggtitle(paste("PCA Plot -", SEX_TO_ANALYZE))
ggsave(file.path(output_dir, "PCA_Plot.pdf"), pca_plot)

# B. Volcano Plot
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj < 0.05), alpha=0.5) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(title=paste("GDM vs Lean -", SEX_TO_ANALYZE), x="Log2 Fold Change", y="-log10 Adjusted P-value")
ggsave(file.path(output_dir, "Volcano_Plot.pdf"), volcano)

# C. Heatmap of Top 30 Genes
top_genes <- head(order(res$padj), 30)
mat <- assay(vst_data)[top_genes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=data.frame(Condition=metadata$Condition, row.names=metadata$SampleID),
         main=paste("Top 30 DEGs -", SEX_TO_ANALYZE),
         filename=file.path(output_dir, "Heatmap_Top30.pdf"))

# Save CSV results
write.csv(as.data.frame(res[order(res$padj),]), 
          file = file.path(output_dir, "DGE_Results.csv"))

cat("\nFinished! Results and Figures saved in:", output_dir, "\n")
