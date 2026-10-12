# --- Integrated DESeq2 Analysis: Male Dual Comparisons ---
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

# --- 1. CONFIGURATION ---
SEX_TO_ANALYZE <- "MALES"
base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
count_file <- file.path(base_dir, "counts_MALES.csv")
meta_file <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")
parent_output <- file.path(base_dir, "MALE_Children_Results")

# --- 2. DATA LOADING & PREP ---
count_data <- read.csv(count_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)
metadata <- as.data.frame(apply(metadata, 2, trimws))

# Standardize Conditions
metadata <- metadata %>%
  rename(SampleID = run, Condition = group, SexCol = sex) %>%
  mutate(Condition = case_when(
    grepl("diabetes", tolower(Condition)) ~ "GDM",
    grepl("lean", tolower(Condition))      ~ "Lean",
    grepl("overweight", tolower(Condition)) ~ "Overweight",
    TRUE ~ "Other"
  )) %>%
  filter(tolower(SexCol) == "male")

# --- 3. COMPARISON FUNCTION ---
# This function runs a 1-vs-1 analysis and saves to the specific directory
run_dual_comparison <- function(full_counts, full_meta, grp1, grp2, target_dir) {
  
  # Create directory if missing
  if (!dir.exists(target_dir)) { dir.create(target_dir, recursive = TRUE) }
  
  # Filter data for only these two groups
  sub_meta <- full_meta %>% filter(Condition %in% c(grp1, grp2))
  sub_meta$Condition <- factor(sub_meta$Condition, levels = c(grp2, grp1)) # grp2 is reference
  
  sub_counts <- full_counts[, intersect(colnames(full_counts), sub_meta$SampleID)]
  sub_meta <- sub_meta[match(colnames(sub_counts), sub_meta$SampleID), ]
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = round(sub_counts), colData = sub_meta, design = ~ Condition)
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep,]
  dds <- DESeq(dds)
  
  res <- results(dds, contrast=c("Condition", grp1, grp2))
  write.csv(as.data.frame(res), file.path(target_dir, paste0("DESeq2_", grp1, "_vs_", grp2, ".csv")))
  
  # 1. Heatmap (Top 30)
  vsd <- vst(dds, blind=FALSE)
  top30 <- head(rownames(res[order(res$padj, na.last = TRUE), ]), 30)
  mat <- assay(vsd)[top30, ]
  mat <- mat - rowMeans(mat)
  
  pheatmap(mat, 
           annotation_col = data.frame(Condition = sub_meta$Condition, row.names = sub_meta$SampleID),
           main = paste("Top 30 DEGs:", grp1, "vs", grp2),
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           filename = file.path(target_dir, paste0("Heatmap_", grp1, "_vs_", grp2, ".pdf")),
           width = 8, height = 9)
  
  # 2. Volcano Plot
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color = (padj < 0.05 & abs(log2FoldChange) > 1)), alpha=0.5) +
    scale_color_manual(values = c("grey", "red")) +
    geom_text_repel(data=head(res_df[order(res_df$padj, na.last = TRUE), ], 10), aes(label=gene)) +
    theme_minimal() + labs(title=paste("Volcano:", grp1, "vs", grp2))
  
  ggsave(file.path(target_dir, paste0("Volcano_", grp1, "_vs_", grp2, ".pdf")), plot=volcano, width=7, height=6)
  
  cat("Completed comparison:", grp1, "vs", grp2, "\n")
}

# --- 4. EXECUTE ALL THREE DUAL COMPARISONS ---
run_dual_comparison(count_data, metadata, "GDM", "Lean", file.path(parent_output, "GDM_vs_Lean"))
run_dual_comparison(count_data, metadata, "GDM", "Overweight", file.path(parent_output, "GDM_vs_Overweight"))
run_dual_comparison(count_data, metadata, "Overweight", "Lean", file.path(parent_output, "Overweight_vs_Lean"))

cat("\n--- ALL MALE DUAL ANALYSES COMPLETE ---")
