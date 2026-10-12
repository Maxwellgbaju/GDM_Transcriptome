# --- Robust Integrated DESeq2 Analysis: All Comparisons ---
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)

# --- 1. CONFIGURATION ---
# Change this to "MALES" when you are ready to run the male cohort
SEX_TO_ANALYZE <- "FEMALES" 

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
count_file <- file.path(base_dir, paste0("counts_", SEX_TO_ANALYZE, ".csv"))
meta_file <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")

# Set main output directory
main_output <- paste0(ifelse(SEX_TO_ANALYZE == "MALES", "MALE", "FEMALE"), "_Children_Results")
if (!dir.exists(main_output)) { dir.create(main_output) }

# --- 2. LOAD & CLEAN DATA ---
count_data <- read.csv(count_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

# STAGE 1: Strip all invisible tabs and spaces from every cell (Crucial for your file)
metadata <- as.data.frame(apply(metadata, 2, trimws))

# STAGE 2: Standardize Column Names and Groups
metadata <- metadata %>%
  rename(SampleID = run, Condition = group, SexCol = sex) %>%
  mutate(Condition = case_when(
    grepl("diabetes", tolower(Condition)) ~ "GDM",
    grepl("lean", tolower(Condition))     ~ "Lean",
    grepl("overweight", tolower(Condition)) ~ "Overweight",
    TRUE ~ "Other"
  )) %>%
  # Handle potential missing sex labels for known female samples
  mutate(SexCol = ifelse(is.na(SexCol) | SexCol == "", "female", SexCol))

# STAGE 3: Filter for chosen sex and match to the Count File
target_sex <- ifelse(SEX_TO_ANALYZE == "FEMALES", "female", "male")
metadata <- metadata %>% filter(tolower(SexCol) == target_sex)

# Ensure only samples present in BOTH files are analyzed
common_samples <- intersect(colnames(count_data), metadata$SampleID)
count_data <- count_data[, common_samples]
metadata <- metadata[match(colnames(count_data), metadata$SampleID), ]

# --- DIAGNOSTIC CENSUS ---
cat("\n--- FINAL SAMPLE CENSUS FOR", SEX_TO_ANALYZE, "---\n")
print(table(metadata$Condition))
cat("Total samples successfully matched:", ncol(count_data), "\n\n")

# --- 3. RUN ALL PAIRWISE COMPARISONS ---
# This loop handles all three combinations
comparisons <- list(c("GDM", "Lean"), c("GDM", "Overweight"), c("Overweight", "Lean"))

for (pair in comparisons) {
  g1 <- pair[1]
  g2 <- pair[2]
  comparison_name <- paste0(g1, "_vs_", g2)
  cat("Processing Comparison:", comparison_name, "...\n")
  
  # Create sub-directory for this comparison
  sub_dir <- file.path(main_output, comparison_name)
  if (!dir.exists(sub_dir)) { dir.create(sub_dir) }
  
  # Subset data for the current pair
  sub_meta <- metadata %>% filter(Condition %in% c(g1, g2))
  sub_meta$Condition <- factor(sub_meta$Condition, levels = c(g2, g1)) # g2 is reference
  sub_counts <- count_data[, sub_meta$SampleID]
  
  # Run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = round(sub_counts), 
                                colData = sub_meta, 
                                design = ~ Condition)
  dds <- dds[rowSums(counts(dds)) > 1, ] # Minimal filtering to preserve samples
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Save Results CSV (Sorted by significance)
  write.csv(as.data.frame(res[order(res$padj),]), 
            file = file.path(sub_dir, "DGE_Results.csv"))
  
  # --- 4. GENERATE FIGURES ---
  vst_data <- vst(dds, blind=FALSE)
  
  # A. PCA Plot
  pca_plot <- plotPCA(vst_data, intgroup="Condition") + 
    theme_bw() + 
    ggtitle(paste("PCA:", comparison_name))
  ggsave(file.path(sub_dir, "PCA_Plot.pdf"), pca_plot)
  
  # B. Volcano Plot
  res_df <- as.data.frame(res)
  volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(aes(color=padj < 0.05), alpha=0.5) +
    theme_minimal() +
    scale_color_manual(values = c("black", "red")) +
    labs(title=paste("Volcano:", comparison_name))
  ggsave(file.path(sub_dir, "Volcano_Plot.pdf"), volcano)
  
  # C. Heatmap of Leading 30 Genes
  top30 <- head(order(res$padj, na.last = TRUE), 30)
  mat <- assay(vst_data)[top30, ]
  mat <- mat - rowMeans(mat) # Center expression
  pheatmap(mat, 
           annotation_col = data.frame(Condition=sub_meta$Condition, row.names=sub_meta$SampleID),
           main = paste("Top 30 Genes:", comparison_name),
           filename = file.path(sub_dir, "Heatmap_Top30.pdf"),
           width = 8, height = 10)
}

cat("\nAnalysis Complete! Folders created in:", main_output, "\n")
