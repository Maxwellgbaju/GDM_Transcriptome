# 1. LOAD ALL NECESSARY LIBRARIES
library(DESeq2)
library(pheatmap)
library(dplyr)
library(org.Hs.eg.db)
library(gridExtra)
library(ggplotify)
library(ggplot2) # Added this to fix the ggsave error

# 2. SETUP & METADATA
metadata <- read.csv("METADATA_For_Sample_Sex.csv", strip.white = TRUE, stringsAsFactors = FALSE)
metadata$group <- tolower(trimws(metadata$group))
metadata$group <- gsub("gestational diabetes", "gdm", metadata$group)
metadata$group <- gsub("overweightfemale", "overweight", metadata$group)
metadata$run <- trimws(metadata$run)

# 3. UPDATED FUNCTION
get_heatmap_plot <- function(count_file, target_sex, comp_index) {
  
  counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
  sample_metadata <- metadata %>% 
    filter(tolower(sex) == tolower(target_sex)) %>%
    filter(run %in% colnames(counts))
  
  counts_filtered <- counts[, sample_metadata$run]
  colData <- data.frame(row.names = sample_metadata$run,
                        Condition = factor(sample_metadata$group, levels = c("lean", "overweight", "gdm")))
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts_filtered), colData = colData, design = ~ Condition)
  dds <- DESeq(dds)
  
  comparisons <- list(c("Condition", "overweight", "lean"), 
                      c("Condition", "gdm", "lean"),
                      c("Condition", "gdm", "overweight"))
  
  comp <- comparisons[[comp_index]]
  res <- results(dds, contrast = comp)
  res_df <- as.data.frame(res)
  
  # Clean IDs and Map Symbols (Fixed the 1:many warning logic)
  clean_ids <- gsub("\\..*", "", rownames(res_df))
  gene_map <- mapIds(org.Hs.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  res_df$symbol <- gene_map
  res_df$display_name <- ifelse(is.na(res_df$symbol), rownames(res_df), res_df$symbol)
  
  # Filtering Top 30
  sig_res <- res_df %>% filter(!is.na(pvalue)) %>% filter(pvalue < 0.05 & abs(log2FoldChange) >= 0.58)
  top_up <- sig_res %>% filter(log2FoldChange > 0) %>% arrange(pvalue) %>% head(15)
  top_down <- sig_res %>% filter(log2FoldChange < 0) %>% arrange(pvalue) %>% head(15)
  top_genes <- rbind(top_up, top_down)
  
  vsd <- vst(dds, blind=FALSE)
  mat <- assay(vsd)[rownames(top_genes), ]
  rownames(mat) <- top_genes$display_name
  mat <- mat - rowMeans(mat)
  
  p <- pheatmap(mat, 
                main = paste(target_sex, ":", comp[2], "vs", comp[3]),
                annotation_col = colData,
                show_colnames = FALSE,
                legend = TRUE, 
                color = colorRampPalette(c("orange", "white", "red"))(100),
                silent = TRUE)
  
  return(as.ggplot(p))
}

# 4. EXECUTION
message("Generating 6 panels for the multipanel figure...")
f1 <- get_heatmap_plot("counts_FEMALES.csv", "Female", 1)
f2 <- get_heatmap_plot("counts_FEMALES.csv", "Female", 2)
f3 <- get_heatmap_plot("counts_FEMALES.csv", "Female", 3)

m1 <- get_heatmap_plot("counts_MALES.csv", "Male", 1)
m2 <- get_heatmap_plot("counts_MALES.csv", "Male", 2)
m3 <- get_heatmap_plot("counts_MALES.csv", "Male", 3)

# 5. STITCHING
final_fig <- grid.arrange(f1, f2, f3, m1, m2, m3, ncol = 3, nrow = 2)

# 6. SAVING (Now with ggplot2 loaded, this will work)
ggsave("Figure_1_Combined_Transcriptome_350dpi.png", 
       plot = final_fig, 
       width = 22, height = 14, units = "in", dpi = 350)

message("SUCCESS! Figure_1_Combined_Transcriptome_350dpi.png is ready for your thesis.")
