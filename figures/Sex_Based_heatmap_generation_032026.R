# Load required libraries
library(DESeq2)
library(pheatmap)
library(dplyr)
library(org.Hs.eg.db) # This is the dictionary you just installed

# 1. SETUP
main_dir <- "DESeq2_All_60"
if (!dir.exists(main_dir)) dir.create(main_dir)
metadata <- read.csv("METADATA_For_Sample_Sex.csv", strip.white = TRUE, stringsAsFactors = FALSE)
metadata$group <- tolower(trimws(metadata$group))
metadata$group <- gsub("gestational diabetes", "gdm", metadata$group)
metadata$group <- gsub("overweightfemale", "overweight", metadata$group) 
metadata$run <- trimws(metadata$run)

# 2. ANALYSIS FUNCTION
run_final_stratified <- function(count_file, target_sex) {
  
  message(paste("\n--- Processing:", target_sex, "---"))
  sex_dir <- file.path(main_dir, target_sex)
  if (!dir.exists(sex_dir)) dir.create(sex_dir)
  
  counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
  colnames(counts) <- trimws(colnames(counts))
  sample_metadata <- metadata %>% 
    filter(tolower(sex) == tolower(target_sex)) %>%
    filter(run %in% colnames(counts))
  
  counts_filtered <- counts[, sample_metadata$run]
  colData <- data.frame(
    row.names = sample_metadata$run,
    Condition = factor(sample_metadata$group, levels = c("lean", "overweight", "gdm"))
  )
  
  dds <- DESeqDataSetFromMatrix(countData = round(counts_filtered), 
                                colData = colData, 
                                design = ~ Condition)
  dds <- DESeq(dds)
  
  comparisons <- list(c("Condition", "gdm", "lean"), 
                      c("Condition", "gdm", "overweight"), 
                      c("Condition", "overweight", "lean"))
  
  for (comp in comparisons) {
    comp_label <- paste(target_sex, comp[2], "vs", comp[3], sep="_")
    res <- results(dds, contrast = comp)
    res_df <- as.data.frame(res)
    
    # --- MAPPING STEP ---
    # 1. Strip the version numbers (e.g. .18)
    clean_ids <- gsub("\\..*", "", rownames(res_df))
    
    # 2. Look up symbols
    res_df$symbol <- mapIds(org.Hs.eg.db,
                           keys = clean_ids,
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
    
    # 3. Create a display name (Symbol if found, otherwise keep ID)
    res_df$display_name <- ifelse(is.na(res_df$symbol), rownames(res_df), res_df$symbol)
    
    write.csv(res_df, file.path(sex_dir, paste0(comp_label, "_Full_Results.csv")))
    
    # FILTER: Nominal p < 0.05 and Log2FC > 0.58
    sig_res <- res_df %>% 
      filter(!is.na(pvalue)) %>% 
      filter(pvalue < 0.05 & abs(log2FoldChange) >= 0.58)
    
    if (nrow(sig_res) >= 2) {
      top_up <- sig_res %>% filter(log2FoldChange > 0) %>% arrange(pvalue) %>% head(15)
      top_down <- sig_res %>% filter(log2FoldChange < 0) %>% arrange(pvalue) %>% head(15)
      top_genes <- rbind(top_up, top_down)
      
      vsd <- vst(dds, blind=FALSE)
      mat <- assay(vsd)[rownames(top_genes), ]
      
      # RENAME ROWS TO SYMBOLS
      rownames(mat) <- top_genes$display_name
      
      mat <- mat - rowMeans(mat) 
      
      # SAVE 350 DPI HEATMAP
      png(filename = file.path(sex_dir, paste0(comp_label, "_Heatmap_350dpi.png")), 
          width = 3200, height = 2800, res = 350)
      
      pheatmap(mat, 
               main = paste("Top 30 DEGs (p < 0.05):", comp_label),
               annotation_col = colData,
               show_colnames = FALSE,
               color = colorRampPalette(c("orange", "white", "red"))(100),
               border_color = NA,
               fontsize_row = 10) 
      dev.off()
      message(paste("Successfully saved Heatmap with Symbols for:", comp_label))
    }
  }
}

# 3. EXECUTE
run_final_stratified("counts_FEMALES.csv", "Female")
run_final_stratified("counts_MALES.csv", "Male")
