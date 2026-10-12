# 1. SETUP ENVIRONMENT
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(org.Hs.eg.db)
})

# 2. CONFIGURATION
target_dir <- "/shared/rc/malariarna/ECFC/gene_count/FEMALE_Children_Results/GDM_vs_Lean"
input_file <- file.path(target_dir, "DGE_Results_GDM_VS_Lean.csv")
gmt_file   <- "/shared/rc/malariarna/ECFC/gene_count/h.all.v2023.1.Hs.symbols.gmt"

# 3. LOAD DATA
if (!file.exists(input_file)) {
    stop(paste("ERROR: Result file not found at:", input_file))
}

cat("Step 1: Loading DGE results...\n")
res <- read.csv(input_file, row.names = 1)

# 4. ID CLEANING AND MAPPING
cat("Step 2: Cleaning IDs and mapping to Gene Symbols...\n")
clean_ensembl_ids <- gsub("\\..*$", "", rownames(res))

gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = clean_ensembl_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

res$Symbol <- gene_symbols
res_mapped <- res[!is.na(res$Symbol), ]
cat(paste("Successfully mapped", nrow(res_mapped), "genes.\n"))

# 5. CREATE RANKED LIST
cat("Step 3: Calculating Rank Scores...\n")
res_mapped$rank_score <- sign(res_mapped$log2FoldChange) * -log10(res_mapped$pvalue)

# Handle duplicates
ranks_df <- res_mapped %>%
  group_by(Symbol) %>%
  summarize(rank_score = rank_score[which.max(abs(rank_score))]) %>%
  arrange(desc(rank_score))

rank_vector <- setNames(ranks_df$rank_score, ranks_df$Symbol)

# 6. RUN GSEA
if (file.exists(gmt_file)) {
    cat("Step 4: Running GSEA with Hallmark Gene Sets...\n")
    pathways <- gmtPathways(gmt_file)
    
    fgseaRes <- fgsea(pathways = pathways, 
                      stats = rank_vector, 
                      minSize = 15, 
                      maxSize = 500)
    
    fgseaRes_ordered <- fgseaRes[order(padj), ]
    
    write.csv(as.data.frame(fgseaRes_ordered), 
              file.path(target_dir, "GSEA_Pathways_Results.csv"), 
              row.names = FALSE)
    
    # 7. GENERATE TOP PATHWAY PLOT
    if (nrow(fgseaRes_ordered) > 0) {
        top_path <- fgseaRes_ordered$pathway[1]
        cat("Step 5: Generating Enrichment Plot for top pathway:", top_path, "\n")
        
        plot_out <- plotEnrichment(pathways[[top_path]], rank_vector) + 
                    labs(title = paste("GDM vs Lean Top Pathway:", top_path),
                         subtitle = paste("NES:", round(fgseaRes_ordered$NES[1], 2), 
                                          "| adj.P:", formatC(fgseaRes_ordered$padj[1], format = "e", digits = 2)))
        
        ggsave(file.path(target_dir, "Top_Pathway_Enrichment_Plot.pdf"), 
               plot_out, width = 8, height = 6)
        
        cat("\nSUCCESS: All files saved in:", target_dir, "\n")
    } else {
        cat("WARNING: No enriched pathways found.\n")
    }
} else {
    cat("ERROR: GMT file not found. Saving .rnk instead.\n")
    write.table(ranks_df, file.path(target_dir, "GDM_vs_Lean.rnk"), 
                sep="\t", quote=F, row.names=F, col.names=F)
}
# End of script
