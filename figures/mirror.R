# --- Load Libraries ---
library(tidyverse)
library(RRHO)

# --- 1. Load your DESeq2 results ---
res_female <- read.csv("DGE_Results_GDM_VS_Lean_Female.csv")
res_male   <- read.csv("DGE_Results_GDM_vs_Lean_Male.csv")

# Improved function to create the Rank List
prepare_rrho_list <- function(df, label) {
  gene_col <- intersect(c("Symbol", "gene_name", "Gene", "ID"), colnames(df))[1]
  if (is.na(gene_col)) gene_col <- colnames(df)[1]
  
  lfc_col <- intersect(c("log2FoldChange", "log2FC", "Log2FC"), colnames(df))[1]
  pval_col <- intersect(c("pvalue", "p.value", "PValue"), colnames(df))[1]
  
  df %>%
    filter(!is.na(!!sym(pval_col)) & !!sym(gene_col) != "") %>%
    mutate(rank_score = sign(!!sym(lfc_col)) * -log10(!!sym(pval_col))) %>%
    select(GeneID = !!sym(gene_col), rank_score) %>%
    group_by(GeneID) %>%
    summarize(rank_score = rank_score[which.max(abs(rank_score))], .groups = 'drop')
}

# Initial processing
female_raw <- prepare_rrho_list(res_female, "Female")
male_raw   <- prepare_rrho_list(res_male, "Male")

# --- 2. CRITICAL STEP: INTERSECT GENES ---
# RRHO requires the two lists to be compared gene-for-gene
common_genes <- intersect(female_raw$GeneID, male_raw$GeneID)
message(paste("Number of common genes found:", length(common_genes)))

female_list <- female_raw %>%
  filter(GeneID %in% common_genes) %>%
  arrange(desc(rank_score)) %>%
  as.data.frame()

male_list <- male_raw %>%
  filter(GeneID %in% common_genes) %>%
  arrange(desc(rank_score)) %>%
  as.data.frame()

# --- 3. Setup Output Directory ---
out_path <- "RRHO_Results_Final"
if (!dir.exists(out_path)) dir.create(out_path)

# --- 4. Run RRHO Analysis ---
# Using a larger stepsize (200) to ensure stability
RRHO_obj <- RRHO(female_list, male_list, 
                 alternative = "enrichment", 
                 stepsize = 200, 
                 labels = c("Female_GDM", "Male_GDM"),
                 outputdir = out_path,
                 plots = TRUE)

message("Success! Heatmaps are saved in the 'RRHO_Results_Final' folder.")
