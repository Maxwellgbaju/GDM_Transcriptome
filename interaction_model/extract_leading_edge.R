suppressMessages({ library(dplyr); library(readr) })
out_dir <- "/shared/rc/malariarna/ECFC/gene_count/Interaction_Model_Output"
gsea_res <- read_csv(file.path(out_dir, "GSEA_Interaction_Ranked_Results.csv"), show_col_types = FALSE)

cat("=== Full 50 Hallmark pathway results (already have this) ===\n")
cat("Total pathways:", nrow(gsea_res), "\n\n")

cat("=== Leading-edge genes: E2F_TARGETS ===\n")
e2f_row <- gsea_res %>% filter(pathway == "HALLMARK_E2F_TARGETS")
cat(e2f_row$leadingEdge, "\n\n")

cat("=== Leading-edge genes: G2M_CHECKPOINT ===\n")
g2m_row <- gsea_res %>% filter(pathway == "HALLMARK_G2M_CHECKPOINT")
cat(g2m_row$leadingEdge, "\n")
