library(dplyr); library(readr)
out_dir <- "/shared/rc/malariarna/ECFC/gene_count/Interaction_Model_Output"
genes_of_interest <- c("E2F1", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CDC6", "CDC25A", "CHEK1")

gdm_file <- list.files(out_dir, pattern = "gestational", full.names = TRUE)[1]
res <- read_csv(gdm_file, show_col_types = FALSE) %>% arrange(pvalue) %>%
  mutate(rank = row_number())

cat("=== E2F1-regulon genes: rank and stats in the GDM interaction test ===\n")
print(res %>% filter(symbol %in% genes_of_interest) %>%
        select(symbol, rank, log2FoldChange, stat, pvalue, padj))
cat("\n(Total genes tested:", nrow(res), ")\n")
