suppressMessages({ library(dplyr); library(readr) })
base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
interaction_dir <- file.path(base_dir, "Interaction_Model_Output")

interaction_res <- read_csv(file.path(interaction_dir, "Interaction_sexfemale_groupgestational_diabetes_FullResults.csv"), show_col_types = FALSE)

cat("=== GDM interaction results file ===\n")
cat("Total rows:", nrow(interaction_res), "\n")
cat("Rows with non-NA pvalue:", sum(!is.na(interaction_res$pvalue)), "\n")
cat("Rows with baseMean > 0:", sum(interaction_res$baseMean > 0, na.rm = TRUE), "\n")
cat("Rows with baseMean > 0 AND non-NA pvalue:", sum(interaction_res$baseMean > 0 & !is.na(interaction_res$pvalue), na.rm = TRUE), "\n\n")

cat("(Manuscript reports 41,780 genes tested for this interaction)\n")
