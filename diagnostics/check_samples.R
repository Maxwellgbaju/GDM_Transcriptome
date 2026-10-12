library(tidyverse)

# --- SETTINGS ---
count_file <- "counts_FEMALES.csv"
meta_file <- "METADATA_For_Sample_Sex.csv"

# 1. Load data
counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(meta_file, stringsAsFactors = FALSE)

# 2. Clean the metadata (strip invisible tabs/spaces)
metadata <- as.data.frame(apply(metadata, 2, trimws))

# 3. Standardize the group names
metadata <- metadata %>%
  mutate(Condition = case_when(
    grepl("diabetes", tolower(group)) ~ "GDM",
    grepl("lean", tolower(group))     ~ "Lean",
    grepl("overweight", tolower(group)) ~ "Overweight",
    TRUE ~ "Other"
  ))

# 4. Match against the count file columns
samples_in_counts <- colnames(counts)
matched_metadata <- metadata %>% filter(run %in% samples_in_counts)

# --- THE REPORT ---
cat("\n==========================================\n")
cat("      DESEQ2 SAMPLE RECOGNITION REPORT\n")
cat("==========================================\n")
cat("Total columns in count file:  ", ncol(counts), "\n")
cat("Total rows in metadata file:   ", nrow(metadata), "\n")
cat("Successfully matched samples:  ", nrow(matched_metadata), "\n\n")

cat("Breakdown of Matched Samples by Group:\n")
print(table(matched_metadata$Condition))

if(nrow(matched_metadata) < 30) {
  cat("\nWARNING: Some samples are still missing!\n")
  cat("Missing IDs from Metadata:\n")
  print(setdiff(metadata$run, samples_in_counts))
} else {
  cat("\nSUCCESS: All 30 samples are recognized and ready.\n")
}
cat("==========================================\n")
