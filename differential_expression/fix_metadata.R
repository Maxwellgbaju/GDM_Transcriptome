library(dplyr)
library(readr)

path_in  <- "/shared/rc/malariarna/ECFC/gene_count/METADATA_For_Sample_Sex.csv"
path_out <- "/shared/rc/malariarna/ECFC/gene_count/METADATA_For_Sample_Sex_corrected.csv"

meta <- read_csv(path_in, show_col_types = FALSE)

cat("=== BEFORE FIX ===\n")
print(table(meta$group, useNA = "ifany"))
cat("\nRows with the typo or missing sex:\n")
print(meta %>% filter(grepl("overweightfemale", group, ignore.case = TRUE) | is.na(sex) | sex == ""))

meta <- meta %>%
  mutate(
    sex = ifelse(grepl("overweightfemale", group, ignore.case = TRUE) & (is.na(sex) | sex == ""), "female", sex),
    group = ifelse(grepl("overweightfemale", group, ignore.case = TRUE), "overweight", group)
  )

cat("\n=== AFTER FIX ===\n")
print(table(meta$group, meta$sex, useNA = "ifany"))

write_csv(meta, path_out)
cat("\nCorrected metadata saved to:", path_out, "\n")
