# 1. SETUP
library(tidyverse)

# Define paths
target_dir <- "/shared/rc/malariarna/ECFC/gene_count/FEMALE_Children_Results/GDM_vs_Overweight"
input_rnk  <- file.path(target_dir, "GDM_vs_Overweight.rnk")
output_rnk <- file.path(target_dir, "GDM_vs_Overweight_FINAL.rnk")

# 2. LOAD AND CLEAN
if (file.exists(input_rnk)) {
    cat("Loading raw .rnk file...\n")
    
    # Read the file (assuming it was saved without a header per your previous script)
    data <- read.table(input_rnk, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    
    # Ensure there are no NAs or Infinite values which web tools hate
    data_clean <- data %>%
        filter(!is.na(V1) & !is.na(V2)) %>%       # Remove NAs
        filter(!is.infinite(V2)) %>%              # Remove Inf
        distinct(V1, .keep_all = TRUE)            # Ensure unique Gene Symbols

    # 3. SAVE FOR WEB TOOLS
    # WebGestalt and GSEA Desktop require NO header and tab-separated values
    write.table(data_clean, 
                file = output_rnk, 
                sep = "\t", 
                row.names = FALSE, 
                col.names = FALSE, 
                quote = FALSE)
    
    cat("SUCCESS: Cleaned file saved for Web Analysis at:\n", output_rnk, "\n")
} else {
    cat("ERROR: .rnk file not found. Check your directory path.\n")
}
