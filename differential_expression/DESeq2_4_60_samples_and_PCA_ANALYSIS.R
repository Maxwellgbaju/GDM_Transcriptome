# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Load the two count files
counts_f <- read.csv("counts_FEMALES.csv", check.names = FALSE)
counts_m <- read.csv("counts_MALES.csv", check.names = FALSE)

# 2. Merge them into a single matrix
# We join by GeneID to ensure rows align perfectly
count_matrix <- full_join(counts_f, counts_m, by = "GeneID")
count_matrix[is.na(count_matrix)] <- 0 # Replace any missing values with 0
rownames(count_matrix) <- count_matrix$GeneID
count_matrix <- count_matrix[, -1] # Remove the GeneID column

# 3. Build the Metadata table automatically
sample_ids <- colnames(count_matrix)

# Create Sex labels based on which file they came from
sex_labels <- ifelse(sample_ids %in% colnames(counts_f), "Female", "Male")

# Create Condition labels based on SRR numerical ranges
# Convert SRR IDs to numbers for range checking
srr_numbers <- as.numeric(gsub("SRR2407", "", sample_ids))

condition_labels <- case_when(
  srr_numbers >= 1876 & srr_numbers <= 1895 ~ "Lean",
  srr_numbers >= 1896 & srr_numbers <= 1915 ~ "Overweight",
  srr_numbers >= 1916 & srr_numbers <= 1935 ~ "GDM"
)

# Combine into colData
metadata <- data.frame(
  row.names = sample_ids,
  Sex = sex_labels,
  Condition = factor(condition_labels, levels = c("Lean", "Overweight", "GDM"))
)

# 4. Run Unified DESeq2
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = metadata,
                              design = ~ Sex + Condition)
dds <- DESeq(dds)

# 5. Generate the Validated PCA Plot (Figure 1B)
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("Sex", "Condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 6. Plotting
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = Condition, shape = Sex)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Lean" = "#44AA99", "Overweight" = "#DDCC77", "GDM" = "#CC6677")) +
  theme_bw() +
  labs(title = "Figure 1B: Validated PCA of ECFC Transcriptomes",
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance"))

# Save result
if (!dir.exists("DESeq2_All_60")) dir.create("DESeq2_All_60")
ggsave("DESeq2_All_60/Figure1B_Validated_PCA.png", plot = pca_plot, width = 8, height = 6)

message("Analysis complete. The PCA now uses the exact Sex/Sample mapping from your files.")
