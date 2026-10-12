suppressPackageStartupMessages({
  library(biomaRt)
  library(decoupleR)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

cat("Starting Validation Pipeline (Version 4.2)...\n")

# 1. Load Data
dge_res <- read.csv("DGE_Results.csv")
colnames(dge_res)[1] <- "ensembl_id"
dge_res$ensembl_clean <- gsub("\\..*", "", dge_res$ensembl_id)

# 2. ID Mapping
cat("Step 2: Mapping IDs via biomaRt...\n")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'ensembl_gene_id',
                 values = dge_res$ensembl_clean,
                 mart = ensembl)
final_data <- merge(dge_res, mapping, by.x = "ensembl_clean", by.y = "ensembl_gene_id")

# 3. Prepare Matrix
cat("Step 3: Preparing matrix...\n")
analysis_ready <- final_data %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) %>%
  group_by(hgnc_symbol) %>%
  slice_max(order_by = abs(stat), n = 1) %>%
  ungroup() %>%
  column_to_rownames("hgnc_symbol") %>%
  select(stat) %>%
  as.matrix()

# 4. Get Network
cat("Step 4: Loading CollecTRI Network...\n")
net <- get_collectri(organism='human', split_complexes=FALSE)

# 5. Run ULM Inference
cat("Step 5: Running ULM inference...\n")
tf_acts <- run_ulm(mat = analysis_ready, net = net, .source='source', .target='target', .mor='mor')

# --- CRITICAL: FIXING THE OUTPUTS ---

# Display Top 10 in Terminal
cat("\n--- TOP 10 REGULATORY DRIVERS FOR MANUSCRIPT ---\n")
print(tf_acts %>% arrange(desc(abs(score))) %>% head(10))

# 6. Generating Visualizations
cat("\nStep 6: Generating Plots...\n")

# A. Volcano Plot
# We use p-values from the tf_acts object
volcano_plot <- ggplot(tf_acts, aes(x = score, y = -log10(p_value))) +
  geom_point(aes(color = abs(score) > 4), alpha = 0.5) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text(data = subset(tf_acts, abs(score) > 5), aes(label = source), vjust = 1.5, size = 3) +
  theme_minimal() +
  labs(title = "TF Activity Volcano Plot", x = "ULM Score", y = "-log10(p-value)")

ggsave("TF_Volcano.png", volcano_plot, width = 7, height = 6)

# B. Top 10 Bar Plot
top10_bar <- tf_acts %>%
  arrange(desc(abs(score))) %>%
  head(10) %>%
  ggplot(aes(x = reorder(source, score), y = score, fill = score)) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_minimal() +
  labs(title = "Top 10 Master Regulators", x = "Transcription Factor", y = "Activity Score")

ggsave("Top10_TFs.png", top10_bar, width = 7, height = 5)

cat("Success! Files saved: 'TF_Volcano.png' and 'Top10_TFs.png'\n")
