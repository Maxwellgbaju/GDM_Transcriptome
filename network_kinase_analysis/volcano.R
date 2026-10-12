suppressPackageStartupMessages({
  library(biomaRt)
  library(decoupleR)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
})

# Silence the specific OmnipathR deprecation warning to clean up terminal
options(warn = -1)

cat("Starting Robust Pipeline v5.3...\n")

# 1. Load Raw Data
if(!file.exists("DGE_Results.csv")) stop("DGE_Results.csv not found!")
dge_res <- read.csv("DGE_Results.csv")
colnames(dge_res)[1] <- "ensembl_id"
dge_res$ensembl_clean <- gsub("\\..*", "", dge_res$ensembl_id)

# Auto-detect columns
p_col <- intersect(colnames(dge_res), c("pvalue", "p_value", "padj", "PValue"))[1]
fc_col <- intersect(colnames(dge_res), c("log2FoldChange", "logFC", "log2FC"))[1]

# 2. ID Mapping
cat("Step 2: Mapping IDs via biomaRt...\n")
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://www.ensembl.org")
mapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'ensembl_gene_id',
                 values = dge_res$ensembl_clean,
                 mart = ensembl)
final_data <- merge(dge_res, mapping, by.x = "ensembl_clean", by.y = "ensembl_gene_id")

# 3. Regulatory Inference
cat("Step 3: Running ULM (CollecTRI)...\n")
net <- get_collectri(organism='human', split_complexes=FALSE)
mat <- final_data %>%
  filter(hgnc_symbol != "" & !is.na(hgnc_symbol)) %>%
  group_by(hgnc_symbol) %>%
  slice_max(order_by = abs(stat), n = 1) %>%
  ungroup() %>%
  column_to_rownames("hgnc_symbol") %>%
  select(stat) %>%
  as.matrix()

tf_acts <- run_ulm(mat = mat, net = net, .source='source', .target='target', .mor='mor')

# 4. Identify E2F1 Targets
e2f1_targets <- net %>% filter(source == "E2F1") %>% pull(target)

# 5. Generate Informative Volcano Plot
cat("Step 5: Generating Informative Volcano Plot...\n")
plot_data <- final_data %>%
  mutate(
    P_VAL = !!sym(p_col),
    LOG_FC = !!sym(fc_col),
    type = case_when(
      hgnc_symbol %in% e2f1_targets & P_VAL < 0.05 ~ "E2F1 Target",
      P_VAL < 0.05 ~ "Significant",
      TRUE ~ "Non-Significant"
    )
  )

volcano <- ggplot(plot_data, aes(x = LOG_FC, y = -log10(P_VAL))) +
  geom_point(aes(color = type), alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("E2F1 Target" = "#E41A1C", "Significant" = "#377EB8", "Non-Significant" = "grey80")) +
  geom_text_repel(data = subset(plot_data, type == "E2F1 Target" & abs(LOG_FC) > 0.4),
                  aes(label = hgnc_symbol), size = 3.5, fontface = "bold", 
                  max.overlaps = 25, box.padding = 0.6) +
  theme_minimal() +
  labs(title = "Placental Transcriptomic Landscape (Female GDM)",
       subtitle = "Red highlights: Validated targets of master regulator E2F1",
       x = "Log2 Fold Change", y = "-Log10 P-value")

ggsave("Final_Informative_Volcano.png", volcano, width = 8, height = 7, dpi = 300)

cat("Success! 'Final_Informative_Volcano.png' saved.\n")
options(warn = 0)
