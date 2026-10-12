# ============================================================
# STANDALONE SCRIPT: Slide figure — E2F Targets vs G2M Checkpoint
# Simple two-bar visual for the sex×GDM interaction test slide.
# Pulls directly from the existing GSEA interaction-ranked output.
# Run via: Rscript Slide_Interaction_Barplot.R
# ============================================================

for (pkg in c("dplyr", "readr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(dplyr); library(readr); library(ggplot2) })

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
out_dir  <- file.path(base_dir, "Interaction_Model_Output")

path_gsea_results <- file.path(out_dir, "GSEA_Interaction_Ranked_Results.csv")

# ---- Load the existing GSEA interaction results ----
gsea_res <- read_csv(path_gsea_results, show_col_types = FALSE)

# ---- Pull just the two pathways needed for the slide ----
target_pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")

slide_data <- gsea_res %>%
  filter(pathway %in% target_pathways) %>%
  mutate(
    pathway_label = gsub("HALLMARK_", "", gsub("_", " ", pathway)),
    pathway_label = tools::toTitleCase(tolower(pathway_label)),
    fdr_label = formatC(padj, format = "e", digits = 1)
  ) %>%
  arrange(desc(NES))

cat("=== Slide data: E2F Targets and G2M Checkpoint (interaction-ranked GSEA) ===\n")
print(slide_data %>% select(pathway_label, NES, pval, padj))

# ---- Simple, clean two-bar plot for a slide ----
p <- ggplot(slide_data, aes(x = reorder(pathway_label, NES), y = NES, fill = pathway_label)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.4) +
  geom_text(aes(label = paste0("FDR = ", fdr_label)),
            hjust = -0.1, size = 5, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("#C2185B", "#8E24AA")) +
  scale_y_continuous(limits = c(0, max(slide_data$NES) * 1.6)) +
  labs(
    title = "Top Sex-Differential Pathways (Interaction Test)",
    subtitle = "Ranked by sex \u00d7 GDM interaction statistic (GSEA, 50 Hallmark pathways)",
    x = NULL,
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13, color = "grey30"),
    axis.text.y = element_text(size = 15, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "Slide_Interaction_TopPathways.png"), p, width = 9, height = 4, dpi = 300)

cat("\nSlide figure saved to:\n")
cat(" -", file.path(out_dir, "Slide_Interaction_TopPathways.png"), "\n")
