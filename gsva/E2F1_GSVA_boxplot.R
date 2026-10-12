# ============================================================
# STANDALONE FIGURE SCRIPT
# Purpose: Plot GSVA pathway activity scores (HALLMARK_E2F_TARGETS)
#          across the four sex/condition groups as a boxplot.
# This is a NEW, ADDITIONAL figure — independent of all previous scripts.
# Run via: Rscript E2F1_GSVA_boxplot.R
# ============================================================

for (pkg in c("dplyr", "readr", "ggplot2", "tidyr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr) })

# ============================================================
# ---- 1. CONFIG ----
# ============================================================

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"

path_gsva     <- file.path(base_dir, "ALL_COUNTS", "GSVA_Pathway_Scores.csv")
path_metadata <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")

out_dir <- file.path(base_dir, "E2F1_GSVA_Boxplot_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_id_col_meta  <- "run"
sex_col_meta        <- "sex"
condition_col_meta  <- "group"
level_gdm    <- "gestational diabetes"
level_lean   <- "lean"
level_female <- "female"
level_male   <- "male"

target_pathway <- "HALLMARK_E2F_TARGETS"

# ============================================================
# ---- 2. LOAD THE GSVA FILE ----
# ============================================================

gsva_raw <- read_csv(path_gsva, show_col_types = FALSE) %>% as.data.frame()

cat("=== GSVA file loaded ===\n")
cat("Dimensions:", nrow(gsva_raw), "rows x", ncol(gsva_raw), "columns\n\n")

# ============================================================
# ---- 3. EXTRACT THE TARGET PATHWAY'S SCORES PER SAMPLE ----
# (Layout confirmed: pathways are ROWS, first column holds pathway names,
#  remaining columns are samples named by SRR accession)
# ============================================================

pathway_col1_name <- names(gsva_raw)[1]  # this is "...1" in your file

pathway_row <- gsva_raw[gsva_raw[[pathway_col1_name]] == target_pathway, ]

if (nrow(pathway_row) == 0) {
  stop(paste0("target_pathway '", target_pathway, "' not found in column 1 of the GSVA file. ",
              "Check spelling exactly against the printed candidate list."))
}

scores_long <- pathway_row %>%
  select(-all_of(pathway_col1_name)) %>%
  pivot_longer(everything(), names_to = "sample_id", values_to = "gsva_score")

cat("=== GSVA scores extracted for:", target_pathway, "===\n")
print(head(scores_long, 10))
cat("n samples with scores:", nrow(scores_long), "\n\n")

# ============================================================
# ---- 4. MERGE WITH METADATA (SEX + CONDITION) ----
# ============================================================

metadata <- read_csv(path_metadata, show_col_types = FALSE) %>% as.data.frame()
metadata[[sample_id_col_meta]] <- trimws(metadata[[sample_id_col_meta]])
metadata[[sex_col_meta]]       <- tolower(trimws(metadata[[sex_col_meta]]))
metadata[[condition_col_meta]] <- tolower(trimws(metadata[[condition_col_meta]]))

scores_long$sample_id <- trimws(scores_long$sample_id)

merged <- scores_long %>%
  inner_join(metadata, by = c("sample_id" = sample_id_col_meta)) %>%
  filter(.data[[condition_col_meta]] %in% c(level_lean, level_gdm)) %>%
  mutate(
    Group = paste0(
      ifelse(.data[[sex_col_meta]] == level_female, "Female", "Male"), "-",
      ifelse(.data[[condition_col_meta]] == level_gdm, "GDM", "Lean")
    ),
    Group = factor(Group, levels = c("Female-Lean", "Female-GDM", "Male-Lean", "Male-GDM"))
  )

cat("=== Merged data — samples per group ===\n")
print(table(merged$Group))
cat("\n")

if (nrow(merged) == 0) {
  stop("No samples matched after merging GSVA scores with metadata — check that sample IDs align (e.g. SRR format) between the two files.")
}

# ============================================================
# ---- 5. STATISTICS: GROUP MEANS + WILCOX TEST ----
# ============================================================

cat("=== Group summary statistics ===\n")
merged %>% group_by(Group) %>%
  summarise(mean_score = mean(gsva_score), sd_score = sd(gsva_score), n = n()) %>%
  print()
cat("\n")

fl <- merged %>% filter(Group == "Female-Lean") %>% pull(gsva_score)
fg <- merged %>% filter(Group == "Female-GDM") %>% pull(gsva_score)
ml <- merged %>% filter(Group == "Male-Lean")  %>% pull(gsva_score)
mg <- merged %>% filter(Group == "Male-GDM")   %>% pull(gsva_score)

cat("Wilcoxon test, Female-GDM vs Female-Lean: p =", wilcox.test(fg, fl)$p.value, "\n")
cat("Wilcoxon test, Male-GDM vs Male-Lean:     p =", wilcox.test(mg, ml)$p.value, "\n\n")

# ============================================================
# ---- 6. PLOT — BOXPLOT WITH INDIVIDUAL POINTS ----
# ============================================================

group_colors <- c("Female-Lean" = "#F2C2D6", "Female-GDM" = "#C2185B",
                   "Male-Lean" = "#B3D9F2", "Male-GDM" = "#1565C0")

p <- ggplot(merged, aes(x = Group, y = gsva_score, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = group_colors) +
  labs(title = "GSVA Pathway Activity Score", subtitle = target_pathway,
       y = "GSVA Score", x = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(file.path(out_dir, "E2F1_GSVA_Boxplot.tiff"), p, width = 7, height = 6, dpi = 600, compression = "lzw")
ggsave(file.path(out_dir, "E2F1_GSVA_Boxplot.png"),  p, width = 7, height = 6, dpi = 300)

write_csv(merged, file.path(out_dir, "GSVA_scores_by_group.csv"))

cat("Figure and data saved to:", out_dir, "\n")
