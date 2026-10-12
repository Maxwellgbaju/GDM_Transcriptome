for (pkg in c("dplyr", "readr", "ggplot2", "tidyr", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}
suppressMessages({ library(dplyr); library(readr); library(ggplot2); library(tidyr); library(patchwork) })

base_dir <- "/shared/rc/malariarna/ECFC/gene_count"
path_gsva     <- file.path(base_dir, "ALL_COUNTS", "GSVA_Pathway_Scores.csv")
path_metadata <- file.path(base_dir, "METADATA_For_Sample_Sex.csv")
out_dir <- file.path(base_dir, "Supp_Figure_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_id_col_meta  <- "run"
sex_col_meta        <- "sex"
condition_col_meta  <- "group"
level_gdm    <- "gestational diabetes"
level_lean   <- "lean"
level_female <- "female"
level_male   <- "male"
target_pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")

gsva_raw <- read_csv(path_gsva, show_col_types = FALSE) %>% as.data.frame()
pathway_col1_name <- names(gsva_raw)[1]

metadata <- read_csv(path_metadata, show_col_types = FALSE) %>% as.data.frame()
metadata[[sample_id_col_meta]] <- trimws(metadata[[sample_id_col_meta]])
metadata[[sex_col_meta]]       <- tolower(trimws(metadata[[sex_col_meta]]))
metadata[[condition_col_meta]] <- tolower(trimws(metadata[[condition_col_meta]]))

all_scores <- list()
for (pw in target_pathways) {
  pathway_row <- gsva_raw[gsva_raw[[pathway_col1_name]] == pw, ]
  scores_long <- pathway_row %>%
    select(-all_of(pathway_col1_name)) %>%
    pivot_longer(everything(), names_to = "sample_id", values_to = "gsva_score") %>%
    mutate(pathway = pw)
  all_scores[[pw]] <- scores_long
}

scores_combined <- bind_rows(all_scores)
scores_combined$sample_id <- trimws(scores_combined$sample_id)

merged <- scores_combined %>%
  inner_join(metadata, by = c("sample_id" = sample_id_col_meta)) %>%
  filter(.data[[condition_col_meta]] %in% c(level_lean, level_gdm)) %>%
  mutate(
    Group = paste0(
      ifelse(.data[[sex_col_meta]] == level_female, "Female", "Male"), "-",
      ifelse(.data[[condition_col_meta]] == level_gdm, "GDM", "Lean")
    ),
    Group = factor(Group, levels = c("Female-Lean", "Female-GDM", "Male-Lean", "Male-GDM")),
    pathway_label = gsub("HALLMARK_", "", gsub("_", " ", pathway))
  )

write_csv(merged, file.path(out_dir, "MaleFemale_GSVA_TwoPathways.csv"))

# Solid, fully saturated colors -- no alpha anywhere, so nothing gets rasterized
group_colors <- c("Female-Lean" = "#F6C9DA", "Female-GDM" = "#C2185B",
                   "Male-Lean" = "#BFDDF2", "Male-GDM" = "#1565C0")

make_panel <- function(df, title_text) {
  ggplot(df, aes(x = Group, y = gsva_score, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.65, color = "black", linewidth = 0.4) +
    geom_point(position = position_jitter(width = 0.15), size = 1.6, shape = 21,
               fill = "grey30", color = "black", stroke = 0.2) +
    scale_fill_manual(values = group_colors) +
    labs(title = title_text, y = "GSVA Score", x = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 25, hjust = 1),
          plot.title = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank())
}

panel_list <- lapply(target_pathways, function(pw) {
  df_sub <- merged %>% filter(pathway == pw)
  title_text <- gsub("HALLMARK_", "", gsub("_", " ", pw))
  make_panel(df_sub, title_text)
})

combined_plot <- wrap_plots(panel_list, ncol = 2) +
  plot_annotation(title = "GSVA Pathway Activity by Fetal Sex and GDM Status",
                   theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

ggsave(file.path(out_dir, "Supp_Fig_Male_Female_GSVA.png"), combined_plot, width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "Supp_Fig_Male_Female_GSVA.pdf"), combined_plot, width = 10, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "Supp_Fig_Male_Female_GSVA.eps"), combined_plot, width = 10, height = 5, device = cairo_ps)

cat("Done. PNG, PDF, and EPS saved to:", out_dir, "\n")
