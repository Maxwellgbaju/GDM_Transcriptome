# ============================================================
# ---- 4. PART B — Top-12 heatmap, forced group clustering (REFINED) ----
# ============================================================

# ... (keep everything above unchanged through building heat_mat_ordered) ...

# Compute gap positions between the four blocks (for visual separation)
group_sizes <- sapply(group_order, function(g) sum(meta_sub$group_label[colnames(heat_mat_ordered) %in% ordered_samples] == g))
gap_positions <- cumsum(group_sizes)[1:3]  # gaps after blocks 1, 2, 3 (not after the last one)

annotation_colors <- list(
  Group = c("Female-Lean" = "#F2C2D6", "Female-GDM" = "#C2185B",
            "Male-Lean"   = "#B3D9F2", "Male-GDM"   = "#1565C0")
)

# High-resolution TIFF, journal-standard 300+ dpi
tiff(file.path(out_dir, "Top12_Heatmap_ForcedGrouping.tiff"),
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
pheatmap(
  heat_mat_ordered,
  cluster_cols     = FALSE,
  cluster_rows     = TRUE,
  gaps_col         = gap_positions,
  annotation_col   = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames    = FALSE,
  color            = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color     = NA,
  main             = "Top 12 Female GDM-vs-Lean DEGs Across Sex/Condition Groups",
  fontsize_row     = 11,
  fontsize         = 10
)
dev.off()

# Keep the PNG too, for quick viewing
png(file.path(out_dir, "Top12_Heatmap_ForcedGrouping.png"), width = 2400, height = 1800, res = 300)
pheatmap(
  heat_mat_ordered, cluster_cols = FALSE, cluster_rows = TRUE,
  gaps_col = gap_positions, annotation_col = annotation_col,
  annotation_colors = annotation_colors, show_colnames = FALSE,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  border_color = NA,
  main = "Top 12 Female GDM-vs-Lean DEGs Across Sex/Condition Groups",
  fontsize_row = 11, fontsize = 10
)
dev.off()
