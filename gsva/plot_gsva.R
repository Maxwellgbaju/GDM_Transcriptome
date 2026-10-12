.libPaths(unique(c(Sys.getenv("R_LIBS_USER"), .libPaths())))
library(pheatmap)

# 1. Load the scores and metadata
gsva_matrix <- read.csv("GSVA_Pathway_Scores.csv", row.names=1)
metadata <- read.csv("metadata.csv", row.names=1)

# 2. Match metadata to the matrix
metadata <- metadata[colnames(gsva_matrix), , drop=FALSE]

# 3. Create the Heatmap
# We clean up the names (remove "HALLMARK_") to make it prettier
rownames(gsva_matrix) <- gsub("HALLMARK_", "", rownames(gsva_matrix))

pheatmap(gsva_matrix, 
         annotation_col = metadata, 
         show_colnames = FALSE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "GSVA Hallmark Pathway Scores",
         filename = "GSVA_Heatmap.png",
         width = 10, height = 12)

message("Heatmap saved as GSVA_Heatmap.png")
