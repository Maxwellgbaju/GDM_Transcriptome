# ==============================================================================
# GSVA Analysis Script (RNA-seq) - FINAL VERSION
# ==============================================================================

# 1. PATH FIXES: Force R to use the Conda library and local C++ libraries
# This avoids the "CXXABI" and "Not Writable" errors
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv("CONDA_PREFIX"), "/lib:", Sys.getenv("LD_LIBRARY_PATH"), sep=""))
.libPaths(unique(c(Sys.getenv("R_LIBS_USER"), .libPaths())))

# 2. LOAD LIBRARIES
message("Loading libraries...")
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db) 
library(GSVA)
library(GSEABase)

# 3. LOAD DATA
message("Loading merged counts and metadata...")
# Assuming your files are named exactly as previously discussed
counts <- read.csv("Merged_Counts_Master_Final.csv", row.names=1)
metadata <- read.csv("metadata.csv", row.names=1)

# Fix Alignment: Ensure metadata rows match count matrix columns exactly
metadata <- metadata[colnames(counts), , drop=FALSE]

if(!all(colnames(counts) == rownames(metadata))) {
    stop("Error: Metadata and Count Matrix sample IDs do not match!")
}

# 4. NORMALIZATION (DESeq2 VST)
message("Performing Variance Stabilizing Transformation (VST)...")
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Condition)
vsd <- vst(dds, blind=FALSE)
norm_counts <- assay(vsd)

# 5. GENE ID CONVERSION (Ensembl to Symbol)
message("Converting Ensembl IDs to Gene Symbols...")

# Remove version numbers (e.g., ENSG000001.1 -> ENSG000001)
clean_ids <- gsub("\\..*", "", rownames(norm_counts))

mapping <- bitr(clean_ids, fromType = "ENSEMBL",
                toType = "SYMBOL",
                OrgDb = org.Hs.eg.db)

# Merge symbols with the normalized counts
norm_df <- as.data.frame(norm_counts)
norm_df$ENSEMBL <- clean_ids
norm_merged <- merge(mapping, norm_df, by = "ENSEMBL")

# Handle duplicate symbols (average them)
norm_final <- aggregate(. ~ SYMBOL, data = norm_merged[,-1], FUN = mean)
rownames(norm_final) <- norm_final$SYMBOL
norm_final_matrix <- as.matrix(norm_final[,-1])

# 6. RUN GSVA
message("Running GSVA calculation...")
# UPDATED FILENAME: Matching your wget download
gene_sets <- getGmt("h.all.v2023.2.Hs.symbols.gmt")

gsva_results <- gsva(norm_final_matrix, 
                     gene_sets, 
                     method="gsva", 
                     kcdf="Gaussian", 
                     min.sz=10, 
                     max.sz=500)

# 7. SAVE RESULTS
message("Saving results to GSVA_Pathway_Scores.csv...")
write.csv(gsva_results, "GSVA_Pathway_Scores.csv")
message("Analysis Complete! You can now look at GSVA_Pathway_Scores.csv")
