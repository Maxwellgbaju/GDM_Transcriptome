# 1. Load the libraries
library(DESeq2)

# 2. Load your files
counts <- read.csv("Merged_Counts_Master_Final.csv", row.names=1)
metadata <- read.csv("metadata.csv", row.names=1)

# --- THE FIX STEP (Run this now to avoid the error) ---
# This forces the metadata rows to be in the exact same order as the count columns
metadata <- metadata[colnames(counts), , drop=FALSE]

# Double check: This should now say TRUE
all(colnames(counts) == rownames(metadata))
# -----------------------------------------------------

# 3. Create the DESeq2 object (This will now work!)
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = metadata,
                              design = ~ Condition)

# 4. Apply Normalization (VST)
# This removes technical noise so the GSVA works better
vsd <- vst(dds, blind=FALSE)
norm_counts <- assay(vsd)

# 5. Save the result
write.csv(norm_counts, "Normalized_Counts_for_GSVA.csv")
