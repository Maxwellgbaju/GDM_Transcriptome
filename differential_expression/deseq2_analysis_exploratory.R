# --- DESeq2 Differential Gene Expression Analysis (Exploratory) ---
#
# This script performs the DGE analysis using an exploratory cutoff:
# Nominal P-value (pvalue) < 0.05 AND |log2FC| > 0.58.
#
# This is designed to find subtle gene expression changes, particularly for
# the GDM groups, that were masked by the strict multiple testing correction.
#
# NOTE: Ensure the required packages (DESeq2, tidyverse, pheatmap, clusterProfiler, pathview) are installed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# *** REVISED INSTALLATION COMMAND (to resolve ggtree/clusterProfiler dependency issues): ***
# BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "pathview", "ggtree", "enrichplot"), force=TRUE)
# install.packages(c("tidyverse", "pheatmap"), repos = "http://cran.us.r-project.org")

# --- AGGRESSIVE FIX: Overriding GLIBCXX Error on Cluster ---
# The GLIBCXX_3.4.30 error means R is trying to use an old system C++ library.
# This block explicitly sets the LD_LIBRARY_PATH environment variable *inside* R,
# pointing it to the correct, newer Conda library path *before* DESeq2 loads.
conda_lib_path <- file.path(Sys.getenv("CONDA_PREFIX"), "lib")
Sys.setenv(LD_LIBRARY_PATH = paste(conda_lib_path, Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
cat(paste0("LD_LIBRARY_PATH set to prioritize Conda libraries: ", Sys.getenv("LD_LIBRARY_PATH"), "\n"))
# --- END AGGRESSIVE FIX ---

# --- FIX for 'no package called’ error when using Rscript ---
# Rscript often fails to find packages installed by Conda/BiocManager.
# This block explicitly adds the Conda R library path to R's search path *before*
# any package loading attempts, which is critical for Rscript execution.
conda_r_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib", "R", "library")
if (dir.exists(conda_r_lib) && !conda_r_lib %in% .libPaths()) {
    .libPaths(c(conda_r_lib, .libPaths()))
    cat("Conda R library path added to search path.\n")
}
# --- END FIX ---

# Load Libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler) 
library(org.Hs.eg.db)    
library(pathview)        # NEW: For detailed KEGG pathway visualization

# --- 1. CONFIGURATION AND FILE LOADING ---

# Define the base directory where the input files are located
# UPDATED PATH: Changed from /shared/rc/malariarna/ECFC/gene_count to /home/msabms/Pregnancy_Project/gene_count
base_dir <- "/home/msabms/Pregnancy_Project/gene_count"

# Define file paths using the full path to load the data
count_file <- file.path(base_dir, "all_groups_unified_counts.tsv")
meta_file <- file.path(base_dir, "all_groups_metadata.tsv")

# The output directory has been changed to a separate folder for exploratory results.
output_dir <- "./exploratory_result"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
    dir.create(output_dir)
}

# Load Count Data
# Note: Skip the first row if it's a general header, but assuming the first line
# is Geneid\tSRR... and it's tab-separated.
count_data <- read.delim(count_file, header = TRUE, sep = "\t", row.names = 1)

# Ensure all counts are integers (DESeq2 requirement)
count_data <- round(count_data)
count_data[is.na(count_data)] <- 0 # Replace any NA introduced by rounding/loading with 0

# Load Metadata
metadata <- read.delim(meta_file, header = TRUE, sep = "\t")

# Check: Ensure the metadata SampleID matches and is in the same order as the
# column names in the count data.
if (!all(metadata$SampleID == colnames(count_data))) {
    # If they don't match, reorder the count data columns based on metadata
    print("Reordering count matrix columns to match metadata...")
    count_data <- count_data[, metadata$SampleID]
}

# --- 2. PRE-FILTERING (Addressing the Galaxy Issue) ---

# Low-count filtering is crucial for statistical power. We filter for genes
# that have at least 10 reads in at least 3 samples (a conservative threshold).
keep <- rowSums(count_data >= 10) >= 3
count_data_filtered <- count_data[keep, ]

cat(paste0("Original number of genes: ", nrow(count_data), "\n"))
cat(paste0("Number of genes after filtering: ", nrow(count_data_filtered), "\n"))

# --- 3. DESeq2 SETUP AND RUN ---

# Set the factor (Condition) as a factor variable
metadata$Condition <- factor(metadata$Condition)

# IMPORTANT: Set 'Lean' as the reference level (control group)
metadata$Condition <- relevel(metadata$Condition, ref = "Lean")

# Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
    countData = count_data_filtered,
    colData = metadata,
    design = ~Condition
)

# Run the DESeq2 analysis pipeline
cat("\nRunning DESeq2 analysis...\n")
dds <- DESeq(dds)
cat("DESeq2 analysis complete.\n")

# Initialize list to collect all significant GeneIDs for the heatmap
all_significant_genes <- c()

# Global variables for Pathview analysis (will be populated in the loop)
gdm_lean_results_df <- NULL
kegg_results <- NULL

# --- 4. EXTRACTING AND SAVING RESULTS AND PLOTTING ---

# Define the three contrasts to run:
contrasts_to_run <- list(
    c("Condition", "GDM", "Lean"),
    c("Condition", "Overweight", "Lean"),
    c("Condition", "GDM", "Overweight")
)

# Set the p-value column to filter on: 'padj' (strict) or 'pvalue' (exploratory)
P_VALUE_COLUMN <- "pvalue" # SET TO NOMINAL P-VALUE FOR EXPLORATORY ANALYSIS

for (contrast in contrasts_to_run) {
    comparison_name <- paste(contrast[2], "_vs_", contrast[3], sep = "")
    
    # Extract results
    res <- results(dds, contrast = contrast, alpha = 0.05)
    
    # Order results by adjusted p-value
    res_ordered <- res[order(res$padj), ]
    
    # Convert to a tidy data frame
    res_df <- as.data.frame(res_ordered) %>%
        rownames_to_column(var = "GeneID")
        
    # --- Store GDM vs Lean results for Pathview analysis later (Section 8) ---
    if (comparison_name == "GDM_vs_Lean") {
        gdm_lean_results_df <- res_df
    }
    
    # --- Save the complete results table ---
    output_path_full <- file.path(output_dir, paste0("DGE_results_", comparison_name, "_FULL.tsv"))
    # Note: We skip saving the full results here to avoid duplication, as they are saved in the primary script.
    
    # --- Filter for Significant Genes ---
    # Define significance as log2FC > |0.58| (1.5-fold) AND pvalue < 0.05
    p_cutoff <- res_df[[P_VALUE_COLUMN]]
    p_label <- paste0("|log2FC| > 0.58 & ", P_VALUE_COLUMN, " < 0.05")
    
    res_sig <- subset(res_df, p_cutoff < 0.05 & abs(log2FoldChange) > 0.58)
    
    # Collect significant gene IDs for later heatmap generation
    all_significant_genes <- c(all_significant_genes, res_sig$GeneID)
    
    output_path_sig <- file.path(output_dir, paste0("DGE_results_", comparison_name, "_SIGNIFICANT_FILTERED_BY_", P_VALUE_COLUMN, ".tsv"))
    write.table(res_sig, output_path_sig, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(paste0("Found ", nrow(res_sig), " significant genes for ", comparison_name, " (using ", P_VALUE_COLUMN, ") and saved to: ", output_path_sig, "\n\n"))
    
    # --- Generate Volcano Plot ---
    # Prepare data for plotting
    res_plot <- res_df %>%
        mutate(
            Significant = case_when(
                p_cutoff < 0.05 & log2FoldChange > 0.58 ~ "Upregulated",
                p_cutoff < 0.05 & log2FoldChange < -0.58 ~ "Downregulated",
                TRUE ~ "Not Significant"
            ),
            # Use the calculated p-value column for the y-axis
            log10_p = -log10(p_cutoff)
        )
    
    # Create the plot
    volcano_plot <- ggplot(res_plot, aes(x = log2FoldChange, y = log10_p)) +
        geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # p-value cutoff line
        geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black") + # log2FC cutoff lines
        labs(
            title = paste("Volcano Plot:", comparison_name, " (Exploratory)"),
            x = expression(log[2]~"Fold Change"),
            y = expression(-log[10]~P-Value~" (", !!P_VALUE_COLUMN, ")"), # Dynamic y-axis label
            caption = paste("Significance:", p_label)
        ) +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    # Save the plot
    plot_path <- file.path(output_dir, paste0("Volcano_Plot_", comparison_name, "_FILTERED_BY_", P_VALUE_COLUMN, ".png"))
    ggsave(plot_path, plot = volcano_plot, width = 8, height = 6, units = "in", dpi = 300)
    cat(paste0("Generated Volcano Plot for ", comparison_name, " and saved to: ", plot_path, "\n\n"))
}

# --- 5. HEATMAP GENERATION (Clustering of Significant Genes) ---

# Get unique list of all significant genes found across all comparisons
unique_sig_genes <- unique(all_significant_genes)

if (length(unique_sig_genes) > 1) {
    cat(paste0("Generating Heatmap for ", length(unique_sig_genes), " unique significant genes...\n"))
    
    # Calculate VST normalized counts (needed for Heatmap)
    vsd <- vst(dds, blind = FALSE)
    vst_counts <- assay(vsd)
    
    # Filter the VST counts to include only the significant genes
    heatmap_matrix <- vst_counts[unique_sig_genes, ]
    
    # Scale the expression data for better visualization (Z-score scaling per gene)
    
    # Prepare the annotation data frame for the columns (samples)
    annotation_col <- metadata %>%
        select(Condition) %>%
        as.data.frame()
    rownames(annotation_col) <- metadata$SampleID
    
    # Define colors for the Condition annotation bar
    ann_colors = list(
        Condition = c(Lean = "darkgreen", Overweight = "orange", GDM = "red")
    )
    
    # Save the Heatmap
    heatmap_path <- file.path(output_dir, paste0("Heatmap_Top_DGE_Genes_FILTERED_BY_", P_VALUE_COLUMN, ".png"))
    
    # Open PNG device
    png(heatmap_path, width = 10, height = 12, units = "in", res = 300)
    
    pheatmap(
        heatmap_matrix,
        cluster_rows = TRUE,
        cluster_cols = TRUE, # Cluster samples
        show_colnames = FALSE, # Hide sample names for clarity
        show_rownames = (length(unique_sig_genes) < 50), # Show gene names only if list is small
        scale = "row", # Scale expression values by row (Z-score for each gene)
        border_color = NA,
        annotation_col = annotation_col,
        annotation_colors = ann_colors,
        main = paste("Hierarchical Clustering of Differentially Expressed Genes (Filtered by", P_VALUE_COLUMN, ")")
    )
    
    # Close PNG device
    dev.off()
    
    cat(paste0("Generated Heatmap and saved to: ", heatmap_path, "\n"))
    
    # --- 7. FUNCTIONAL ENRICHMENT ANALYSIS (GO and KEGG) ---
    cat("\n--- 7. FUNCTIONAL ENRICHMENT ANALYSIS (GO and KEGG) ---\n")

    if (length(unique_sig_genes) > 0) {
        
        # 7a. Convert Gene Symbols to Entrez IDs (required by enrichGO and enrichKEGG)
        # Note: 'drop = TRUE' removes genes that couldn't be mapped
        entrez_ids <- bitr(
            unique_sig_genes, 
            fromType = "SYMBOL", 
            toType = "ENTREZID", 
            OrgDb = org.Hs.eg.db,
            drop = TRUE
        )
        
        if (nrow(entrez_ids) > 0) {
            
            # --- GO Enrichment (Biological Process) ---
            cat("\n--- Running GO Enrichment (Biological Process) ---\n")
            go_results <- enrichGO(
                gene = entrez_ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = 'ENTREZID',
                ont = "BP", # Biological Process
                pAdjustMethod = "BH", # Benjamini-Hochberg correction
                qvalueCutoff = 0.05
            )
            
            # Save the GO results
            go_df <- as.data.frame(go_results)
            output_path_go <- file.path(output_dir, paste0("GO_Enrichment_Results_", P_VALUE_COLUMN, ".tsv"))
            write.table(go_df, output_path_go, sep = "\t", quote = FALSE, row.names = FALSE)
            cat(paste0("Saved GO Enrichment table to: ", output_path_go, "\n"))
            
            # Generate GO Dot Plot for visualization (Top 10 enriched terms)
            if (nrow(go_df) > 0) {
                dotplot_go <- dotplot(go_results, showCategory = 10) + 
                    ggtitle(paste0("GO Biological Process Enrichment (Filtered by ", P_VALUE_COLUMN, ")"))
                
                plot_path_dot <- file.path(output_dir, paste0("GO_DotPlot_", P_VALUE_COLUMN, ".png"))
                ggsave(plot_path_dot, plot = dotplot_go, width = 10, height = 7, units = "in", dpi = 300)
                cat(paste0("Generated GO Dot Plot and saved to: ", plot_path_dot, "\n"))
            } else {
                 cat("No GO terms enriched at q-value < 0.05.\n")
            }
            
            # --- KEGG Enrichment ---
            cat("\n--- Running KEGG Enrichment ---\n")
            kegg_results <- enrichKEGG(
                gene = entrez_ids$ENTREZID,
                organism = 'hsa', # Homo sapiens
                pvalueCutoff = 0.05
            )

            # Save the KEGG results
            kegg_df <- as.data.frame(kegg_results)
            output_path_kegg <- file.path(output_dir, paste0("KEGG_Enrichment_Results_", P_VALUE_COLUMN, ".tsv"))
            write.table(kegg_df, output_path_kegg, sep = "\t", quote = FALSE, row.names = FALSE)
            cat(paste0("Saved KEGG Enrichment table to: ", output_path_kegg, "\n"))

            # Generate KEGG Dot Plot for visualization (Top 10 enriched terms)
            if (nrow(kegg_df) > 0) {
                 dotplot_kegg <- dotplot(kegg_results, showCategory = 10) + 
                    ggtitle(paste0("KEGG Pathway Enrichment (Filtered by ", P_VALUE_COLUMN, ")"))
                
                plot_path_kegg_dot <- file.path(output_dir, paste0("KEGG_DotPlot_", P_VALUE_COLUMN, ".png"))
                ggsave(plot_path_kegg_dot, plot = dotplot_kegg, width = 10, height = 7, units = "in", dpi = 300)
                cat(paste0("Generated KEGG Dot Plot and saved to: ", plot_path_kegg_dot, "\n"))
            } else {
                 cat("No KEGG terms enriched at p-value < 0.05.\n")
            }
            
            # --- 8. KEGG PATHWAY MAPPING WITH PATHVIEW (GDM_vs_Lean) ---
            if (nrow(kegg_df) > 0 && !is.null(gdm_lean_results_df)) {
                cat("\n--- 8. KEGG PATHWAY MAPPING (PATHVIEW) ---\n")
                
                # 8a. Prepare the log2FC vector for pathview
                # Create a named vector: names are Entrez IDs, values are log2FC
                lfc_vector <- gdm_lean_results_df %>%
                    select(GeneID, log2FoldChange) %>%
                    # Map Gene Symbols to Entrez IDs using bitr again for all DGE genes
                    # Note: Using the whole gene list is better for pathview coloring
                    right_join(entrez_ids, by = c("GeneID" = "SYMBOL")) %>% 
                    drop_na(ENTREZID, log2FoldChange) %>%
                    # Ensure LFC values are numeric and use EntrezID as the name
                    {
                        lfc <- .$log2FoldChange
                        names(lfc) <- .$ENTREZID
                        lfc
                    }
                
                # 8b. Select the top N significant KEGG pathway IDs
                top_kegg_ids <- head(kegg_df$ID, 3) # Visualize the top 3 pathways
                
                cat(paste0("Mapping log2FC (GDM vs Lean) onto the top ", length(top_kegg_ids), " KEGG pathways...\n"))
                
                # Set output directory for pathview files
                # Pathview generates files named like 'hsa05200.pathview.png' in the current working directory,
                # so we must temporarily change the working directory to save them in our desired output folder.
                original_wd <- getwd()
                setwd(output_dir)
                
                tryCatch({
                    pathview(
                        gene.data = lfc_vector, 
                        pathway.id = top_kegg_ids, 
                        species = "hsa", 
                        out.suffix = paste0("GDM_vs_Lean_", P_VALUE_COLUMN),
                        keys.align = "y", # Aligns keys/labels nicely
                        kegg.native = TRUE # Uses KEGG's native format
                    )
                    cat(paste0("Generated detailed KEGG map images (PNG/PDF) for pathway IDs: ", paste(top_kegg_ids, collapse = ", "), "\n"))
                }, error = function(e) {
                    cat(paste("An error occurred during pathview generation:", conditionMessage(e), "\n"))
                }, finally = {
                    setwd(original_wd) # Always change back to the original working directory
                })

            } else {
                 cat("Skipping Pathview: No KEGG pathways enriched or GDM_vs_Lean results missing.\n")
            }
            # --- END OF PATHVIEW SECTION ---

        } else {
            cat("Skipping GO/KEGG Enrichment: Failed to map significant Gene Symbols to Entrez IDs. Check if 'SYMBOL' is the correct ID type.\n")
        }
    } else {
        cat("Skipping GO/KEGG Enrichment: Fewer than 1 significant gene found.\n")
    }
}

# --- 6. OPTIONAL: Save VST Normalized Counts for QC/Plotting ---

# Note: We skip saving VST counts again as they are identical to the previous run.

# Final summary of output files
cat("\n--- EXPLORATORY ANALYSIS COMPLETE ---\n")
cat(paste0("Output files saved in the '", output_dir, "' directory.\n"))

