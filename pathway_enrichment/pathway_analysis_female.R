# ==============================================================================
# Pathway Analysis: GDM vs. Lean (FEMALES ONLY) - Publication Quality
# ==============================================================================
.libPaths(unique(c(Sys.getenv("R_LIBS_USER"), .libPaths())))
library(limma)
library(ggplot2)

# 1. LOAD DATA
metadata <- read.csv("metadata.csv") 
gsva_matrix <- as.matrix(read.csv("GSVA_Pathway_Scores.csv", row.names=1))

# 2. MATCH IDs AND SUBSET FEMALES
metadata[,1] <- make.names(metadata[,1])
rownames(metadata) <- metadata[,1]
female_samples <- rownames(metadata[trimws(tolower(metadata$Sex)) == "female", ])
final_female_ids <- intersect(female_samples, colnames(gsva_matrix))

# Create the subset matrix and metadata
gsva_female <- gsva_matrix[, final_female_ids]
meta_female <- metadata[final_female_ids, ]

# 3. SAVE THE FEMALE-ONLY INPUT DATA
write.csv(gsva_female, "Female_Only_GSVA_Scores_Input.csv")

# 4. STATISTICAL ANALYSIS
meta_female$Condition <- gsub(" ", "_", trimws(tolower(meta_female$Condition)))
group <- as.factor(meta_female$Condition)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(gsva_female, design)
contrast_matrix <- makeContrasts(gestational_diabetes - lean, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 5. RESULTS TABLE GENERATION
stats_table <- topTable(fit2, number=Inf, adjust.method="BH")
stats_table$Pathway <- gsub("HALLMARK_", "", rownames(stats_table))
stats_table$Pathway <- gsub("_", " ", stats_table$Pathway)
stats_table <- stats_table[stats_table$Pathway != "SPERMATOGENESIS", ]

# Save the full statistics for females
write.csv(stats_table, "Female_Only_Stats_Full.csv")

# 6. PLOT GENERATION
top_10 <- head(stats_table[order(stats_table$P.Value), ], 10)

plot <- ggplot(top_10, aes(x=reorder(Pathway, logFC), y=logFC, fill=logFC > 0)) +
  geom_bar(stat="identity", color="black", width=0.7) +
  coord_flip() +
  scale_fill_manual(values = c("skyblue3", "indianred2"), 
                    labels = c("Downregulated", "Upregulated")) +
  labs(title="Pathway Activity in Female GDM vs. Lean",
       subtitle="Positive logFC indicates higher activity in Gestational Diabetes",
       x="Biological Pathway", y="Activity Score (logFC)", fill="Direction") +
  theme_minimal() +
  theme(axis.text.y = element_text(face="bold"))

# 7. EXPORT FORMATS
ggsave("Female_GDM_Pathways_Final.png", plot, width=10, height=7, dpi=300)
ggsave("Female_GDM_Pathways_Final.eps", plot, width=10, height=7, device="eps")

message(paste("--- SUCCESS: Found", length(final_female_ids), "female samples ---"))
message("Files generated: .png, .eps, and two .csv files.")
