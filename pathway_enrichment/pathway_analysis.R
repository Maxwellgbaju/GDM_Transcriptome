# ==============================================================================
# Pathway Analysis: Gestational Diabetes vs. Lean (FEMALES ONLY) - FINAL FIX
# ==============================================================================

.libPaths(unique(c(Sys.getenv("R_LIBS_USER"), .libPaths())))
library(limma)
library(ggplot2)

# 1. LOAD DATA
metadata <- read.csv("metadata.csv") # Load without row.names first to handle ID cleaning
gsva_matrix <- as.matrix(read.csv("GSVA_Pathway_Scores.csv", row.names=1))

# 2. MATCH SAMPLE IDs (Handle the "R-naming" quirk)
# If GSVA column names start with 'X', but metadata doesn't, we fix it
meta_ids <- metadata[,1]
gsva_ids <- colnames(gsva_matrix)

if(!all(meta_ids %in% gsva_ids)) {
  message("Attempting to fix Sample ID naming mismatch...")
  metadata[,1] <- make.names(metadata[,1])
}
rownames(metadata) <- metadata[,1]

# 3. SUBSET FOR FEMALES
# trimws() handles any invisible whitespace; tolower() handles 'Female' vs 'female'
female_samples <- rownames(metadata[trimws(tolower(metadata$Sex)) == "female", ])

# Filter the matrix to keep ONLY female columns that actually exist in the data
final_female_ids <- intersect(female_samples, colnames(gsva_matrix))
message(paste("Successfully identified", length(final_female_ids), "female samples for analysis."))

if(length(final_female_ids) == 0) {
    stop("ERROR: 0 female samples found. Please check if Sample IDs match between files!")
}

gsva_female <- gsva_matrix[, final_female_ids]
meta_female <- metadata[final_female_ids, ]

# 4. CLEAN CONDITION NAMES & STATS
meta_female$Condition <- gsub(" ", "_", trimws(tolower(meta_female$Condition)))
group <- as.factor(meta_female$Condition)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(gsva_female, design)
contrast_matrix <- makeContrasts(gestational_diabetes - lean, levels=design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 5. RESULTS & EXCLUDE SPERMATOGENESIS
stats_table <- topTable(fit2, number=Inf, adjust.method="BH")
stats_table$Pathway <- gsub("HALLMARK_", "", rownames(stats_table))
stats_table$Pathway <- gsub("_", " ", stats_table$Pathway)

# Exclude the non-specific Spermatogenesis signal for a cleaner story
stats_table <- stats_table[stats_table$Pathway != "SPERMATOGENESIS", ]

# 6. PLOT TOP 10
top_10 <- head(stats_table[order(stats_table$P.Value), ], 10)

plot <- ggplot(top_10, aes(x=reorder(Pathway, logFC), y=logFC, fill=logFC > 0)) +
  geom_bar(stat="identity", color="black", width=0.7) +
  coord_flip() +
  scale_fill_manual(values = c("skyblue3", "indianred2"), labels = c("Lower in GDM", "Higher in GDM")) +
  labs(title="Top 10 Pathway Changes (FEMALES ONLY)",
       subtitle="Gestational Diabetes vs. Lean | Spermatogenesis excluded",
       x="Pathway", y="Activity Difference (logFC)", fill="Direction") +
  theme_minimal() +
  theme(axis.text.y = element_text(face="bold"))

ggsave("Female_Only_Pathways_REAL.png", plot, width=10, height=7)
write.csv(stats_table, "Female_Only_Stats_Full.csv")
message("DONE! Results saved to Female_Only_Pathways_REAL.png")
