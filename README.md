# GDM-ECFC Downstream Analysis Scripts

Scripts supporting the manuscript: "Sex-stratified transcriptomic analysis identifies E2F1-regulon activation and an unscheduled proliferative signature in fetal endothelial progenitors exposed to gestational diabetes."

Raw sequencing data: NCBI SRA, BioProject PRJNA952510 (60 samples, SRR24071876-SRR24071935).
Raw-data processing scripts (SRA to FASTQ, FastQC, Trimmomatic, HISAT2, SAM to BAM, Picard) are deposited separately, see the companion GDM_Pipeline_Scripts repository (steps 01-06).

This repository contains the downstream statistical and regulatory-network analysis performed after read quantification.

## Gene Quantification (documentation, not a saved script)

The featureCounts command was executed directly on the cluster and was not retained as a saved script file. Based on the Methods description and the GENCODE v44 GTF file present in the working environment, the command used was equivalent to:

featureCounts -a gencode.v44.annotation.gtf -o gene_counts.txt -T <threads> <sample>.dedup.bam

(uniquely mapped reads only, as stated in Methods; exact flag values may vary slightly from this reconstruction.)

## Folder Structure

| Folder | Contents | Manuscript Section |
|---|---|---|
| differential_expression/ | Pooled and sex-stratified DESeq2 analyses across all three metabolic groups; metadata correction script | Sections 3.2-3.4, Fig. 2, Fig. 3 |
| pathway_enrichment/ | GSEA (WebGestalt/Reactome) and pathway-analysis scripts; FDR summary table generator | Section 3.5, Fig. 4, Supplementary Table S2 |
| gsva/ | Gene Set Variation Analysis scripts (female-specific and male comparator) | Section 3.7.1, Fig. 6, Supplementary Fig. S1 |
| interaction_model/ | Formal sex x GDM interaction testing (DESeq2 factorial model, GSEA on interaction statistic, GSVA factorial model), PCA/PERMANOVA, male-female correlation analysis | Section 3.2 (PCA), Section 3.10, Table 3 |
| network_kinase_analysis/ | decoupleR transcription-factor activity inference; kinase-network expression chart | Section 3.7.2-3.7.3, Fig. 7D |
| figures/ | Supplementary heatmap and slide-figure generation scripts | Supplementary figures |
| diagnostics/ | Gene-count reconciliation, software version checks, sample verification | Methods (Bioinformatics subsection), response to reviewer comments |
| utilities/ | Metadata generation, sample-accession checking, count-file merging | Methods |

## Known Exclusions

- E2F1_target_heatmap_figure.R and top12_heatmap_only.R were exploratory and do not correspond to any figure in the final manuscript; not deposited.
- WGCNA-related scripts were part of an earlier, separate analytical approach and are not part of the current manuscript; not deposited.
- Per-sample/per-contrast duplicate scripts (e.g., repeated GSEA or Picard submissions differing only in sample ID) are represented here by a single, labeled "representative" script; the same command was applied once per sample/contrast across the full cohort.

## Note on File Naming

One script in the raw-processing pipeline (Picard_97.sh, not included in this repository) was found during preparation to contain an internal reference to the wrong sample ID (SRR24071896 instead of SRR24071897), apparently from an incomplete copy-paste edit. This was verified against the final counts matrix: SRR24071896 and SRR24071897 have distinct, non-identical count data (Pearson r = 0.95, consistent with two genuine biological samples rather than a duplicated run), confirming Picard_97.sh was an abandoned draft that was never executed under this name and did not affect the reported results.

## Software Environment

Primary differential-expression analysis: R v4.3.1, DESeq2 v1.40.2.
Sex x condition interaction testing and downstream pathway analyses: R v4.5.2, DESeq2 v1.50.2, limma v3.66.0, fgsea v1.36.2, msigdbr v26.1.1, decoupleR. Full sessionInfo() output is provided in diagnostics/sessionInfo_output.txt.

## Citation

If you use these scripts, please cite the associated manuscript (citation to be finalized upon publication).
