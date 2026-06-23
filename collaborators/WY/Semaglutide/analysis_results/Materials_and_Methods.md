# Materials and Methods

## 16S rRNA Amplicon Sequencing and Analysis

### Sample Processing and Sequencing

Microbial DNA was extracted from lung, gut (intestinal/BO), oral, and fecal samples using a commercial DNA extraction kit following the manufacturer's instructions. The bacterial 16S rRNA gene V3–V4 hypervariable region was amplified using primers 341F/806R, and paired-end sequencing (2 × 250 bp) was performed on an Illumina MiSeq platform.

### Bioinformatic Processing

Raw sequencing reads were processed using the DADA2 pipeline (v1.x). Briefly, adapter sequences were removed using Cutadapt, followed by quality filtering and truncation (filterAndTrim). Error rates were learned independently for each sample (learnErrors), and amplicon sequence variants (ASVs) were inferred via the dada function. Paired-end reads were merged (mergePairs) and chimeric sequences were removed (removeBimeraDenovo). Taxonomic classification was performed using the SILVA 138 reference database (assignTaxonomy), annotating ASVs to the genus level. The resulting genus-level abundance table comprised 911 genera across 47 samples.

### Diversity Analysis

Alpha diversity metrics (Observed ASVs, Shannon index, and Simpson index) were calculated from rarefied ASV count tables. Beta diversity was assessed using Bray–Curtis dissimilarity matrices and visualized by principal coordinates analysis (PCoA). Differences in community structure among groups were tested using PERMANOVA (vegan::adonis2, 999 permutations).

---

## Dysbiosis Index Analysis

Microbial community dysbiosis was quantified using the dysbiosisR R package. The WT_Control group was designated as the healthy reference community. For each sample, the median Aitchison distance to all reference samples was calculated as the MedianCLV (Median Centered Log-Ratio Variance) dysbiosis index, where higher values indicate greater deviation from the healthy community composition.

---

## Transcriptomic Sequencing and Analysis

### RNA Extraction, Library Preparation, and Sequencing

Total RNA was extracted from lung, gut, and oral tissue samples using the TRIzol reagent (Thermo Fisher Scientific) according to the manufacturer's protocol. RNA integrity was verified by Bioanalyzer (RIN ≥ 7). Strand-specific RNA-seq libraries were constructed following poly(A) selection or ribosomal RNA depletion, and sequenced on an Illumina platform (paired-end, 150 bp). The experimental design comprised four groups (WT_Control, WT_Treated, DB_Control, and DB_Treated) with three biological replicates per group per tissue site, yielding 36 samples in total (3 sites × 4 groups × 3 replicates).

### Read Alignment and Quantification

Raw FASTQ reads were quality-trimmed using Trim Galore. Trimmed reads were aligned to the mouse reference genome (GRCm38/mm10) using HISAT2 (v2.x) with default splice-aware alignment parameters. Gene-level read counts were quantified using featureCounts (Subread v2.x) against the Ensembl genome annotation (GTF), generating a count matrix of 54,532 genes × 36 samples.

### Differential Expression Analysis

Raw count data were imported into DESeq2 (v1.x) in R and normalized using the median-of-ratios method. Variance-stabilizing transformation (VST) was applied for visualization purposes. Two complementary testing frameworks were employed:

**Multi-group testing (LRT)**: A likelihood ratio test (LRT) was performed comparing a full model (~Group) to a reduced intercept-only model (~1), identifying genes whose expression varied significantly across the four groups (group-variable genes; FDR < 0.05).

**Pairwise testing (Wald test)**: Wald tests were applied to four pre-specified contrasts: DB_Control vs. WT_Control (disease effect), DB_Treated vs. DB_Control (treatment effect), DB_Treated vs. WT_Control (residual deviation after treatment), and WT_Treated vs. WT_Control (drug effect in healthy animals). Differentially expressed genes (DEGs) were defined as FDR < 0.05 and |log2FoldChange| > 1. P-values were adjusted for multiple testing using the Benjamini–Hochberg (BH) procedure throughout.

### Expression Trend Clustering and Rescue Gene Identification

Among LRT-significant genes, VST-normalized mean expression values across three experimental states (WT_Control, DB_Control, and DB_Treated) were used as feature vectors for k-means clustering (k = 6). Rescue clusters were defined as those satisfying: (i) a substantial shift in DB_Control relative to WT_Control, and (ii) a reversal toward WT_Control values in DB_Treated, capturing the WT → DB → DB_Treated recovery trajectory.

A gene-level rescue score was additionally computed for each LRT-significant gene as:

$$\text{Rescue Score} = -\frac{\text{LFC}_{DB\_Treated/DB\_Control}}{\text{LFC}_{DB\_Control/WT\_Control}}$$

A rescue score of 1.0 indicates complete normalization, values between 0 and 1 indicate partial recovery, and values > 1 indicate over-correction. This score was used as the ranking metric for rescue-focused gene set enrichment analysis (GSEA).

### Functional Enrichment Analysis

DEGs from the DB_Control vs. WT_Control contrast (FDR < 0.05, |LFC| > 0.5) were subjected to over-representation analysis (ORA) for Gene Ontology biological processes (GO-BP) and KEGG pathways using clusterProfiler (v4.x) with BH correction. Gene set enrichment analysis (GSEA) was performed using GO-BP, KEGG, and MSigDB Hallmark gene sets, with genes pre-ranked by Wald test statistic. Pathways related to inflammatory response, epithelial barrier integrity, cytokine signaling, and oxidative stress were prioritized for interpretation. Mouse gene annotation was obtained from org.Mm.eg.db.

---

## Host–Microbiome Joint Analysis

### CLR Transformation and Rescue Genus Identification

Genus-level relative abundances were transformed using the centered log-ratio (CLR) transformation: CLR(x) = log(x / g(x)), where g(x) denotes the geometric mean of all genera within a sample, to remove the compositional constraint of relative abundance data. For each genus, three-group CLR means were computed, and a genus-level rescue score was defined analogously to the gene-level metric:

$$\text{Rescue Score}_{genus} = -\frac{\Delta_{treat}}{\Delta_{disease}} = -\frac{\overline{CLR}_{DB\_Treated} - \overline{CLR}_{DB\_Control}}{\overline{CLR}_{DB\_Control} - \overline{CLR}_{WT\_Control}}$$

Genera satisfying the following criteria were classified as rescue genera: (i) |Δ_disease| ≥ 0.3 (sufficient disease-associated shift in CLR abundance); (ii) rescue score ≥ 0.3 (at least 30% reversal by treatment); and (iii) treatment-associated change in the opposite direction to the disease-associated change.

### Host Module Eigengene Computation

For each rescue gene cluster, a module eigengene was calculated as the first principal component (PC1) scores derived from PCA of the VST expression matrix of member genes, representing the overall transcriptional activity of that module across samples. When the PC1 orientation was inconsistent with the expected disease direction, scores were multiplied by −1 to ensure biological interpretability. This approach is mathematically equivalent to the module eigengene concept used in weighted gene co-expression network analysis (WGCNA).

### Eigengene–Genus Correlation Analysis

For each tissue site, Spearman rank correlations were computed between host module eigengenes and genus CLR values across matched samples (Hmisc::rcorr). BH correction was applied simultaneously to all eigengene–genus pairs within each site (FDR < 0.05 as significance threshold).

### Gene-Level Host–Microbiome Correlation in Oral Samples

In the oral cohort, Spearman rank correlations were computed between individual rescue gene VST expression values (34 genes) and rescue genus CLR values (13 genera) across matched samples (n ≈ 11), yielding 442 gene–genus pairs. BH correction was applied across all pairs; gene–genus pairs with uncorrected P < 0.05 were reported as exploratory significant associations, with FDR < 0.20 as a supplementary lenient threshold given the limited sample size.

### Host–Microbiome Interaction Network Construction

A bipartite network was constructed using significant gene–genus pairs (P < 0.05) as edges, with host rescue genes and microbial rescue genera as two distinct node classes. Network analysis was performed using the igraph package, and visualizations were generated with ggraph. Node size was scaled by degree (number of significant connections), and edge width was proportional to the absolute Spearman correlation coefficient. Positive and negative correlations were distinguished by edge color and line type. A force-directed layout (Fruchterman–Reingold algorithm) was applied to reveal the overall network topology. A bipartite layout was additionally used to emphasize the host–microbiome interface.

### Three-State Network Analysis

To visualize the dynamic reorganization of the host–microbiome interaction network across disease and treatment states, the network topology (node positions and edges) derived above was held constant, and nodes were independently colored according to their normalized deviation from WT_Control for each state (WT_Control, DB_Control, and DB_Treated): delta_norm = (group mean − WT mean) / max|deviation|. A diverging blue–white–red color scale was applied (blue: below WT level; white: equivalent to WT; red: above WT level), enabling direct visual comparison of network state across three groups on a fixed background.

### Association Between Host Module Eigengene and Dysbiosis Index

Spearman rank correlation was computed between host module eigengenes and the MedianCLV dysbiosis index across individual oral samples to test whether host transcriptional module activity was associated with the degree of microbial community dysbiosis at the individual sample level. This analysis was intended to link the host–microbiome co-regulatory network to treatment-associated restoration of oral microbial homeostasis.

---

## Statistical Analysis and Software

All statistical analyses were conducted in R (v4.x). Key packages included DESeq2 (differential expression analysis), clusterProfiler (enrichment analysis), org.Mm.eg.db (gene annotation), DADA2 (16S amplicon processing), vegan (diversity analysis), dysbiosisR (dysbiosis index), Hmisc (correlation analysis), igraph and ggraph (network analysis and visualization), pheatmap (heatmap), and ggplot2 (statistical visualization). Unless stated otherwise, multiple testing correction was performed using the Benjamini–Hochberg (BH) procedure to control the false discovery rate (FDR).
