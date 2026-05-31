# NGS-analysis-scripts-
This repository contains the scripts for my Master's capstone project titled "Exploring alternative splicing isoforms in Ovarian Cancer" using the dataset GSE102073 taken from GEO and published as "Molecular analysis of high-grade serous ovarian carcinoma with and without associated serous tubal intra-epithelial carcinoma. Ducie J, Dao F, Considine M, Olvera N et al. Nat Commun 2017 Oct 17;8(1):990."
This paper was selected on for the sequence data having a high read depth coverage of 40 million reads, necessary for an analysis for alternative splicing.

It contains the following:

salmon_download.sh - a bash script to download the fastq files for this project from the Sequence Read Archive 

salmon_quant.sh - a script to transcript gene quantification using the pseudoaligner Salmon against the indexed GENCODE Human transcriptome

star_aligner.sh - a bash script to align reads to the current release of the GENCODE Human genome set using STAR

dexseq_count.py - a bash script that carries out the DEXSEQ differential splicing 

Ovarian_cancer_project.rmd - an R Markdown document that documents the importing of the transcriptome data into R to be used with DESeq2 for differential expression of transcripts, differential splicing to be analyzed with DEXSEQ plus isoform switching analysis with IsoformSwitchAnalyzeR, Gene Set Enrichment Analysis (GSEA), and WGCNA analysis to be later exported into Cytoscape.  

