# NGS-analysis-scripts-
This repository contains the scripts for my Master's capstone project titled "Exploring alternative splicing isoforms in Ovarian Cancer" using the dataset GSE102073 taken from GEO. The associated paper is titled "Molecular analysis of high-grade serous ovarian carcinoma with and without associated serous tubal intra-epithelial carcinoma. Ducie J, Dao F, Considine M, Olvera N et al. Nat Commun 2017 Oct 17;8(1):990."
This paper was selected because the sequence data had a depth coverage of 40 million reads, necessary to analyze for differences in alternative splicing of transcripts. The comparison groups were Platinum Sensitive and Resistant Ovarian Cancer and Ovarian Cancer with and without Serous Tubal Epithelial Carcinoma lesions (STICs). The following analyses were carried out: differential expression of transcripts, differential splicing events and isoform switching.

It contains the following:

salmon_download.sh - a bash script to download the paired fastq files of the raw sequence data from the Sequence Read Archive (SRA) 

salmon_quant.sh - a script to quantify transcript isoforms using the pseudoaligner Salmon against the indexed GENCODE Human transcriptome

star_aligner.sh - a bash script to align the paired fastq reads to GENCODE Human reference genome using STAR and outputs unsorted bam files

dexseq_count.py - a Python script that runs Samtools on the unsorted bam files and then quantifies them using HTSeq

Ovarian_cancer_project.rmd - an R Markdown document that documents the following: 

The importing of the transcript isoform counts into R to be used with DESeq2 for analyzing differential expression of transcripts

The DESeq2 analysis with Volcano plots and hierarchical clustering 

The importing of the counted reads to be analyzed with DEXSEQ for analyzing differential splicing by testing for differential exon usage

Isoform switching analysis with IsoformSwitchAnalyzeR using the Salmon transcript counts, DEXSEQ and DRIMSeq

Gene Set Enrichment Analysis (GSEA), GO.db analysis and DAVID analysis 

Whole Genome Network Correlation (WGCNA) generation of network data to be exported into Cytoscape


