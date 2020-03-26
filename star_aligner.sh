#!/bin/bash

#A script to run the STAR aligner on fastq paired RNA-Seq reads

file="/home/deepa/Bioinformatics/RBIF120/ovarian_sequence_data"
genomeDir="/home/deepa/Bioinformatics/RBIF120/index/"
genomeFile="/home/deepa/Bioinformatics/RBIF120/gencode.v33.annotation.gtf"
FastaFile="/home/deepa/Bioinformatics/RBIF120/GRCh38.primary_assembly.genome.fa"
outfile="/home/deepam/Bioinformatics/RBIF120"
Star_location="/home/deepa/miniconda3/bin"
aligned_location="/home/deepa/Bioinformatics/RBIF120/aligned"

#Generate the genome using the specified FASTA file and .gtf file with STAR

$Star_location/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $genomeDir \
    --genomeFastaFiles $FastaFile --sjdbGTFfile $genomeFile --sjdbOverhang 100 
    
for fn in ${file}/SRR5884{094..178}; #for each folder in the directory
do
echo ${fn} "alignment in progress" #print statement
samp=`basename ${fn}` #using the sample name
#mkdir $aligned_location/$samp
# #echo ${samp}
# #echo $fn/${samp}_1.fastq.gz
    $Star_location/STAR --runThreadN 12 --runMode alignReads --genomeDir $genomeDir \ #run STAR aligner
    --readFilesIn $fn/${samp}_1.fastq.gz $fn/${samp}_2.fastq.gz \ #use paired files
    --readFilesCommand zcat \ # allows STAR to read gzipped FASTQ files
    --outFileNamePrefix $aligned_location/$samp \ #specify output directory
    --outSAMtype BAM Unsorted #output unsorted BAMS
sleep .25 
done







