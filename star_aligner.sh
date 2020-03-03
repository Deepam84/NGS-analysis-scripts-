#!/bin/bash

file="/home/deepa/Bioinformatics/RBIF120/ovarian_sequence_data"
genomeDir="/home/deepa/Bioinformatics/RBIF120/index/"
genomeFile="/home/deepa/Bioinformatics/RBIF120/gencode.v33.annotation.gtf"
FastaFile="/home/deepa/Bioinformatics/RBIF120/GRCh38.primary_assembly.genome.fa"
outfile="/home/deepam/Bioinformatics/RBIF120"
Star_location="/home/deepa/miniconda3/bin"
aligned_location="/home/deepa/Bioinformatics/RBIF120/aligned"

#$Star_location/STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $genomeDir \
#    --genomeFastaFiles $FastaFile --sjdbGTFfile $genomeFile --sjdbOverhang 100 
    
for fn in ${file}/SRR5884113;
do
echo ${fn} "alignment in progress"
samp=`basename ${fn}`
#mkdir $aligned_location/$samp
# #echo ${samp}
# #echo $fn/${samp}_1.fastq.gz
    $Star_location/STAR --runThreadN 12 --runMode alignReads --genomeDir $genomeDir \
    --readFilesIn $fn/${samp}_1.fastq.gz $fn/${samp}_2.fastq.gz \
    --readFilesCommand zcat \
    --outFileNamePrefix $aligned_location/$samp \
    --outSAMtype BAM Unsorted
sleep .25
done







