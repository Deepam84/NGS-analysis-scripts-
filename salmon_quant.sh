#!/bin/bash

#A script to run the salmon transcript/gene quantification step on fastq read pairs

data="ovarian_sequence_data" #specify directory
#proj_acc_name="SRR5884"


for fn in ${data}/SRR5884{094..178}; #loop through each sample run filename
do
#echo $fn
samp=`basename ${fn}` 
echo "Processing sample ${samp}" #print sample name
salmon quant -i gencode_index -l A \ #quantify transcripts/genes with gencode index
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant #output folder with counts for each sample
done 

#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz