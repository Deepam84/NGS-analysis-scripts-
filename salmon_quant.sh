#!/bin/bash

data="ovarian_sequence_data"
#proj_acc_name="SRR5884"


for fn in ${data}/SRR5884{094..178};
do
#echo $fn
samp=`basename ${fn}`
echo "Processing sample ${samp}"
salmon quant -i gencode_index -l A \
         -1 ${fn}/${samp}_1.fastq.gz \
         -2 ${fn}/${samp}_2.fastq.gz \
         -p 8 --validateMappings -o quants/${samp}_quant
done 

#ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz