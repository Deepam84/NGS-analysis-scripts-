#!/bin/bash

#script to download a set of  RNA-Seq paired read data 

#data variable specifies the name of file to be created
#sequence specifies the names of sequences found by SRA 

data="ovarian_sequence_data"
# subfolder="SRR588"
# proj_acc_name="SRR5884"
seq=`seq -w 094 178`
seq2=`seq -w 000 009`

# if the specified directory doesn't exist then create it
if [ ! -d "./$data" ]
then
    mkdir ./$data
fi

cd $data
#go through each accession and download the data


for i in $seq;
do
    if [ -d "./SRR5884${i}" ]
    then
        echo "folder exists, skipping" #skip folder if it already exists
    else
        echo "creating directory SRR5884${i}"
        mkdir SRR5884${i}; #make directory if it doesn't 
    fi
    cd SRR5884${i};
    for j in $seq2; do #search through each ftp repo download fastq files 
        if [  -f "SRR5884${i}_1.fastq.gz" ]; then
            echo "SRR5884${i}_1.fastq.gz exists....skipping" 
        else    
            echo "downloading....."
            wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/${j}/SRR5884${i}/SRR5884${i}_1.fastq.gz" #download fastq_1 if found
        fi
        if [  -f "SRR5884${i}_2.fastq.gz" ]; then
            echo "SRR5884${i}_2.fastq.gz exists....skipping" 
        else
            echo "downloading......"
            wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/${j}/SRR5884${i}/SRR5884${i}_2.fastq.gz" #download fastq_2 if found 
        fi
    sleep .25

    done
    cd ..; 
    sleep .25
done
cd ..
