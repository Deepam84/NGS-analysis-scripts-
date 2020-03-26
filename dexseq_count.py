#!/usr/bin/python

#A script that runs samtools on unsorted bams in a given directory then quantifies
#The reads using the HTSeq protocol found in the Bioconducotr package DEXSeq

import os, sys #import modules
import re

inputdir = "/home/deepam84/Bioinformatics/RBIF120/star_analysis/aligned/"
R_location = "/mnt/d/Users/deepa/OneDrive/Documents/R/win-library/3.6/DEXSeq/python_scripts/dexseq_count.py"
Index_file = "/home/deepam84/Bioinformatics/RBIF120/gencode.v33.annotation.gff"
outputdir = "/home/deepam84/Bioinformatics/RBIF120/star_analysis/counted/"


for bam in os.listdir(inputdir): #loop through unsorted bams directory files
    print(bam)
    if bam.endswith("Aligned.out.bam"):
        name = re.sub(r"Aligned.out.bam","", bam ,0, re.MULTILINE)   #strip out suffix 
        #name_sorted = name + "sorted.bam"
        print(name)
        samtools_command = "samtools sort -n " + inputdir + bam + " " + inputdir + name + ".sorted.bam" # samtools command to sort reads by name
        print(samtools_command) 
        os.system(samtools_command) #run samtools command 

for bam in os.listdir(inputdir): #loop through sorted bams in directory
    if bam.endswith(".sorted.bam.bam"):
        print(bam)
        name = re.sub(r".sorted.bam.bam","",bam,0,re.MULTILINE) #strip out suffix
        DexSeq_count = "python3" + " " + R_location + " " + Index_file + " " + inputdir + bam + " " + outputdir + name + ".txt -p yes -f bam" #DEXSeq counting
        print(DexSeq_count)
        os.system(DexSeq_count) #run DEXSeq/HTSeq counting code located in the R folder

