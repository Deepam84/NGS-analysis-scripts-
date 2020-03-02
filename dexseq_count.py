#!/usr/bin/python

import os, sys
import re

inputdir = "/home/deepam84/Bioinformatics/RBIF120/star_analysis/aligned/"
R_location = "/mnt/d/Users/deepa/OneDrive/Documents/R/win-library/3.6/DEXSeq/python_scripts/dexseq_count.py"
Index_file = "/home/deepam84/Bioinformatics/RBIF120/gencode.v33.annotation.gff"
outputdir = "/home/deepam84/Bioinformatics/RBIF120/star_analysis/counted/"


for bam in os.listdir(inputdir):
    print(bam)
    if bam.endswith("Aligned.out.bam"):
        name = re.sub(r"Aligned.out.bam","", bam ,0, re.MULTILINE)
        #name_sorted = name + "sorted.bam"
        print(name)
        samtools_command = "samtools sort -n " + inputdir + bam + " " + inputdir + name + ".sorted.bam"
        print(samtools_command)
        os.system(samtools_command)

for bam in os.listdir(inputdir):
    if bam.endswith(".sorted.bam.bam"):
        print(bam)
        name = re.sub(r".sorted.bam.bam","",bam,0,re.MULTILINE)
        DexSeq_count = "python3" + " " + R_location + " " + Index_file + " " + inputdir + bam + " " + outputdir + name + ".txt -p yes -f bam"
        print(DexSeq_count)
        os.system(DexSeq_count)

