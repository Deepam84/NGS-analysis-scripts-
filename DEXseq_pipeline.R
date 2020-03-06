if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DEXSeq")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library("DEXSeq")
library("Biostrings")

pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)
system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

dir3 <- "D:/Users/deepa/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/deepam84/Bioinformatics/RBIF120/star_analysis/counted"

dir4 <- "D:/Users/deepa/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/deepam84/Bioinformatics/RBIF120"

countFiles = list.files(dir3, pattern=".txt$", full.names=TRUE)
basename(countFiles)

flattenedFile = list.files(dir4, pattern="gff$", full.names=TRUE)
basename(flattenedFile)

samples <- read.csv(file.path("D:","Users","deepa","Downloads", "SraRunTable.txt"), header = TRUE)
samples

sampleTable = data.frame(samples[c(1,18,23,30)])

sampleTable[sampleTable$platinumstatus == 'Resistant' ,]
sampleTable[sampleTable$platinumstatus == 'Sensitive' ,]

sampleTable[sampleTable$stic == 'Yes' ,]
sampleTable[sampleTable$stic == 'No' ,]

selected_stic_samples <- data.frame(sampleTable[c(18,28,33,48,55,56,66,75,81,82),])
countFiles_plat = countFiles[c(18,19,20,21,25,27,60,67,72,85)]

basename((countFiles_plat))
library("DEXSeq")

dxd1 = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ Run + exon + stic:exon,
  flattenedfile=flattenedFile )

dxd2 = DEXSeqDataSetFromHTSeq(
  countFiles_plat,
  sampleData=selected_plat_samples,
  design= ~ Run + exon + platinumstatus:exon,
  flattenedfile=flattenedFile )

colData(dxd1)
colData(dxd2)

head( counts(dxd1), 5 )
head( counts(dxd2), 5 )

head( featureCounts(dxd1), 5 )
head( featureCounts(dxd2), 5 )

head( rowRanges(dxd1), 3 )
head( rowRanges(dxd2), 3 )


dxd1 = estimateSizeFactors(dxd1) # 1 = stic 
dxd1t = estimateSizeFactors(dxd1t) # 1 = stic 
dxd2 = estimateSizeFactors(dxd2 ) # 2 = platinum

saveRDS(dxd1,"dexseq_stic.rds")
saveRDS(dxd2,"dexseq_platinum.rds")

dxd1 <- dxd1[rowSums(counts(dxd1)) > 0,]
dxd1t <- dxd1t[rowSums(counts(dxd1t)) > 0,]

dxd2 <- dxd2[rowSums(counts(dxd2)) > 0,]

dxd1 <- dxd1[idx1,]
dxd2 <- dxd2[idx2,]

dim( counts(dxd1) )
dim( counts(dxd1t) )
dim( counts(dxd2) )

dxd1 = estimateDispersions( readRDS("dexseq_stic.rds"))
dxd1t = estimateDispersions( readRDS("dexseq_stic_test.rds"))
dxd2 = estimateDispersions( readRDS("dexseq_platinum.rds"),quiet=FALSE,BPPARAM=SnowParam(3))

dxd2 = estimateDispersions( dxd2 )

dxd1 = estimateDispersionsTRT( dxd1 )

plotDispEsts( dxd1 )
plotDispEsts( dxd2 )

estimateDispersionsTRT
estimateDispersionsTRT

dxd1 = testForDEU( dxd1 )
dxd2 = testForDEU( dxd2 )

dxd2 = estimateExonFoldChanges( dxd2, fitExpToVar="condition")

dxr2 = DEXSeqResults( dxd2 )
dxd1 = estimateExonFoldChanges( dxd1, fitExpToVar="stic")
dxd2 = estimateExonFoldChanges( dxd2, fitExpToVar="platinumstatus")

dxr1 = DEXSeqResults( dxd1 )
dxr1
dxr2 = DEXSeqResults( dxd2 )
dxr2

table ( dxr2$padj < 0.1 )
plotMA( dxr2, cex=0.8 )
table ( tapply( dxr2$padj < 0.1, dxr2$groupID, any ) )

dxr2Sig <- subset(dxr2, padj < 0.1) # strongest downregulation 
head(dxr2Sig[ order(dxr2Sig$log2fold_Sensitive_Resistant), ])

head(dxr2Sig [ order(dxr2Sig$log2fold_Sensitive_Resistant, decreasing=TRUE), ]) # strongest upregulation 

resOrdereddxr2 <- dxr2[order(dxr2$padj),]
head(resOrdereddxr2)

plotDEXSeq( dxr2, "ENSG00000026025", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

resSig2 <- subset(res2s, padj < 0.1)
head(resSig2[ order(resSig2$log2FoldChange), ])

resSig3 <- subset(res3p, padj < 0.1)
head(resSig3[ order(resSig3$log2FoldChange), ])

head(resSig1[ order(resSig1$log2FoldChange, decreasing=TRUE), ]) # strongest upregulation 
head(resSig2[ order(resSig2$log2FoldChange, decreasing=TRUE), ])
head(resSig3[ order(resSig3$log2FoldChange, decreasing=TRUE), ])

resOrdered1 <- res1ps[order(res1ps$padj),]
head(resOrdered1)

resOrdered2 <- res2s[order(res2s$padj),]
head(resOrdered2)

resOrdered3 <- res3p [order(res3p$padj),]
head(resOrdered3)

plotDEXSeq( dxr2, "FBgn0010909", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
