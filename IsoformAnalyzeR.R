
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("IsoformSwitchAnalyzeR")

BiocManager::install("Biostrings")

library("IsoformSwitchAnalyzeR")
library("Biostrings")


samples <- read.csv(file.path("D:","Users","deepa","Downloads", "SraRunTable.txt"), header = TRUE)
samples
samples_run <- samples$Run
samples_run <- as.vector(samples_run)
condition_stic <- as.character(samples$stic)
condition_plat <- as.character(samples$platinumstatus)

dir2 <- "D:/Users/deepa/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/deepam84/Bioinformatics/RBIF120/salmon_analysis"
dir5 <- "D:/Users/deepa/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/deepam84/Bioinformatics/RBIF120/star_analysis"
dir6 <- "D:/Users/deepa/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/deepam84/Bioinformatics/RBIF120/salmon_analysis/quants"

salmonQuant <- importIsoformExpression(
  parentDir = (dir6),
  addIsofomIdAsColumn = TRUE)

# colnames(salmonQuant$abundance)[-1]  <- samples_run
# colnames(salmonQuant$counts)[-1]  <- samples_run


myDesign_stic <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = c("Yes", "No",  "No",  "Yes", "No",  "Yes", "No",  "No",  "No",  "No",  "Yes" ,"No" , "No",  "No" , "No",  "No",  "No",  "Yes", "Yes", "Yes", "Yes", "Yes", "Yes", "No",  "Yes", "Yes", "Yes",
   "Yes", "No",  "No",  "No",  "Yes", "Yes", "Yes", "Yes", "No" , "Yes", "No" , "Yes", "Yes", "No",  "No" , "Yes", "Yes", "No",  "No",  "No",  "Yes", "No",  "No" , "No",  "No",  "No",  "No", 
   "Yes", "Yes", "No",  "No",  "Yes", "Yes", "No",  "No",  "Yes", "Yes", "Yes", "Yes", "Yes", "No",  "Yes", "No",  "Yes", "Yes", "No",  "Yes", "Yes", "Yes", "Yes", "Yes", "No",  "No",  "Yes",
   "Yes", "No",  "No",  "No")
)

myDesign_plat <- data.frame(
  sampleID = colnames(salmonQuant$abundance)[-1],
  condition = c( "Not available", "Not available", "Not available", "Not available", "Not available", "Not available", "Not available", "Not available", "Not available",
                 "Not available", "Not available", "Not available", "Not available", "Sensitive",     "Not available", "Sensitive",     "Sensitive",     "Sensitive" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Resistant",     "Sensitive",     "Resistant" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Resistant",     "Sensitive",     "Resistant",     "Not available", "Sensitive" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Resistant",     "Sensitive",     "Sensitive",     "Sensitive" ,   
                 "Resistant",     "Sensitive",     "Sensitive",     "Resistant",     "Sensitive",     "Resistant",     "Resistant",     "Sensitive",     "Resistant" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Resistant",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive" ,   
                 "Sensitive",     "Sensitive",     "Sensitive",     "Sensitive" )
)

# myDesign2 
# 
# head(salmonQuant$abundance, 2)

#moo <- importGTF(file.path(dir5,"gencode.v33.annotation.gtf"))

aSwitchList_stic <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign_stic,
  isoformExonAnnoation = file.path(dir5,"gencode.v33.annotation.gtf"),
  isoformNtFasta       = file.path(dir2,"gencode.fa.gz"),
  showProgress = FALSE
)

aSwitchList2_plat <- importRdata(
  isoformCountMatrix   = salmonQuant$counts,
  isoformRepExpression = salmonQuant$abundance,
  designMatrix         = myDesign_plat,
  isoformExonAnnoation = file.path(dir5,"gencode.v33.annotation.gtf"),
  isoformNtFasta       = file.path(dir2,"gencode.fa.gz"),
  showProgress = FALSE
)

# aSwitchList
# aSwitchList2_plat

aSwitchListFiltered_stic <- preFilter(
  switchAnalyzeRlist = aSwitchList_stic,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

aSwitchListFiltered_plat <- preFilter(
  switchAnalyzeRlist = aSwitchList2_plat,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE
)

aSwitchListFiltered_stic <- preFilter(aSwitchListFiltered_stic) # preFilter to remove lowly expressed features
aSwitchListFiltered_plat <- preFilter(aSwitchListFiltered_plat) # preFilter to remove lowly expressed features


aSwitchListAnalyzed_plat <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered_plat,
  reduceToSwitchingGenes=TRUE
)

aSwitchListAnalyzed_stic <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchListFiltered_plat,
  reduceToSwitchingGenes=TRUE
)

switch_analysis_stic <- isoformSwitchAnalysisPart1(aSwitchListFiltered_stic)
switch_analysis_plat <- isoformSwitchAnalysisPart1(aSwitchListFiltered_plat)

extractSwitchSummary(switch_analysis_plat)
# 
# aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
#   switchAnalyzeRlist = aSwitchListFiltered,
#   reduceToSwitchingGenes=TRUE
#)
# switch_analysis <- isoformSwitchAnalysisPart1(aSwitchListFiltered,switchTestMethod='DRIMSeq')
# switch_analysis_alt <- isoformSwitchAnalysisPart1(aSwitchListFiltered)
# switch_analysis_alt2 <- isoformSwitchAnalysisPart1(aSwitchListFiltered2)
# 
# switch_analysis_stic <- isoformSwitchAnalysisCombined(aSwitchListFiltered_stic,switchTestMethod='DRIMSeq')
# switch_analysis_plat <- isoformSwitchAnalysisCombined(aSwitchListFiltered_plat,switchTestMethod='DRIMSeq')
# extractSwitchSummary( exampleSwitchList )

# exampleSwitchListAnalyzed <- analyzeORF(
#   exampleSwitchListAnalyzed,
#   orfMethod = "longest",
#   # genomeObject = Hsapiens, # not necessary since sequences are already stored in the switchAnalyzeRlist
#   showProgress=FALSE
# )
# 
# getCDS()
# 
# exampleSwitchListAnalyzed <- analyzeAlternativeSplicing(
#   switchAnalyzeRlist = exampleSwitchListAnalyzed,
#   quiet=TRUE
# )
# 
# data("switch_analysis_stic")
# switch_analysis_stic
# 
# extractSwitchSummary(
#   switch_analysis_stic,
#   filterForConsequences = TRUE
# ) 
