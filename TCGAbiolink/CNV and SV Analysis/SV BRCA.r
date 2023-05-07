library(TCGAbiolinks)
library(TCGAmutations)
library(maftools)

setwd('E:/TCGA Data Analysis/CNV and SNV Analysis')

#Create TCGA object and read maf data
tcga_brca=tcga_load(study = "BRCA") 
BRCA = read.maf(maf = tcga_brca@data, clinicalData = tcga_brca@clinical.data)


#Take a overview of the data
plotmafSummary(maf = BRCA, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = T)

#Make oncoplot of top 10 mutated genes in BRCA
oncoplot(maf = BRCA, top = 10)

#Let's make an oncoprint of BRCA1 and BRCA2 genes
oncostrip(maf = BRCA, genes = c('BRCA1','BRCA2'))
dev.off()

#Let's add mutation type strip
library(grid)
grid.rect(gp=gpar(fill=0),vp=viewport(layout.pos.col=500,layout.pos.row=300))
oncostrip(maf = BRCA, genes = c('BRCA1','BRCA2'),clinicalFeatures = 'stage_event_pathologic_stage',width(1000))

#Let's check the total transition and transersion mutation frequency
titv = titv(maf = BRCA, plot = T)


#lollipop plot for BRCA1 gene in breast cancer.
lollipopPlot(maf = BRCA, gene = 'BRCA1', AACol = 'HGVSp_Short', showMutationRate = TRUE)

#Compare mutation load with other TCGA cohort
BRCA.mutload = tcgaCompare(maf = BRCA, cohortName = 'BRCA')

#Visualizing top co-occuring mutated genes
somaticInteractions(maf = BRCA, top = 25, pvalue = c(0.05, 0.1))

#Detecting cancer driver genes that are frequently mutated
BRCA.sig = oncodrive(maf = BRCA, AACol = 'HGVSp_Short', minMut = 5, pvalMethod = 'zscore')
head(BRCA.sig)
plotOncodrive(BRCA.sig)

#Visualizing multiple genes
multi_genes = c("TP53", "MCM2", "BRCA1", "BRCA2", "DNMT3B", "PTEN", "STK11", "IDH1", "IDH2", "FLT3")
oncoplot(maf = BRCA, genes = multi_genes)

