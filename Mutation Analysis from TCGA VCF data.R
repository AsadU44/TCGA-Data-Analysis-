library(TCGAmutations)
library(maftools)
library(ggplot2)

#Working Directory
setwd('D:/CancerData')

#Fetch Mutation and Clinical Data
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

#Plot Summary
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

ggsave('Mutation1.jpg', height = 16, width = 20, units = 'cm', bg='white')

#Generate OncoPlot
oncoplot(maf = laml, top=10)

oncostrip(maf = laml, genes = c('TP53','CDH1', 'FLG')) #Selected Genes

#Detect Cancer Driver Gene
LAML.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = LAML.sig, fdrCutOff = 0.1, useFraction = TRUE)

#Summary of Transition and Transversion
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)


#Lollipop Diagram for a single gene 
lollipopPlot(
  maf = laml,
  gene = 'TP53',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
)


#Survival analysis based on Expression of Single Mutatied Gene
mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', 
            Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", 
             Status = "Overall_Survival_Status") #Group of Genes



#Selected Multiple Gene Mutation Status
aml_genes = c("TP53", "WT1", "PHF6", "DNMT3A", "DNMT3B", "TET1", "TET2", "IDH1", "IDH2", "FLT3")

aml_genes_vaf = subsetMaf(maf = laml, genes = aml_genes, fields = "i_TumorVAF_WU", 
                          mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
colnames(aml_genes_vaf)[2] = "VAF"
head(aml_genes_vaf)  #Extract variant allele frequency

laml.mutsig = system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.mutsig = data.table::fread(input = laml.mutsig)[,.(gene, q)]
laml.mutsig[,q := -log10(q)] #transoform to log10
head(laml.mutsig) #Assess significant mutated genes


oncoplot(
  maf = laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  rightBarLims = c(0, 20)
)

library(TCGAmutations)
library(maftools)
library(ggplot2)

#Working Directory
setwd('D:/CancerData')

#Fetch Mutation and Clinical Data
laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools') 
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 
laml = read.maf(maf = laml.maf, clinicalData = laml.clin)

#Plot Summary
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

ggsave('Mutation1.jpg', height = 16, width = 20, units = 'cm', bg='white')

#Generate OncoPlot
oncoplot(maf = laml, top=10)

oncostrip(maf = laml, genes = c('TP53','CDH1', 'FLG')) #Selected Genes

#Detect Cancer Driver Gene
LAML.sig = oncodrive(maf = laml, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = LAML.sig, fdrCutOff = 0.1, useFraction = TRUE)

#Summary of Transition and Transversion
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)


#Lollipop Diagram for a single gene 

lollipopPlot(
  maf = laml,
  gene = 'TP53',
  AACol = 'Protein_Change',
  showMutationRate = TRUE,
  
)


#Survival analysis based on Expression of Single Mutatied Gene
mafSurvival(maf = laml, genes = 'DNMT3A', time = 'days_to_last_followup', 
            Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvGroup(maf = laml, geneSet = c("DNMT3A", "FLT3"), time = "days_to_last_followup", 
             Status = "Overall_Survival_Status") #Group of Genes



#Selected Multiple Gene Mutation Status
aml_genes = c("TP53", "WT1", "PHF6", "DNMT3A", "DNMT3B", "TET1", "TET2", "IDH1", "IDH2", "FLT3")

aml_genes_vaf = subsetMaf(maf = laml, genes = aml_genes, fields = "i_TumorVAF_WU", 
                          mafObj = FALSE)[,mean(i_TumorVAF_WU, na.rm = TRUE), Hugo_Symbol]
colnames(aml_genes_vaf)[2] = "VAF"
head(aml_genes_vaf)  #Extract variant allele frequency

laml.mutsig = system.file("extdata", "LAML_sig_genes.txt.gz", package = "maftools")
laml.mutsig = data.table::fread(input = laml.mutsig)[,.(gene, q)]
laml.mutsig[,q := -log10(q)] #transoform to log10
head(laml.mutsig) #Assess significant mutated genes


oncoplot(
  maf = laml,
  genes = aml_genes,
  leftBarData = aml_genes_vaf,
  leftBarLims = c(0, 100),
  rightBarData = laml.mutsig,
  rightBarLims = c(0, 20)
)

#Further Reading
https://www.bioconductor.org/packages/devel/bioc/vignettes/maftools/inst/doc/maftools.html#references