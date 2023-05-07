library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

#Get project information
GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

#Get specific project information
TCGAbiolinks:::getProjectSummary("TCGA-BRCA")

#Download TCGA Expression Data
query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts")

#Visualize Results
BRCA_res = getResults(query_TCGA) # make results as table
head(BRCA_res) # data of the first 6 patients.
colnames(BRCA_res) # columns present in the table

#Check sample information
head(BRCA_res$sample_type) # first 6 types of tissue.

#Total sample type of each category
summary(factor(BRCA_res$sample_type)) # summary of distinct tissues types present in this study

[There are 100 controls  and 742 cancer samples. For simplicity, we will ignore the small class of recurrent solid tumors]

barcodes<- c("TCGA-BH-A0B6",
             "TCGA-GM-A2D9",
             "TCGA-BH-A0HA",
             "TCGA-BH-A18P",
             "TCGA-AR-A2LE",
             "TCGA-BH-A0W7",
             "TCGA-E2-A1IN",
             "TCGA-E2-A1LH",
             "TCGA-GM-A2DH",
             "TCGA-GM-A2DO",
             "TCGA-A8-A095",
             "TCGA-BH-A0BP",
             "TCGA-E2-A152",
             "TCGA-BH-A0H5",
             "TCGA-BH-A209",
             "TCGA-BH-A0WA",
             "TCGA-BH-A0BW",
             "TCGA-Z7-A8R6",
             "TCGA-BH-A0B8",
             "TCGA-E2-A14Z",
             "TCGA-A8-A06O",
             "TCGA-BH-A0AV",
             "TCGA-E2-A1IH",
             "TCGA-AR-A255",
             "TCGA-BH-A0B0",
             "TCGA-BH-A0HY",
             "TCGA-OL-A66J",
             "TCGA-AR-A2LR",
             "TCGA-AR-A1AK",
             "TCGA-AR-A1AY",
             "TCGA-GM-A2DI",
             "TCGA-GM-A2DD",
             "TCGA-B6-A40B",
             "TCGA-BH-A0BR",
             "TCGA-BH-A0DX",
             "TCGA-BH-A0BL",
             "TCGA-BH-A18K",
             "TCGA-AO-A03V",
             "TCGA-E2-A154",
             "TCGA-A2-A3XZ",
             "TCGA-AR-A24S",
             "TCGA-GM-A2DL",
             "TCGA-AR-A1AJ",
             "TCGA-EW-A1J6",
             "TCGA-BH-A1FD",
             "TCGA-A2-A0YF",
             "TCGA-E2-A15C",
             "TCGA-EW-A1IY",
             "TCGA-A7-A0CD",
             "TCGA-AR-A1AP",
             "TCGA-AR-A24N",
             "TCGA-A8-A0AD",
             "TCGA-E2-A15O",
             "TCGA-BH-A0BG",
             "TCGA-S3-AA14",
             "TCGA-E2-A1II",
             "TCGA-GM-A2DK",
             "TCGA-BH-A0DO",
             "TCGA-AR-A1AX",
             "TCGA-BH-A1FG",
             "TCGA-A1-A0SE",
             "TCGA-BH-A18S",
             "TCGA-BH-A0H6",
             "TCGA-BH-A0BO",
             "TCGA-AR-A24P",
             "TCGA-BH-A1ET",
             "TCGA-BH-A1EU",
             "TCGA-E2-A156",
             "TCGA-B6-A1KI",
             "TCGA-B6-A0X0",
             "TCGA-E2-A15J",
             "TCGA-BH-A0BQ",
             "TCGA-A2-A259",
             "TCGA-E2-A1IF",
             "TCGA-A2-A0EP",
             "TCGA-BH-A0C3",
             "TCGA-E2-A1IJ",
             "TCGA-BH-A0H3",
             "TCGA-E2-A1IO",
             "TCGA-A2-A0YI",
             "TCGA-AR-A252",
             "TCGA-A8-A08A",
             "TCGA-B6-A402",
             "TCGA-AO-A03M",
             "TCGA-E2-A14U",
             "TCGA-A7-A3IY",
             "TCGA-AO-A03U",
             "TCGA-E2-A15F",
             "TCGA-E2-A14S",
             "TCGA-A1-A0SB",
             'TCGA-BH-A18V',
             'TCGA-BH-A1FJ',
             'TCGA-E2-A1LH',
             'TCGA-BH-A1EW', 
             'TCGA-BH-A0BQ',
             'TCGA-E9-A1N6',
             'TCGA-AC-A2FB',
             'TCGA-BH-A0DQ',
             'TCGA-BH-A0B5',
             'TCGA-E9-A1R7',
             'TCGA-AC-A23H',
             'TCGA-BH-A1FG',
             'TCGA-BH-A0DO',
             'TCGA-BH-A1FD',
             'TCGA-E9-A1RB',
             'TCGA-E9-A1RF',
             'TCGA-BH-A1FR',
             'TCGA-E9-A1NG',
             'TCGA-A7-A13F',
             'TCGA-BH-A1EV',
             'TCGA-BH-A0E0',
             'TCGA-BH-A0B3',
             'TCGA-E9-A1N5',
             'TCGA-BH-A18K',
             'TCGA-BH-A18S',
             'TCGA-BH-A204')
#Redo Design
query.exp.BRCAi <- GDCquery(project = "TCGA-BRCA",
                            data.category = "Gene expression",
                            data.type = "Gene expression quantification",
                            platform = "Illumina HiSeq", 
                            file.type  = "results", 
                            sample.type = c("Primary Tumor", "Solid Tissue Normal"),
                            legacy = TRUE,
                            barcode=barcodes)


setwd('D:/CancerData')

#Download TCGA data
GDCdownload(query.exp.BRCAi)

#Load RNASeq
tcga_data = GDCprepare(query.exp.BRCAi)
dim(tcga_data)

# chain functions to save time and space
colnames(colData(tcga_data))
metadata<- as.data.frame(rowData(tcga_data))
coldata<- as.data.frame(colData(tcga_data))


#Check Different features
table(tcga_data@colData$vital_status)
table(tcga_data@colData$tumor_stage)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$age_at_index)


#Obtain Count Matrix
dim(assay(tcga_data))
head(assay(tcga_data)[,1:10])
head(rowData(tcga_data))

#Save and Load 
saveRDS(object = tcga_data,
        file = "TCGA_BRCA.RDS",
        compress = FALSE)

tcga_data = readRDS(file = "TCGA_BRCA.RDS")



#DEG Analysis
#Create group
clinical_data = colData(tcga_data)
group = factor(clinical_data$definition)


#Assign reference group
group = relevel(group, ref="Solid Tissue Normal")

#Form design matrix
design = model.matrix(~group)
head(design)

#Create DEGlist object for LIMMA-VOOM
dge = DGEList( counts=assay(tcga_data),
  samples=colData(tcga_data),
  genes=as.data.frame(rowData(tcga_data)))


# filtering
keep = filterByExpr(dge,design) # defining which genes to keep
dge = dge[keep,,keep.lib.sizes=FALSE] # filtering the dge object
rm(keep) #  use rm() to remove objects from memory if you don't need them anymore

#TMM (trimmed mean of M-values) Normalization 
dge = calcNormFactors(dge,method="TMM")
v = voom(dge,design,plot=TRUE)


#Create MDS Plot
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(dge, labels = group, col=col.group)


#Creating boxplot of sample distribution
dge2 <- dge
dge2$samples$norm.factors <- 1
dge2$counts[,1] <- ceiling(dge2$counts[,1]*0.05)
dge2$counts[,2] <- dge2$counts[,2]*5

par(mfrow=c(1,2))
lcpm <- cpm(dge2, log=TRUE)
boxplot(lcpm, las=2, col='red', main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

dge2 <- calcNormFactors(dge2)  
dge2$samples$norm.factors

lcpm <- cpm(dge2, log=TRUE)
boxplot(lcpm, las=2, col='blue', main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

#Applying Fit model
fit = lmFit(v, design)
fit = eBayes(fit)
topGenes = topTable(fit, coef=1, sort.by="p")
print(topGenes)

DEGs<- topTable(fit, coef = 2, n=Inf, adjust="BH")
write.csv(DEGs, file="BRCA_DEGs.csv")

# Creating Limma Pipeline
limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number=Inf, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}

#Limma res file
limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)

# Save the data as a file, if you need it later, you can just load this file
# instead of having to run the whole pipeline again
saveRDS(object = limma_res,
        file = "limma_res.RDS",
        compress = FALSE)

#Make PCA Plot
plot_PCA = function(voomObj, condition_variable){
  group = factor(voomObj$targets[, condition_variable])
  pca = prcomp(t(voomObj$E))
  # Take PC1 and PC2 for the plot
  plot(pca$x[,1:2],col=group, pch=19)
  # include a legend for points
  legend("bottomleft", inset=.01, levels(group), pch=19, col=1:length(levels(group)))
  return(pca)
}

#Visualize PCA Plot
res_pca = plot_PCA(limma_res$voomObj, "definition")


