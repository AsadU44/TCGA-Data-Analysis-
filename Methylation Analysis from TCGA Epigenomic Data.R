library(MethylMix)
library(doParallel)

setwd('D:/CancerData')
cancerSite <- "BRCA"
targetDirectory <- paste0(getwd(), "/")
GetData(cancerSite, targetDirectory)

#Download methylation data
METdirectories <- Download_DNAmethylation(cancerSite, targetDirectory)

# Process
METProcessedData <- Preprocess_DNAmethylation(cancerSite, METdirectories)

#Save methylation processed data
saveRDS(METProcessedData, file = paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))

# Downloading gene expression data
GEdirectories <- Download_GeneExpression(cancerSite, targetDirectory)

# Processing gene expression data
GEProcessedData <- Preprocess_GeneExpression(cancerSite, GEdirectories)

# Saving gene expression processed data
saveRDS(GEProcessedData, file = paste0(targetDirectory, "GE_", cancerSite, "_Processed.rds"))

# Clustering probes to genes methylation data
METProcessedData <- readRDS(paste0(targetDirectory, "MET_", cancerSite, "_Processed.rds"))
res <- ClusterProbes(METProcessedData[[1]], METProcessedData[[2]])

# Putting everything together in one file
toSave <- list(METcancer = res[[1]], METnormal = res[[2]], GEcancer = GEProcessedData[[1]],
               GEnormal = GEProcessedData[[2]], ProbeMapping = res$ProbeMapping)
saveRDS(toSave, file = paste0(targetDirectory, "data_", cancerSite, ".rds"))

stopCluster(cl)

#Design
MethylMixResults$MethylationDrivers
MethylMixResults$MixtureStates

# Plot the most famous methylated gene for glioblastoma
plots <- MethylMix_PlotModel("MGMT", MethylMixResults, METcancer)
plots$MixtureModelPlot