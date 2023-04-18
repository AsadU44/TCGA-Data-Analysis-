# Load packages
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

#Set Directory
setwd('D:/CancerData')

#Read Data
tcga_data = readRDS(file = "tcga_data.RDS")
limma_res = readRDS(file = "limma_res.RDS")

# extract clinical data
clinical = tcga_data@colData
dim(clinical)

#"Primary solid Tumor" cases for survival
clin_df = clinical[clinical$definition == "Primary solid Tumor",
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "tumor_stage")]

# create a new boolean variable that has TRUE for dead patients
clin_df$deceased = clin_df$vital_status == "Dead"

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who are still alive
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)

#KM Plot
Surv(clin_df$overall_survival, clin_df$deceased)
Surv(clin_df$overall_survival, clin_df$deceased) ~ clin_df$gender

# fit a survival model
fit = survfit(Surv(overall_survival, deceased) ~ gender, data=clin_df)
print(fit)

#Final Kaplan Meier plot
ggsurvplot(fit, data=clin_df)
