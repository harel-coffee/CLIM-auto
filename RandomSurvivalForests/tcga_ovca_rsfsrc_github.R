# Libraries

library(cgdsr)
library(RColorBrewer)
library(stringr)
library(randomForestSRC)

# Import TCGA data

mycgds <- CGDS("http://www.cbioportal.org/")

study <- getCancerStudies(mycgds)
ov_study <- "ov_tcga_pan_can_atlas_2018"
case_list <- getCaseLists(mycgds,ov_study)
complete_case <- case_list[6,1]

profiles <- getGeneticProfiles(mycgds,ov_study)[,c(1:2)]
clinical <- getClinicalData(mycgds,complete_case)
rnaseq2_profile <- profiles[9,1]
mut_profile <- profiles[7,1]
cna_profile <- profiles[6,1]

# User-defined gene list
{
  sgoc_genes <- read.csv(file = "sgoc_genes.txt", sep = ",", header = F, stringsAsFactors = F)
  del_genes <- read.csv(file = "del_genes.txt", sep = ",", header = F, stringsAsFactors = F)
}

# Import TCGA data into data.frame
{
  import_genes <- unique(sgoc_genes$V1)
  data <- getProfileData(mycgds, import_genes, rnaseq2_profile, complete_case)
  cna_data <- getProfileData(mycgds, del_genes$V1, cna_profile, complete_case)
}

# Retain samples common in both mut and rnaseq datasets
{
data <- data[colSums(!is.na(data)) > 0]
common_cases <- intersect(rownames(data), rownames(clinical))
data <- data[common_cases,]
clinical <- clinical[common_cases,]
cna_data <- cna_data[common_cases,]
}

status <- str_replace(clinical$OS_STATUS, ":DECEASED", "")
status <- str_replace(status, ":LIVING", "")

# Create object with relevant genes and clinical data that will be factors for the rfsrc function
{
  rm(ov_tcga)
  ov_tcga <- data
  ov_tcga$surv <- as.numeric(status)
  ov_tcga$time <- clinical$OS_MONTHS
  ov_tcga$age <- clinical$AGE
}

# Change keep_ov to select cohorts "all", "19pdel", "19pnondel"

keep_ov <- which(cna_data$STK11 < -0.4) # this condition is for 19pdel
# Run RFS on TCGA cohort
{
ov_tcga.obj <- rfsrc(Surv(time,surv) ~ ., ov_tcga[keep_ov,])
ov_tcga.vimp <- vimp(ov_tcga.obj)
print(vimp(ov_tcga.obj)$importance)

# Plot important variables
imp<-which(ov_tcga.vimp$importance>0.001)
imp_names <- ov_tcga.vimp$xvar.names[imp]
print(imp_names)
plot.variable(ov_tcga.obj, imp_names, partial = TRUE, smooth.lines = TRUE)
}

