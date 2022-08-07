# This code follows Seurat Pipeline in analyzing Single cell data 
#https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#If you have not done so before, Install Seurat 

# Enter commands in R (or R studio, if installed)
#install.packages('Seurat')

#load the library
{
library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)
library(stringi)
library(data.table)
library(ggplot2)
}
#set the working directory
setwd("D:/RProjects/CibersortxOvarian/Data/GSE118828_RAW")


#Get the names of the raw data files
RawDataFileNames <- list.files("D:/RProjects/CibersortxOvarian/Data/GSE118828_RAW")
TrimmeFilenAme <- stri_sub(RawDataFileNames, from=1, length=10)


#load the raw data for each sample

FirstDatSet <- as.data.frame(t(read.csv(RawDataFileNames[1], header=T,sep=',', row.names = 1)))
OV.combined <- CreateSeuratObject(counts = FirstDatSet, min.cells = 0, min.features = 0)



for (i in 2:length(RawDataFileNames)){
  NextDataSet <- as.data.frame(t(read.csv(RawDataFileNames[i], header=T,sep=',', row.names = 1)))
  NextDataSet <- CreateSeuratObject(counts = NextDataSet, min.cells = 0, min.features = 0)
  lISTA = paste(TrimmeFilenAme[i],"Data",sep="")
  
  OV.combined <- merge(OV.combined, y = NextDataSet, add.cell.ids = c("", lISTA), project = "OV_SigMatrix")
 
}

head(OV.combined)


#Create the Seurat Object

OV.combined

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

#Identify the level of mitochondrion genes in each of the cells. 
#The lower the percentage of mitochondrion genes, the better the health of the cell
OV.combined[["percent.mt"]] <- PercentageFeatureSet(OV.combined, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(OV.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#We can also visualize as feature scatter plots.
#This will help us decide on the QC metrics

plot1 <- FeatureScatter(OV.combined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(OV.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#Now remove outrageous/unwanted cells from the analysis
OV.combined <- subset(OV.combined, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

#Now perform the normalization

OV.combined <- NormalizeData(OV.combined, normalization.method = "LogNormalize", scale.factor = 10000)

#Now identify variable features

OV.combined <- FindVariableFeatures(OV.combined, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(OV.combined), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(OV.combined)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot2
#Now we perform scaling

all.genes <- rownames(OV.combined)
OV.combined <- ScaleData(OV.combined, features = all.genes)

#Perform dimensional reduction
OV.combined <- RunPCA(OV.combined, features = VariableFeatures(object = OV.combined))

# Examine and visualize PCA results a few different ways
print(OV.combined[["pca"]], dims = 1:5, nfeatures = 10)

VizDimLoadings(OV.combined, dims = 1:2, reduction = "pca")

DimPlot(OV.combined, reduction = "pca")

DimHeatmap(OV.combined, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(OV.combined, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the dimensionality of the datasets

OV.combined <- JackStraw(OV.combined, num.replicate = 100)
OV.combined <- ScoreJackStraw(OV.combined, dims = 1:20)
JackStrawPlot(OV.combined, dims = 1:15)
ElbowPlot(OV.combined)

#Now cluster the cells

OV.combined <- FindNeighbors(OV.combined, dims = 1:20)
OV.combinedres05 <- FindClusters(OV.combined, resolution = 0.5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
OV.combinedres05 <- RunUMAP(OV.combinedres05, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(OV.combinedres05, reduction = "umap", label = TRUE,label.size = 8)


saveRDS(OV.combinedres05, file = "OV.combinedres05clusted.rds")

# find markers for every cluster compared to all remaining cells, report only the positive ones
OV.combinedres05.markers <- FindAllMarkers(OV.combinedres05, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

Top10 <- OV.combinedres05.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#--Export the Gene Markers----------------------------------------------------------------------------------------
setwd("D:/RProjects/CibersortxOvarian/Output")
write.csv(x = OV.combinedres05.markers, file = "AllClusterOV.markers.csv", quote = FALSE)
write.csv(x = Top10, file = "AllCluster.markersTop10OV.csv", quote = FALSE,row.names=FALSE)

#Now we can identify clusters based on the markers provided in the article

Epithelial.MArkers <- c("KRT17",
                      "KRT6A",
                      "KLK10",
                      "KLK7",
                      "KLK8",
                      "KRT4",
                      "EPCAM",
                      "MMP7",
                      "CLDN3",
                      "WNT7A")

Fibroblasts.Markers <- c("PDPN"	,
                      "DCN"	,
                      "THY1"	,
                      "SNAI2"	,
                      "ALDH1A2"	,
                      "DCN"	,
                      "MMP2"	,
                      "CFI"	,
                      "DKK3"	,
                      "C1R")

Erythrocyte.Markers <- c("HBB"	,
                         "HBA2"	,
                         "HBA1"	,
                         "GATA1"	,
                         "KLF1")

Dendritic.Markers <- c("CD1C"	,
                         "CD1E"	,
                         "CCR7 "	,
                         "CD83",
                       "HLA-DRA"	,
                       "HLA-DPA1 "	,
                       "HLA-DQA1",
                       "CSF2RA"	,
                       "HLA-DQB2 "	,
                       "HLA-DQA2 "	,
                       "HLA-DPB1"	)

Macrophages.Markers <- c("CD14"	,
                         "AIF1"	,
                         "CSF1R"	,
                         "CD163",
                         "C1QB"	,
                         "C1QA"	,
                         "FCGR1A"	)

BCells.Markers <- c("CD19"	,
                         "CD79A"	,
                         "CD79B"	)

TCells.Markers <- c("CD3D"	,
                    "CD2"	,
                    "GZMB"	)
#We can't identify cluster 7 and 10

OV.combinedres05 <- RenameIdents(OV.combinedres05,`0` = "Tcells",`1` =  "Epithelial",
                                 `2` = "Fibroblasts",`3` =  "Macrophages",
                                 `4` = "Fibroblasts",`5` =  "Fibroblasts",
                                 `6` = "Epithelial",`8` =  "Fibroblasts",
                                 `9` = "Epithelial",`11` =  "Bcells",
                                 `12` = "Epithelial",`13` =  "Fibroblasts",
                                 `14` = "Epithelial")


VlnPlot(OV.combinedres05, features = Epithelial.MArkers)
VlnPlot(OV.combinedres05, features = Fibroblasts.Markers)
VlnPlot(OV.combinedres05, features = Erythrocyte.Markers)
VlnPlot(OV.combinedres05, features = Dendritic.Markers)
VlnPlot(OV.combinedres05, features = Macrophages.Markers)
VlnPlot(OV.combinedres05, features = BCells.Markers)
VlnPlot(OV.combinedres05, features = TCells.Markers)

DimPlot(OV.combinedres05, reduction = "umap", label = TRUE,label.size = 8)

#Now lets subset the data into individual cell types.

#-------------------------------------------------------------------------------#
Epithelial <-subset(OV.combinedres05, idents = "Epithelial")
Epithelial
Epithelial.small<-subset(Epithelial,subset = EPCAM > 3)
Epithelial.small

Epithelial.counts <- as.matrix(GetAssayData(Epithelial, slot = "counts"))
Epithelial.small.counts <- as.matrix(GetAssayData(Epithelial.small, slot = "counts"))


EpithelialName<-as.character(Epithelial@active.ident)
colnames(Epithelial.counts) <- EpithelialName

Epithelial.smallName<-as.character(Epithelial.small@active.ident)
colnames(Epithelial.small.counts) <- Epithelial.smallName
#-------------------------------------------------------------------------------#
Tcells <-subset(OV.combinedres05, idents = "Tcells")
Tcells
Tcells.small<-subset(Tcells,subset = CD3D > 2.5)
Tcells.small

Tcells.counts <- as.matrix(GetAssayData(Tcells, slot = "counts"))
Tcells.small.counts <- as.matrix(GetAssayData(Tcells.small, slot = "counts"))


TcellsName<-as.character(Tcells@active.ident)
colnames(Tcells.counts) <- TcellsName

Tcells.smallName<-as.character(Tcells.small@active.ident)
colnames(Tcells.small.counts) <- Tcells.smallName
#------------------------------------------------------------------------------#
Fibroblasts <-subset(OV.combinedres05, idents = "Fibroblasts")
Fibroblasts
Fibroblasts.small<-subset(Fibroblasts,subset = DCN > 4.7)
Fibroblasts.small

Fibroblasts.counts <- as.matrix(GetAssayData(Fibroblasts, slot = "counts"))
Fibroblasts.small.counts <- as.matrix(GetAssayData(Fibroblasts.small, slot = "counts"))


FibroblastsName<-as.character(Fibroblasts@active.ident)
colnames(Fibroblasts.counts) <- FibroblastsName

Fibroblasts.smallName<-as.character(Fibroblasts.small@active.ident)
colnames(Fibroblasts.small.counts) <- Fibroblasts.smallName
#-------------------------------------------------------------------------------#
Bcells <-subset(OV.combinedres05, idents = "Bcells")
Bcells
Bcells.small<-subset(Bcells,subset = CD79A > 0)
Bcells.small

Bcells.counts <- as.matrix(GetAssayData(Bcells, slot = "counts"))
Bcells.small.counts <- as.matrix(GetAssayData(Bcells.small, slot = "counts"))


BcellsName<-as.character(Bcells@active.ident)
colnames(Bcells.counts) <- BcellsName

Bcells.smallName<-as.character(Bcells.small@active.ident)
colnames(Bcells.small.counts) <- Bcells.smallName
#-------------------------------------------------------------------------------#
Macrophages <-subset(OV.combinedres05, idents = "Macrophages")
Macrophages
Macrophages.small<-subset(Macrophages,subset = AIF1 > 2.2)
Macrophages.small

Macrophages.counts <- as.matrix(GetAssayData(Macrophages, slot = "counts"))
Macrophages.small.counts <- as.matrix(GetAssayData(Macrophages.small, slot = "counts"))


MacrophagesName<-as.character(Macrophages@active.ident)
colnames(Macrophages.counts) <- MacrophagesName

Macrophages.smallName<-as.character(Macrophages.small@active.ident)
colnames(Macrophages.small.counts) <- Macrophages.smallName
#-------------------------------------------------------------------------------#
#Now merge the data to form the signature matrix

#-------------------------------------------------------------------------------#
#large Signature Matrix
l <- list(Epithelial.counts, Tcells.counts,Fibroblasts.counts,Bcells.counts,Macrophages.counts)
CombinedSigMatrix <-Epithelial.counts
for(i in 2:length(l)) {
  CombinedSigMatrix <- cbind(CombinedSigMatrix, l[[i]])
}
nrow(CombinedSigMatrix)
ncol(CombinedSigMatrix)
CombinedSigMatrix[1:5,100:107]
#-------------------------------------------------------------------------------#

#Small Signature Matrix
l.small <- list(Epithelial.small.counts, Tcells.small.counts, Fibroblasts.small.counts,
                Bcells.small.counts, Macrophages.small.counts)
CombinedSigMatrix.small <-Epithelial.small.counts
for(i in 2:length(l.small)) {
  CombinedSigMatrix.small <- cbind(CombinedSigMatrix.small, l.small[[i]])
}
nrow(CombinedSigMatrix.small)
ncol(CombinedSigMatrix.small)
CombinedSigMatrix.small[1:5,100:107]

#-------------------------------------------------------------------------------#
#Export the Signature Matrix

write.table(x = CombinedSigMatrix, file = "CombinedSigMatrix_03_28_2021.csv", quote = FALSE, col.names = NA, sep = " ")
write.table(x = CombinedSigMatrix.small, file = "CombinedSigMatrix.small_03_28_2021.csv", quote = FALSE, col.names = NA, sep = " ")

#Prepare the Mixture Matrix

MixtureMatrix <- read.csv("data_RNA_Seq_v2_expression_median.txt", header=T, sep = "")
MixtureMatrix <- na.omit(MixtureMatrix)
MixtureMatrix <-MixtureMatrix[,-2]
colnames(MixtureMatrix)[1] <- "Gene"
write.table(x = MixtureMatrix, file = "MixtureMatrix_03_28_2021.txt", quote = FALSE, sep = " ",row.names = FALSE)



