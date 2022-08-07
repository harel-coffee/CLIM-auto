#install.packages("Rmagic")
#install.packages('R.utils')
#------------------------------------------------------------------------------#
#Installing PhateR

if (!suppressWarnings(require(devtools))) install.packages("devtools")
reticulate::py_install("phate", pip=TRUE)
devtools::install_github("KrishnaswamyLab/phateR")
#------------------------------------------------------------------------------#
{
  if (!require(viridis)) install.packages("viridis")
  if (!require(ggplot2)) install.packages("ggplot2")
  if (!require(readr)) install.packages("readr")
  if (!require(phateR)) install.packages("phateR")
}

{
  library(Rmagic)
  ## Loading required package: Matrix
  library(ggplot2)
  library(cowplot)
  ## Warning: package 'ggplot2' was built under R version 3.5.3
  library(readr)
  library(viridis)
  ## Loading required package: viridisLite
  library(phateR)
  library(data.table)
  library("ggpubr")
}

setwd("D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Data")
bmmsc <- read_tsv("GSE146026_Izar_HGSOC_ascites_10x_log.tsv.gz", quote ="\"" ,strings)
bmmsc <- fread("GSE146026_Izar_HGSOC_ascites_10x_log.tsv.gz")

{
rawdata_oV <- as.data.frame(t(bmmsc))
rawdata_oV[1:10,1:10]
rawdata_oV <- rawdata_oV[,8:ncol(rawdata_oV)]
rawdata_oV[1:10,1:10]
colnames(rawdata_oV) <- rawdata_oV[1,]
rawdata_oV <- rawdata_oV[-1,]
rawdata_oV[1:10,1:10]
rownameslist <-rownames(rawdata_oV)
colnameslist <-colnames(rawdata_oV)
}


rawdata_oV_Num <- matrix(nrow = nrow(rawdata_oV), ncol = ncol(rawdata_oV))


for(column in 1:ncol(rawdata_oV)){
  rawdata_oV_Num[, column] <- as.numeric(rawdata_oV[,column])
}
rownames(rawdata_oV_Num)  <-  rownameslist
colnames(rawdata_oV_Num)  <-  colnameslist
rawdata_oV_Num <- as.data.frame(rawdata_oV_Num)



# keep genes expressed in at least 10 cells
keep_cols <- colSums(rawdata_oV_Num > 0) > 10
rawdata_oV_Num <- rawdata_oV_Num[,keep_cols]


# look at the distribution of library sizes
ggplot() +
  geom_histogram(aes(x=rowSums(rawdata_oV_Num)), bins=50) +
  geom_vline(xintercept = 1000, color='red')

# keep cells with at least 1000 UMIs
keep_rows <- rowSums(rawdata_oV_Num) > 1000
rawdata_oV_Num <- rawdata_oV_Num[keep_rows,]

rawdata_oV_Num <- library.size.normalize(rawdata_oV_Num)
rawdata_oV_Num <- sqrt(rawdata_oV_Num)

#rawdata_oV_Num_MAGIC <- magic(rawdata_oV_Num, genes=c("KRT18", "VIM", "MTHFD2","UQCR11"))

rawdata_oV_Num_MAGIC3 <- magic(rawdata_oV_Num, 
                               genes=c(Epithelial.MArkers,Fibroblasts.Markers,
                                       Macrophages.Markers,TCells.Markers,
                                       BCells.Markers,Compare.Genes))
                              

ggplot(rawdata_oV_Num) +
  geom_point(aes(KRT18, VIM, color=MTHFD2)) +
  scale_color_viridis(option="B")

ggplot(rawdata_oV_Num) +
  geom_point(aes(EPCAM, PDPN, color=UQCR11)) +
  scale_color_viridis(option="B")

ggplot(rawdata_oV_Num_MAGIC2) +
  geom_point(aes(EPCAM, PDPN, color=MTHFD2)) +
  scale_color_viridis(option="B")

ggplot(rawdata_oV_Num_MAGIC2) +
  geom_point(aes(EPCAM, PDPN, color=UQCR11)) +
  scale_color_viridis(option="B")

rawdata_oV_Num_MAGIC3 <- magic(rawdata_oV_Num, 
                               genes=c(Epithelial.MArkers,Fibroblasts.Markers,
                                       Macrophages.Markers,TCells.Markers,
                                       BCells.Markers,Compare.Genes), 
                               t=20, init =rawdata_oV_Num_MAGIC3)

ggplot(rawdata_oV_Num_MAGIC4) +
  geom_point(aes(MTHFD2, UQCR11, color=UQCR11)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 


ggplot(rawdata_oV_Num_MAGIC2) +
  geom_point(aes(EPCAM, PDPN, color=UQCR11)) +
  scale_color_viridis(option="B")

ggplot(rawdata_oV_Num_MAGIC2) +
  geom_point(aes(EPCAM, PDPN, color=MTHFD2)) +
  scale_color_viridis(option="B")

ggplot(rawdata_oV_Num_MAGIC3) +
  geom_point(aes(MTHFD2, UQCR11, color=GZMB)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

ggplot(rawdata_oV_Num_MAGIC3) +
  geom_point(aes(MTHFD2, UQCR11, color=EPCAM)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

ggplot(rawdata_oV_Num_MAGIC3) +
  geom_point(aes(MTHFD2, UQCR11, color=MTHFD2)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

ggplot(rawdata_oV_Num_MAGIC3) +
  geom_point(aes(MTHFD2, UQCR11, color=blues9)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())

ggplot(rawdata_oV_Num_MAGIC3) +
  geom_point(aes(MTHFD2, UQCR11, color=PDPN)) +
  scale_color_viridis(option="B")+
cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 



ggplot(rawdata_oV_Num_MAGIC2) +
  geom_point(aes(EPCAM, PDPN, color=UQCR11)) +
  scale_color_viridis(option="B")

rawdata_oV_Num_PHATE <- phate(rawdata_oV_Num, knn=5, decay=100, t=20)

clusters <-cluster_phate(rawdata_oV_Num_PHATE, k = 5, seed = NULL)
{
ggplot(rawdata_oV_Num_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2,color=factor(clusters))) + 
  scale_color_viridis(discrete = TRUE, option="C") +
  scale_fill_viridis(discrete = TRUE) +
cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed()
ggsave(filename = "PHATEplot.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}
#labs(color="EPCAM")
#color=rawdata_oV_Num_MAGIC$result$MTHFD2)

ggplot(rawdata_oV_Num_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2,color=clusters)) + 
  scale_color_gradientn( colours = c("cyan","grey", "blue", "red","green" ))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 
#labs(color="EPCAM")
#color=rawdata_oV_Num_MAGIC$result$MTHFD2)

ggplot(rawdata_oV_Num_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2,color=rawdata_oV_Num_MAGIC3$result$MTHFD2)) + 
  scale_color_gradient2(low = "white", mid = "grey", high = "blue" )+
 labs(color="MTHFD2")+
cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 
# color=rawdata_oV_Num_MAGIC$result$MTHFD2)

#Now we can identify clusters based on the markers found in the literature
{
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

Compare.Genes <- c("MTHFD2","UQCR11")
}

genesCombine <- c(Epithelial.MArkers,Fibroblasts.Markers,
        Macrophages.Markers,TCells.Markers,BCells.Markers)

#---Spearman's correlation-----------------------------------------------------#
{
  AllDatatoPlot <- rawdata_oV_Num_MAGIC3$result$UQCR11
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGIC3$result$MTHFD2 )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGIC3$result$EPCAM )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGIC3$result$CD14 )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGIC3$result$PDPN )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGIC3$result$CD79A )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGIC3$result$GZMB )
  colnames(AllDatatoPlot) <- c("UQCR11", "MTHFD2","EPCAM","CD14","PDPN","CD79A","GZMB")
  
  AllDatatoPlot <- as.data.frame(AllDatatoPlot)
  AllDatatoPlot <- cbind(AllDatatoPlot,clusters) 
}
AllDatatoPlot[AllDatatoPlot == 4] <- 0

#clusters1 <-cluster_phate(rawdata_oV_Num_PHATE, k = 5, seed = NULL)

ggscatter(AllDatatoPlot, x = "UQCR11", 
          y = "MTHFD2" , color = "PDPN",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "UQCR11", ylab = "MTHFD2") +
  scale_color_gradientn( colours = c("black","blue", "green", "red" ))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 


EpithelialOnly <- AllDatatoPlot[AllDatatoPlot$clusters == 3,]
EpithelialOnly <- EpithelialOnly[EpithelialOnly$EPCAM > 0.07,]


ggscatter(EpithelialOnly, x = "UQCR11", 
          y = "MTHFD2" , color = "EPCAM",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "UQCR11", ylab = "MTHFD2") +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 
{
ggscatter(EpithelialOnly, x = "UQCR11", 
          y = "MTHFD2" , color = "EPCAM",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "UQCR11", ylab = "MTHFD2") +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed(ratio = 1.5)
ggsave(filename = "EpithelialPearson'scorrelationNOTrendline.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}


FibroblastsOnly <- AllDatatoPlot[AllDatatoPlot$clusters == 1,]
ggscatter(FibroblastsOnly, x = "UQCR11", 
          y = "MTHFD2" , color = "MTHFD2",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "UQCR11", ylab = "MTHFD2") +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

{
ggscatter(AllDatatoPlot, x = "UQCR11", 
          y = "MTHFD2" , color = "clusters",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MTHFD2", ylab = "UQCR11") +
  scale_color_gradientn( colours = c("black","blue", "green", "red" ))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed(ratio = 2.3)
ggsave(filename = "AllcellsPearson'scorrelationWithTrendline.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
  }
{
  ggscatter(AllDatatoPlot, x = "UQCR11", 
            y = "MTHFD2" , color = "clusters",
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "MTHFD2", ylab = "UQCR11") +
    scale_color_gradientn( colours = c("black","blue", "green", "red" ))+
    cowplot::theme_cowplot() + 
    theme(axis.line  = element_blank()) +
    coord_fixed(ratio = 2.3)
  ggsave(filename = "AllcellsPearson'scorrelationNOTrendline.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}

{
ggplot(rawdata_oV_Num_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2,color=rawdata_oV_Num_MAGIC3$result$MTHFD2)) + 
  scale_color_gradient2(low = "white", mid = "grey", high = "red" )+
  labs(color="MTHFD2")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank())

}