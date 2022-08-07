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


#---Run Magic on the raw datasets----------------------------------------------#

rawdata_oV_Num_MAGICPathways <- magic(rawdata_oV_Num, 
                               genes=c(Epithelial.MArkers,Fibroblasts.Markers,
                                       Macrophages.Markers,TCells.Markers,
                                       BCells.Markers,Compare.Genes,NAD.genes,
                                       CysMet1.genes, CysMet2.genes,
                                       CysMet3.genes,ETC.genes,ETC.UQCR11.genes,
                                       Onecarbon1.genes,Onecarbon2.genes))
rawdata_oV_Num_MAGICPathways <- magic(rawdata_oV_Num, 
                                      genes=c(Epithelial.MArkers,Fibroblasts.Markers,
                                              Macrophages.Markers,TCells.Markers,
                                              BCells.Markers,Compare.Genes,NAD.genes,
                                              CysMet1.genes, CysMet2.genes,
                                              CysMet3.genes,ETC.genes,ETC.UQCR11.genes,
                                              Onecarbon1.genes,Onecarbon2.genes),
                                              t=20, init = rawdata_oV_Num_MAGICPathways)

#---Run Magic on the raw datasets----------------------------------------------#
rawdata_oV_Num_MAGICPathwaysref2 <- magic(rawdata_oV_Num, 
                                      genes=c(Epithelial.MArkers,Fibroblasts.Markers,
                                              Macrophages.Markers,TCells.Markers,
                                              BCells.Markers,Compare.Genes,NAD.genes,
                                              CysMet1.genes, CysMet2.genes,
                                              CysMet3.genes,ETC.genes,ETC.UQCR11.genes,
                                              Onecarbon1.genes,Onecarbon2.genes,
                                              SGT1.genes,SGT2.genes))
rawdata_oV_Num_MAGICPathwaysref2 <- magic(rawdata_oV_Num, 
                                      genes=c(Epithelial.MArkers,Fibroblasts.Markers,
                                              Macrophages.Markers,TCells.Markers,
                                              BCells.Markers,Compare.Genes,NAD.genes,
                                              CysMet1.genes, CysMet2.genes,
                                              CysMet3.genes,ETC.genes,ETC.UQCR11.genes,
                                              Onecarbon1.genes,Onecarbon2.genes,SGT1.genes,SGT2.genes),
                                      t=20, init = rawdata_oV_Num_MAGICPathwaysref2)

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=UQCR11)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

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
#Generate marker plots
{
Plot1 <- ggplot(rawdata_oV_Num_PHATE) +
  geom_point(aes(x=PHATE1, y=PHATE2,color=rawdata_oV_Num_MAGICPathways$result$SNAI2)) + 
  scale_color_gradient2(low = "white", mid = "grey", high = "red" )+
  labs(color="SNAI2")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
  coord_fixed()
ggsave(Plot1, filename = "SNAI2_Fibroblasts.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}

clustersref2 <- clusters
clustersref2[clustersref2 == 4] <- 0
{
ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=clustersref2)) +
    scale_color_gradientn( colours = c("black","blue", "green", "red" ))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed(ratio = 2/5)
ggsave(filename = "MAGICPlot_ClusternamedRef2.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}
  #ggplot(df, aes(x=x, y=y)) + geom_point() + coord_fixed()
# color=rawdata_oV_Num_MAGIC$result$MTHFD2)

#----Pathway Genes lists-------------------------------------------------------#
{
NAD.genes <- c("NAMPT",
               "CD38",
               "BST1",
               "NNMT",
               "NT5E",
               "ENPP1",
               "AOX1")

CysMet1.genes <- c("LDHC",
                   "MAT1A",
                   "BHMT",
                   "DNMT1",
                   "DNMT3B",
                   "DNMT3A",
                   "TRDMT1")

CysMet2.genes <- c("AMD1",
                   "LDHB",
                   "GOT2",
                   "MPST",
                   "AHCY",
                   "ADI1",
                   "SRM")

CysMet3.genes <- c("SMS",
                   "APIP",
                   "LDHA",
                   "GOT1",
                   "IL4I1",
                   "SDS",
                   "ENOPH1",
                   "MAT2B")

ETC.genes <- c("ATP5F1",
               "NDUFS4",
               "COX4I1",
               "NDUFS7",
               "NDUFA6",
               "ATP5H",
               "ATP5G1",
               "NDUFB5",
               "NDUFS6",
               "ATP5G3",
               "COX5A",
               "NDUFS3",
               "COX5B",
               "NDUFA2",
               "UQCRQ",
               "COX7B",
               "NDUFA1",
               "UQCR11",
               "COX8A",
               "NDUFS8",
               "ATP5L",
               "ATP5J2",
               "NDUFB2",
               "COX6A1",
               "ATP5O",
               "ATP5J",
               "NDUFB3",
               "NDUFB4",
               "NDUFC1",
               "NDUFA3",
               "COX6B1",
               "NDUFA11",
               "NDUFB7",
               "NDUFB6",
               "COX7A2",
               "NDUFAB1",
               "NDUFB10",
               "ATP5I",
               "NDUFC2",
               "NDUFA8",
               "NDUFB9",
               "COX6C")

ETC.UQCR11.genes <- c("COX5B",
                      "NDUFA2",
                      "UQCRQ",
                      "COX7B",
                      "NDUFA1",
                      "UQCR11",
                      "COX8A",
                      "NDUFS8")

Onecarbon1.genes <- c("ATIC",
                      "MTHFD1L",
                      "SHMT2",
                      "MTHFD2",
                      "GART",
                      "MTHFD1",
                      "DHFR",
                      "TYMS",
                      "MTHFD2L",
                      "MTHFS",
                      "MTFMT")

Onecarbon2.genes <- c("ALDH1L1",
                      "FTCD",
                      "AMT",
                      "MTR",
                      "MTHFR")
SGT1.genes <- c("SHMT2",
                "PHGDH",
                "PSPH",
                "CHDH",
                "GATM",
                "MAOB",
                "MAOA",
                "BHMT",
                "DMGDH",
                "ALAS2")

SGT2.genes <- c("PSAT1",
                "SHMT1",
                "DLD",
                "GLDC",
                "CBS",
                "PIPOX")

}
#----Pathway Genes lists-------------------------------------------------------#

All.average.genes <- rawdata_oV_Num_MAGICPathways$result

#--NAD------
NAD.average.genes <- All.average.genes[,c("NAMPT","CD38","BST1","NNMT","NT5E","AOX1")]
NAD.average.genes$NADmeans <- apply(NAD.average.genes, 1, mean)

NAD.Pathway.Score <- NAD.average.genes$NADmeans
{
ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=NAD.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed(ratio = 2/5)
ggsave(filename = "NAD_PathwayScore.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}

#--One Carbon Metabolism 
OCM.average.genes <- All.average.genes[,Onecarbon1.genes]
OCM.average.genes$OCMmeans <- apply(OCM.average.genes, 1, mean)

OCM.Pathway.Score <- OCM.average.genes$OCMmeans

{
ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=OCM.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed(ratio = 2/5)
ggsave(filename = "OCM_PathwayScore.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}

#--ETC-UQCR11 Metabolism 
ETC_U.average.genes <- All.average.genes[,ETC.UQCR11.genes]
ETC_U.average.genes$ETC_Umeans <- apply(ETC_U.average.genes, 1, mean)

ETC_U.Pathway.Score <- ETC_U.average.genes$ETC_Umeans

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=ETC_U.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

#--ETC Metabolism 
ETC.average.genes <- All.average.genes[,ETC.genes]
ETC.average.genes$ETCmeans <- apply(ETC.average.genes, 1, mean)

ETC.Pathway.Score <- ETC.average.genes$ETCmeans

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=ETC.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

#--CysMet1 Metabolism 
CysMet1.average.genes <- All.average.genes[,c( "DNMT1",
                                              "DNMT3A",
                                              "TRDMT1")]
CysMet1.average.genes$CysMet1means <- apply(CysMet1.average.genes, 1, mean)

CysMet1.Pathway.Score <- CysMet1.average.genes$CysMet1means

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=CysMet1.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

#--CysMet2 Metabolism 
CysMet2.average.genes <- All.average.genes[,CysMet2.genes]
CysMet2.average.genes$CysMet2means <- apply(CysMet2.average.genes, 1, mean)

CysMet2.Pathway.Score <- CysMet2.average.genes$CysMet2means

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=CysMet2.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

#--CysMet3 Metabolism 
CysMet3.average.genes <- All.average.genes[,CysMet3.genes]
CysMet3.average.genes$CysMet3means <- apply(CysMet3.average.genes, 1, mean)

CysMet3.Pathway.Score <- CysMet3.average.genes$CysMet3means

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=CysMet3.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

#--OCM2 Metabolism 
OCM2.average.genes <- All.average.genes[,c("AMT",
                                           "MTR",
                                           "MTHFR")
                                        ]
OCM2.average.genes$OCM2means <- apply(OCM2.average.genes, 1, mean)

OCM2.Pathway.Score <- OCM2.average.genes$OCM2means

ggplot(rawdata_oV_Num_MAGICPathways) +
  geom_point(aes(MTHFD2, UQCR11, color=OCM2.Pathway.Score)) +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 

#--KEGG_SGT2 Metabolism 
All.averageRef2.genes <- rawdata_oV_Num_MAGICPathwaysref2$result
SGT2.average.genes <- All.averageRef2.genes[,c("PSAT1",
                                               "SHMT1",
                                               "DLD",
                                               "CBS")]
SGT2.average.genes$SGT2means <- apply(SGT2.average.genes, 1, mean)

SGT2.Pathway.Score <- SGT2.average.genes$SGT2means

{
  ggplot(rawdata_oV_Num_MAGICPathways) +
    geom_point(aes(MTHFD2, UQCR11, color=SGT2.Pathway.Score)) +
    scale_color_viridis(option="B")+
    cowplot::theme_cowplot() + 
    theme(axis.line  = element_blank()) +
    coord_fixed(ratio = 2/5)
  ggsave(filename = "SGT2_PathwayScore.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}


#--KEGG_SGT1 Metabolism 
All.averageRef2.genes <- rawdata_oV_Num_MAGICPathwaysref2$result
SGT1.average.genes <- All.averageRef2.genes[,c("SHMT2",
                                              "PHGDH",
                                              "PSPH",
                                              "GATM",
                                              "MAOB",
                                              "MAOA",
                                              "ALAS2")]
SGT1.average.genes$SGT1means <- apply(SGT1.average.genes, 1, mean)

SGT1.Pathway.Score <- SGT1.average.genes$SGT1means

{
  ggplot(rawdata_oV_Num_MAGICPathways) +
    geom_point(aes(MTHFD2, UQCR11, color=SGT1.Pathway.Score)) +
    scale_color_viridis(option="B")+
    cowplot::theme_cowplot() + 
    theme(axis.line  = element_blank()) +
    coord_fixed(ratio = 2/5)
  ggsave(filename = "SGT1_PathwayScore.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}



#---Spearman's correlation-----------------------------------------------------#
{
  AllDatatoPlot <- rawdata_oV_Num_MAGICPathways$result$UQCR11
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGICPathways$result$MTHFD2 )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGICPathways$result$EPCAM )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGICPathways$result$CD14 )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGICPathways$result$PDPN )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGICPathways$result$CD79A )
  AllDatatoPlot <- cbind(AllDatatoPlot,rawdata_oV_Num_MAGICPathways$result$GZMB )
  colnames(AllDatatoPlot) <- c("UQCR11", "MTHFD2","EPCAM","CD14","PDPN","CD79A","GZMB")
  
  AllDatatoPlot <- as.data.frame(AllDatatoPlot)
  AllDatatoPlot <- cbind(AllDatatoPlot,clusters) 
  }
AllDatatoPlot[AllDatatoPlot == 4] <- 0

#clusters1 <-cluster_phate(rawdata_oV_Num_PHATE, k = 5, seed = NULL)
{
ggscatter(AllDatatoPlot, x = "MTHFD2", 
           y = "UQCR11" , color = "clusters",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MTHFD2", ylab = "UQCR11") +
  scale_color_gradientn( colours = c("black","blue", "green", "red" ))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) +
coord_fixed(ratio = 2/5)
ggsave(filename = "PearsonCorrAllcells_Trendline.pdf", device = "pdf", path = "D:/RProjects/OV_Ascites_ScRNASeq_Analysis/Output")
}


EpithelialOnly <- AllDatatoPlot[AllDatatoPlot$clusters == 3,]
EpithelialOnly <- EpithelialOnly[EpithelialOnly$EPCAM > 0.07,]


ggscatter(EpithelialOnly, x = "UQCR11", 
          y = "MTHFD2" , color = "UQCR11",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "UQCR11", ylab = "MTHFD2") +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 



FibroblastsOnly <- AllDatatoPlot[AllDatatoPlot$clusters == 1,]
ggscatter(FibroblastsOnly, x = "UQCR11", 
          y = "MTHFD2" , color = "MTHFD2",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "UQCR11", ylab = "MTHFD2") +
  scale_color_viridis(option="B")+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 


ggscatter(AllDatatoPlot, x = "MTHFD2", 
          y = "UQCR11" , color = "clusters",
          add = "reg.line", add.params = list(color = "magenta"), conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MTHFD2", ylab = "UQCR11") +
  scale_color_gradientn( colours = c("black","blue", "green", "red" ))+
  cowplot::theme_cowplot() + 
  theme(axis.line  = element_blank()) 



