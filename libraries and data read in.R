###loading the libraries and ther essential data

library(FNN)
library(data.table) 
library(googlesheets4)
library(vioplot)
library(bitops)
library(ggplot2)
library(dplyr)
library(fastICA)
library(uwot)
library(RColorBrewer)
library(gridExtra)  
library(tidyr)
library(umap)
library(spatstat)
library(readr)
library(dplyr)
library(ggplot2)
library(car)   # for leveneTest()
library(FSA)   # for Dunn if you ever need it
library(ggvenn)
library(patchwork)
library(destiny)
library(princurve)
library(slingshot)
library(SingleCellExperiment)
library(writexl)
library(plotly)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(ggrepel)   
library(car)




load("listdatacysift2_Phen_Ki_BT_TCF1_tDBSCAN_BTinteract.RData")

markers<- c("Hoechst", "Ki67", "AMCAR", "HMWCK", "CD19", "SMA", "CD20", "CD11b", "CD68", "CD163", "CD4", "CD3d", "CD8a", "TCF1", "FOXP3", "PD1", "CD57", "CD11c", "GranzymeB", "CD15", "HLADR", "CD103", "CD31", "pTBK1", "HLAA", "CD24", "CD44", "CD206")
markers<- markers[length(markers):1]
mpflag<- numeric(length(markers))
for (i in 0:(length(markers)-1)){
  mpflag[i+1]<- 2^i
}
mpflag<- mpflag[length(mpflag):1]

