setwd('D:/Projects/NIH/')
load('moominGEM.RData')
library(readxl)
library(openxlsx)
library(dplyr)
library(VarfromPDB)
library(RISmed)
library(stringi)

fl <- 'G:/Alzheimer Disease/GEM Simulations/Mouse  - Alzheimer/MOOMIN/Results/resMOOMIN.xlsx'
hGem <- 'G:/Alzheimer Disease/GEM Simulations/Mouse  - Alzheimer/Human-GEM.xlsx'
wb <- loadWorkbook(paste0(fl))
sheetNames <- names(wb)

# df1 <- readxl::read_xlsx(path = paste0(hGem)) %>% dplyr::select(ID, SUBSYSTEM)
mmEnrich_UP <- list()
mmEnrich_DW <- list()

mmEnrich_UP.1 <- list()
mmEnrich_DW.1 <- list()

for(i in 1:length(sheetNames)){
  df <- readxl::read_xlsx(path = paste0(fl), sheet = sheetNames[i])
  df <- df[df$`Output Color Code`!=6,]
  df.down <- df %>% filter(df$`Output Color Code`<0)
  df.down.1 <- df.down
  df.down <- table(df.down$subSystem) %>% as.data.frame()
  colnames(df.down) <- c('Pathway',paste0(sheetNames[i]))
  
  df.up <- df %>% filter(df$`Output Color Code`>0)
  df.up.1 <- df.up
  df.up <- table(df.up$subSystem) %>% as.data.frame()
  colnames(df.up) <- c('Pathway',paste0(sheetNames[i]))
  
  mmEnrich_DW[[paste0(sheetNames[i],'_DOWN')]] <- df.down[order(df.down[,2], decreasing = T),] %>% head(20)
  mmEnrich_UP[[paste0(sheetNames[i],'_UP')]] <- df.up[order(df.up[,2], decreasing = T),] %>% head(20)
  
  mmEnrich_DW.1[[paste0(sheetNames[i],'_DOWN')]] <- df.up.1
  mmEnrich_UP.1[[paste0(sheetNames[i],'_UP')]] <- df.down.1 
}

library(plyr)
library(tidyr)
library(tidyverse)
mmEnrich_DWN2 <- join_all(mmEnrich_DW,by = 'Pathway', type = 'full', match = 'all') %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Pathway')
mmEnrich_DWN2[is.na(mmEnrich_DWN2)] <- 0
  
mmEnrich_UP2 <- join_all(mmEnrich_UP,by = 'Pathway', type = 'full', match = 'all') %>%
  remove_rownames() %>%
  column_to_rownames(var = 'Pathway')
mmEnrich_UP2[is.na(mmEnrich_UP2)] <- 0

BR <- c('MIC','MIC','MIC','MIC','CX','HP','BR','MIC','MIC','MIC','BR','HP','HP','CX',
        'CX','CX','CX','HP','HP','HP','FCX','DG','MIC','MIC')
BR <- sapply(colnames(mmEnrich_UP2), function(x){strsplit(x,'_',fixed = T)[[1]][2]})#[[1]]
GB <- c('APP/PS1','APP/PS1','TAU','TAU','APP/PS1','APP/PS1','APP/PS1','FAD-APP','E3E4-FAD','E3E4','E3E4','B-CELL','B-CELL','B-CELL',
        'FAD','FAD','FAD','FAD','FAD','FAD','APP/PS1','FNDC5','R47H','R47H-TAU')
BR2 <- as.integer(factor(BR))
GB2 <- as.integer(factor(GB))
GB_Path1 <- cor(t(mmEnrich_UP2),GB2)
GB_Path2 <- cor(t(mmEnrich_DWN2),GB2)

BR_Path1 <- cor(t(mmEnrich_UP2),BR2)
BR_Path2 <- cor(t(mmEnrich_DWN2),BR2)
res.DF1 <- data.frame(Pathway = row.names(),) 
# BR2 <- as.numeric(levels(BR))[BR]
BiocManager::install('ComplexHeatmap')
library(circlize)
library(ComplexHeatmap)
library(pheatmap)

# By Brain Region
annotation_col <- data.frame(
  'Brain Region' =BR)
rownames(annotation_col) <- colnames(mmEnrich_UP2)


pheatmap(mmEnrich_UP2,
         show_rownames=T,
         show_colnames=FALSE,
         # kmeans_k = 6,
         # cluster_cols = T,
         # cluster_rows = F,
         annotation_col=annotation_col,
         scale = "row",
         clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

# By Genetic Background
annotation_col <- data.frame(
  'Genetic Background' =GB,
  'Brain Region' =BR)
rownames(annotation_col) <- colnames(mmEnrich_UP2)


pheatmap(mmEnrich_UP2,
         show_rownames=T,
         show_colnames=FALSE,
         # kmeans_k = 6,
         # cluster_cols = T,
         # cluster_rows = F,
         annotation_col=annotation_col,
         scale = 'none',
         clustering_method="ward.D2",
         clustering_distance_cols="euclidean")

# col_fun <- colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
# Heatmap(as.matrix(mmEnrich_UP2), name = 'Frequency', col = col_fun,column_km = 10)
# 
moominRes <- function(lsDF){
  resList <- list()
  for(i in 1:length(lsDF)){
    df <- lsDF[[i]]
    df$Metabolites <- sapply(df$`Metabolic Reaction`, 
                                      function(x){strsplit(x,'=',fixed = T)[[1]][1][1]})
    df$Metabolites <- gsub('<','',df$Metabolites)
    df$Metabolites <- gsub('\\+',',',df$Metabolites)
    affPaths <- unique(df$subSystem)
    res.Data <- data.frame(matrix(data = NA, ncol = 5,nrow = length(affPaths)))
    colnames(res.Data) <- c('Subsystem','Reactions','Metabolites','Genes','Enrichment')
    for(j in 1:length(affPaths)){
      # i <- 1
      df.1 <- df %>% filter(.$subSystem==affPaths[j])
      res.Data[j,] <- c(affPaths[j],
                        toString(df.1$`Reaction Name`),
                        toString(df.1$Metabolites),
                        gsub('NA,','',toString(unique(df.1$`Leading Gene`)),ignore.case = F),
                        nrow(df.1))
    }
    resList[[i]] <- res.Data[order(as.numeric(res.Data$Enrichment),decreasing = T),]
  }
  names(resList) <- names(lsDF)
  return(resList)
}
mmEnrich_UP.2 <- moominRes(mmEnrich_UP.1)
mmEnrich_DW.2 <- moominRes(mmEnrich_DW.1)
save.image('moominGEM.RData')
