setwd('D:/Projects/NIH/Mouse/Data/Comp/')
load('MethodologyComparison.RData')
library(dplyr)
library(tidyverse)
library(readxl)

BarnardStats_Mouse <- read_excel("BarnardStats_Mouse.xlsx", 
                                 sheet = "GSE129296_MIC_APPS1_12M")
moomin_Mouse <- read_excel("moominResults_1.xlsx", 
                              sheet = "GSE129296_MIC_APPS1_12M")

fishera_Mouse <- read_excel("GSE129296_MIC_APPS1_12M.xlsx", 
                                      sheet = "Fishers Test Results")
ElasticNet_MetabFeatures_iMAT <- read_excel("D:/Projects/NIH/ElasticNet_MetabFeatures_iMAT.xlsx")

# Select significant cases only based on p-value
BarnardStats_Mouse <- BarnardStats_Mouse %>% filter(.$BPValue<0.05)
fishera_Mouse <- fishera_Mouse %>% filter(.$`Fishers P-Value`<0.05)

pathways.MM <- table(moomin_Mouse$Subsystems) %>%
  as.data.frame()
pathways.FS <- table(fishera_Mouse$Subsystem) %>%
  as.data.frame()
pathways.BD <- table(BarnardStats_Mouse$Subsystem) %>%
  as.data.frame()
common_rxns_all <- intersect(moomin_Mouse$Reaction,
                         intersect(fishera_Mouse$Reaction, BarnardStats_Mouse$Reaction))
common_rxns_bf <- intersect(fishera_Mouse$Reaction, BarnardStats_Mouse$Reaction)
common_rxns_bm <- intersect(moomin_Mouse$Reaction,BarnardStats_Mouse$Reaction)
common_rxns_fm <- intersect(moomin_Mouse$Reaction,fishera_Mouse$Reaction)

common_rxns_bf <- fishera_Mouse %>% filter(.$Reaction==common_rxns_bf)
common_rxns_bm <- moomin_Mouse %>% filter(.$Reaction==common_rxns_bm)

# Plot Venn Diagrams on predicted reactions
BiocManager::install('ggVennDiagram')
library(ggVennDiagram)
library(ggplot2)
rxnList <- list(Barnards = BarnardStats_Mouse$Reaction,
                Fishers = fishera_Mouse$Reaction,
                MOOMIN = moomin_Mouse$Reaction)
ggVennDiagram(rxnList, label_alpha = 0,
              category.names = c("Barnards Exact Test","Fishers Exact Test","MOOMIN"))

# Load all annotation data
mouseCompartments <- read_excel('D:/Projects/NIH/mouseComps.xlsx')
HumanGEMetabolites <- read_excel('D:/Projects/NIH/Human-GEM.xlsx', 
                             sheet = 'METS')
HumanGEM <- read_excel('D:/Projects/NIH/Human-GEM.xlsx')

extractMetData <- function(humanGEM,
                           listOfReactions){
  paths.actMat <- table(humanGEM$SUBSYSTEM[which(humanGEM$ID %in% listOfReactions)]) %>% as.data.frame()
  mets.actMat <- humanGEM[which(humanGEM$ID %in% listOfReactions),] %>%
    select(ID,EQUATION,SUBSYSTEM)
  mets.actMat$MetabolitesSubstrate <- sapply(mets.actMat$EQUATION, 
                                             function(x){strsplit(x,'=',fixed = T)[[1]][1][1]})
  mets.actMat$MetabolitesProduct <- sapply(mets.actMat$EQUATION, 
                                           function(x){strsplit(x,'=',fixed = T)[[1]][1][1]})
  mets.actMat$MetabolitesSubstrate <- gsub('<','',mets.actMat$MetabolitesSubstrate)
  mets.actMat$MetabolitesSubstrate <- gsub('>','',mets.actMat$MetabolitesSubstrate)
  mets.actMat$MetabolitesSubstrate <- gsub('\\+',',',mets.actMat$MetabolitesSubstrate)
  
  mets.actMat$MetabolitesProduct <- gsub('<','',mets.actMat$MetabolitesProduct)
  mets.actMat$MetabolitesProduct <- gsub('>','',mets.actMat$MetabolitesProduct)
  mets.actMat$MetabolitesProduct <- gsub('\\+',',',mets.actMat$MetabolitesProduct)
  
  affPaths <- unique(mets.actMat$SUBSYSTEM)
  res.Data <- data.frame(matrix(data = NA, ncol = 5,nrow = length(affPaths)))
  colnames(res.Data) <- c('Reaction','Subsystem','Metabolites_Substrate','Metabolites_Product','Frequency')
  for(i in 1:length(affPaths)){
    df <- mets.actMat %>% filter(.$SUBSYSTEM==affPaths[i])
    res.Data[i,] <- c(toString(df$ID),
                      affPaths[i],
                      toString(df$MetabolitesSubstrate),
                      toString(df$MetabolitesProduct),
                      nrow(df))
  }
  res.Data <- res.Data[order(as.numeric(res.Data$Frequency), decreasing = T),]
  return(res.Data)
}
# For the GSE129296_MIC_APPS1_12M dataset
MetabData_Fishers <- extractMetData(humanGEM = HumanGEM,
                                    listOfReactions = fishera_Mouse$Reaction)
MetabData_Barnards <- extractMetData(humanGEM = HumanGEM,
                                    listOfReactions = BarnardStats_Mouse$Reaction)
MetabData_MOOMIN <- extractMetData(humanGEM = HumanGEM,
                                    listOfReactions = moomin_Mouse$Reaction)
# Remember, ElasticNet is general -- All samples included!
MetabData_ElasticNet <- extractMetData(humanGEM = HumanGEM,
                                   listOfReactions = ElasticNet_MetabFeatures_iMAT$rxnName)

# The most affected METABOLITE
affectedMetabo <- function(
    res.Data,
    Subsystem,
    Method
    ){
  mets <- res.Data %>% filter(res.Data$Subsystem == paste0(Subsystem))
  
  # Substrates
  metsSub <- strsplit(mets$Metabolites_Substrate,',')[[1]]
  for(i in 1:length(metsSub)){
    metsSub[i] <- stringr::str_remove(metsSub[i],'(?<=\\[).*(?=\\])')
    metsSub[i] <- gsub('\\[','',metsSub[i])
    metsSub[i] <- gsub('\\]','',metsSub[i])
  }
  metsSub <- table(metsSub) %>%
    as.data.frame()
  colnames(metsSub) <- c('Compartment','Frequency')
  metsSub <- metsSub[as.character(metsSub$Compartment)!=' ',]
  metsSub <- metsSub[order(metsSub$Frequency,decreasing = T),]
  
  png(filename = paste0('Plots/',Subsystem,Method,'Substrate.png'),
      width = 6.5, height = 5.5, res = 600, units = 'in')
  print(piePlot(metsSub,paste0('Metabolites in ',Subsystem)))
  dev.off()
  
  # Products
  metsPro <- strsplit(mets$Metabolites_Product,',')[[1]]
  for(i in 1:length(metsPro)){
    metsPro[i] <- stringr::str_remove(metsPro[i],'(?<=\\[).*(?=\\])')
    metsPro[i] <- gsub('\\[','',metsPro[i])
    metsPro[i] <- gsub('\\]','',metsPro[i])
  }
  metsPro <- table(metsPro) %>%
    as.data.frame()
  colnames(metsPro) <- c('Compartment','Frequency')
  metsPro <- metsPro[as.character(metsPro$Compartment)!=' ',]
  metsPro <- metsPro[order(metsPro$Frequency,decreasing = T),]
  
  png(filename = paste0('Plots/',Subsystem,Method,'Product.png'),
      width = 6.5, height = 5.5, res = 600, units = 'in')
  print(piePlot(metsPro,paste0('Metabolites in ',Subsystem)))
  dev.off()
  
  return(list(Sustrates = metsSub,
              Products = metsPro))
}
# Metabolites in a subsystem
# MOOMIN
MOOMIN_metabolites <- affectedMetabo( res.Data = MetabData_MOOMIN,
                                      Subsystem = 'Transport reactions',
                                      Method = 'MOOMIN')
# Fisher's
Fishers_metabolites <- affectedMetabo( res.Data = MetabData_Fishers,
                                       Subsystem = 'Transport reactions',
                                       Method = 'Fishers')
# Barnard's
Barnards_metabolites <- affectedMetabo( res.Data = MetabData_Barnards,
                                        Subsystem = 'Transport reactions',
                                        Method = 'Barnards')
# ElasticNet
ElsticNet_metabolites <- affectedMetabo(res.Data = MetabData_ElasticNet,
                                        Subsystem = 'Transport reactions',
                                        Method = 'ElasticNet')

# The most affected cellular COMPARTMENT
affectedCompart <- function( res.Data,
                             Subsystem,
                             Method,
                             mouseComps
){
  comps <- res.Data %>% filter(res.Data$Subsystem == paste0(Subsystem))
  
  # Substrates
  compsSubs <- strsplit(comps$Metabolites_Substrate,',')[[1]]
  for(j in 1:length(compsSubs)){
    compsSubs[j] <- stringr::str_extract(compsSubs[j],'(?<=\\[).*(?=\\])')
  }
  compsSubs <- compsSubs[!is.na(compsSubs)]
  compsSubs <- table(mouseComps$Names[match(compsSubs,mouseComps$Abbreviations)]) %>%
    as.data.frame()
  colnames(compsSubs) <- c('Compartment','Frequency')
  png(filename = paste0('Plots/Compartments_',Subsystem,Method,'_Substrate.png'),
      width = 6.5, height = 5.5, res = 600, units = 'in')
  print(piePlot(compsSubs,paste0('Metabolites in ',Subsystem)))
  dev.off()
  
  # Product
  compsProd <- strsplit(comps$Metabolites_Product,',')[[1]]
  for(j in 1:length(compsProd)){
    compsProd[j] <- stringr::str_extract(compsProd[j],'(?<=\\[).*(?=\\])')
  }
  compsProd <- compsProd[!is.na(compsProd)]
  compsProd <- table(mouseComps$Names[match(compsProd,mouseComps$Abbreviations)]) %>%
    as.data.frame()
  colnames(compsProd) <- c('Compartment','Frequency')
  png(filename = paste0('Plots/Compartments_',Subsystem,Method,'_Product.png'),
      width = 6.5, height = 5.5, res = 600, units = 'in')
  print(piePlot(compsProd,paste0('Metabolites in ',Subsystem)))
  dev.off()
  return(list(SubstrateCompartments = compsSubs,
              ProductCompartments = compsProd))
}

# Compartments in a subsystem
# MOOMIN
MOOMIN_compartment <- affectedCompart(res.Data = MetabData_MOOMIN,
                                      Subsystem = 'Transport reactions',
                                      Method = 'MOOMIN',
                                      mouseComps = mouseCompartments)
# Fisher's
Fishers_compartment <- affectedCompart( res.Data = MetabData_Fishers,
                                       Subsystem = 'Transport reactions',
                                       Method = 'Fishers',
                                       mouseComps = mouseCompartments)
# Barnard's
Barnards_compartment <- affectedCompart( res.Data = MetabData_Barnards,
                                        Subsystem = 'Transport reactions',
                                        Method = 'Barnards',
                                        mouseComps = mouseCompartments)
# ElasticNet
ElsticNet_compartment <- affectedCompart(res.Data = MetabData_ElasticNet,
                                        Subsystem = 'Transport reactions',
                                        Method = 'ElasticNet',
                                        mouseComps = mouseCompartments)
# Checking for Differential Expression of Leading Genes: Violin plots - Control vs AD

library(openxlsx)
geneList <- moomin_Mouse$LeadingGene
geneList <- geneList[!is.na(geneList)]

boxPlotter <- function(df, XLab,geneListJ){
  p <- ggplot(df, aes(x=Condition, y=Expression, fill=Condition)) +
    geom_violin() + 
    geom_boxplot(width = .2) +
    # geom_boxplot(position=position_dodge(1)) +
    ggtitle(paste0(geneListJ)) +
    xlab(paste0(XLab))
    theme_classic()
  return(p)
}
fl <- 'D:/Projects/NIH/Mouse/Trancriptome/'
f1.files <- list.files(path = 'D:/Projects/NIH/Mouse/Trancriptome/', pattern = '.xlsx')
for(i in 1:length(f1.files)){
  # i <- 1
  wb <- loadWorkbook(paste0(fl,f1.files[i]))
  stNames <- names(wb)
  cons <- readxl::read_xlsx(path = paste0(fl,f1.files[i]), sheet = stNames[1])
  dis <- readxl::read_xlsx(path = paste0(fl,f1.files[i]), sheet = stNames[2])
  for(j in 1:length(geneList)){
    # j <- 14
    if(geneList[j] %in% cons$GeneID){
      expre.con <- cons %>% filter(.$GeneID == geneList[j])
      expre.dis <- dis %>% filter(.$GeneID == geneList[j])
      if(t.test(expre.con[,-1],expre.dis[,-1], alternative = 'two.sided')$p.value < 0.05){
        ttestRes <- t.test(expre.con[,-1],expre.dis[,-1], alternative = 'two.sided')$p.value
        geneFC <- mean(as.numeric(expre.dis[,-1]))/mean(as.numeric(expre.con[,-1])) 
        geneFC <- ifelse(geneFC>1,geneFC,-1/geneFC)
        xLab <- paste0('p-value: ',round(ttestRes,digits = 4),' FC ',round(geneFC,digits = 2))
        df.con <- expre.con[,-1]
        colnames(df.con) <- rep('control',ncol(df.con))
        df.dis <- expre.dis[,-1]
        colnames(df.dis) <- rep('AD',ncol(df.dis))
        df <- cbind(df.con,df.dis) %>%
          pivot_longer(cols = starts_with('AD')|starts_with('con'),
                       names_to = 'Condition',
                       values_to = 'Expression')
        fname <- gsub('.xlsx','',f1.files[i])
        geneListJ <- paste0(geneList[j],'_',fname)
        png(filename = paste0('D:/Projects/NIH/Mouse/Data/Comp/Plots/',geneList[j],'_Gene_',fname,'.png'),
            width = 4.5, height = 4.5, res = 600, units = 'in')
        print(boxPlotter(df,xLab,geneListJ))
        dev.off()
      }
    }
  }
}
save.image('MethodologyComparison.RData')
