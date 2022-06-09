BiocManager::install('limma')
BiocManager::install('rmarkdown')
BiocManager::install('BioNet')
BiocManager::install('ClassComparison')
BiocManager::install('EBSeq')
library(BioNet)
library(limma)
library(rmarkdown)
library(ClassComparison)
library(tidyverse)
library(tidyr)
library(dplyr)
library(readxl)
library(EBSeq)
library(limma)
library(openxlsx)
library(stringi)
# library(preprocessCore)
persDEG <- function(cons, dis){
  require(dplyr)
  expDesign <- c(rep(1, ncol(cons)-1), rep(2, ncol(dis)-1))
  exprData <- log2(abs(as.matrix(cbind(cons[,-1],dis[,-1]))))
  exprData <- data.frame(GeneID = cons$GeneID,exprData)
  exprData <- aggregate(.~GeneID,exprData,max)
  rownames(exprData) <- exprData$GeneID
  exprData <- exprData[,-1]
  design <- model.matrix(~0 + factor(expDesign))#Adjust this based on the exprimental setup
  colnames(design) <- c('control', 'disease')#
  contrast <- makeContrasts(disease - control, levels = design)
  fit <- lmFit(exprData, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')[, c(1,4)]
  fc <- degDF[,1]
  for(j in 1:length(fc)){
    fc[j] <- ifelse(fc[j]>0, 2^fc[j],-1*2^abs(fc[j]))
  }
  Sizes <- rep(1,ncol(exprData))
  PPDE <- EBTest(Data=as.matrix(exprData), Conditions=as.factor(expDesign),
                  sizeFactors = Sizes, maxround = 10)
  PPDE <- GetPPMat(PPDE)
  
  return(data.frame(GeneID = rownames(exprData), PPDE = PPDE[,2], FC = fc))
}

fl <- 'D:/Projects/NIH/Mouse/Trancriptome/'
f1.files <- list.files(path = 'D:/Projects/NIH/Mouse/Trancriptome/', pattern = '.xlsx')
for(i in 1:length(f1.files)){
  i <- 1
  wb <- loadWorkbook(paste0(fl,f1.files[i]))
  stNames <- names(wb)
  cons <- readxl::read_xlsx(path = paste0(fl,f1.files[i]), sheet = stNames[1])
  dis <- readxl::read_xlsx(path = paste0(fl,f1.files[i]), sheet = stNames[2])
  cons.rem <- which(cons$GeneID==';' | cons$GeneID==''| cons$GeneID=='1'| cons$GeneID=='1;')
  cons <- cons[-cons.rem,]
  dis <- dis[-cons.rem,]
  # dgeRes <- persDEG(cons, dis)
  
  expDesign <- c(rep('control', ncol(cons)-1), rep('disease', ncol(dis)-1))
  exprData <- log2(abs(as.matrix(cbind(cons[,-1],dis[,-1])))+1)
  exprData <- data.frame(GeneID = cons$GeneID,exprData)
  exprData <- aggregate(.~GeneID,exprData,max)
  rownames(exprData) <- exprData$GeneID
  exprData <- exprData[,-1]
  design <- model.matrix(~0 + factor(expDesign))#Adjust this based on the exprimental setup
  colnames(design) <- c('control', 'disease')#
  contrast <- makeContrasts(disease - control, levels = design)
  fit <- lmFit(exprData, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  degDF <- topTable(fit2, n = 'Inf', sort.by = 'none')[, c(1,4)]
  fc <- degDF[,1]
  for(j in 1:length(fc)){
    fc[j] <- ifelse(fc[j]>0, 2^fc[j],-1*2^abs(fc[j]))
  }
  Sizes <- rep(1,ncol(exprData))
  exprData.1 <- exprData
  # exprData.1[] <- sapply(exprData.1, function(x){as.numeric(x)})
  exprData.1 <- apply(as.matrix.noquote(exprData),  # Using apply function
                      2,
                      as.numeric)
  row.names(exprData.1) <- rownames(exprData)
  colnames(exprData.1) <- NULL
  # exprData.1 <- matrix(data = as.numeric(exprData[,c(1:ncol(exprData))]),
  #                      ncol = ncol(exprData), nrow = nrow(exprData))
  PPDE <- EBTest(Data=exprData.1, Conditions= as.factor(expDesign),
                 sizeFactors = Sizes, maxround = 10)
  PPDE <- GetPPMat(PPDE)
  
}


# df1 <- readxl::read_xlsx(path = paste0(hGem)) %>% dplyr::select(ID, SUBSYSTEM)
mmEnrich_UP <- list()


for(i in 1:length(sheetNames)){
  # i <- 1
  df <- readxl::read_xlsx(path = paste0(fl), sheet = sheetNames[i])
}
data(GeneMat)
str(GeneMat)

Sizes <- rep(1,ncol(GeneMat))
EBOut <- EBTest(Data=GeneMat, Conditions=as.factor(rep(c("C1","C2"),each=5)),
               sizeFactors = Sizes, maxround = 10)
PP <- GetPPMat(EBOut)
FC <- EBSeq::PostFC(EBOut)


GSE129296_MIC_APPS1_3M <- read_excel("C:/Users/odong/Downloads/GSE129296_MIC_APPS1_3M.xlsx")
PPDE <- Bum(GSE129296_MIC_APPS1_3M$PPDE)
countSignificant(PPDE, 0.00, by='Emp')

eBayes()
PPDE.1 <- ClassComparison::selectSignificant(PPDE, alph = 0.5)

pvals <- data.frame(GSE129296_MIC_APPS1_3M$GeneID,GSE129296_MIC_APPS1_3M$PPDE)
colnames(pvals) <- c('GeneID','PVals')
pvals <- aggregate(.~GeneID,pvals,min)[-1,]
pvals.1 <- data.frame(pvals$PVals)
colnames(pvals.1) <- 'PVals'  
rownames(pvals.1) <- pvals$GeneID
PPDE <- BioNet::fitBumModel(pvals.1,starts = 20,plot = T)
summary.bum(PPDE)
PPDE$pvalues == pvals.1
PPDE.Lik <- likelihoodBum(PPDE)
selectSignificant()