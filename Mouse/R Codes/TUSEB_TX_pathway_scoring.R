

rm_data.list <- readRDS(file = 'G:\\ParkinsonsDatasets\\Combined Datasets\\Orhan\\mouseratfcresultsusingtwometabolicmodelgenesrata\\mouse_diff.RDS')
rr_data.list <- readRDS(file = 'G:\\ParkinsonsDatasets\\Combined Datasets\\Orhan\\mouseratfcresultsusingtwometabolicmodelgenesrata\\rat_diff.RDS')

rm_data.list.1 <- list(rm_data.list$GSE52584_ST_LRRK2,rm_data.list$GSE60080_ST_MPTP_1wk,
                       rm_data.list$GSE60080_ST_MPTP_8wk,rm_data.list$GSE89562_NST_MP,
                       rm_data.list$GSE4758_WB_aSYN, rm_data.list$GSE4788_SN_MPTP_MME,
                       rm_data.list$GSE4788_SN_MPTP_MML, rm_data.list$GSE7707_FC_MPTP,
                       rm_data.list$GSE7707_ST_MPTP, rm_data.list$GSE7707_MB_MPTP,
                       rm_data.list$GSE8030_ST_METH, rm_data.list$GSE8030_ST_MPTP,
                       rm_data.list$GSE13033_CO_HtrA2_KO, rm_data.list$GSE17542_2day_SN,
                       rm_data.list$GSE17542_10day_SN, rm_data.list$GSE17542_VTA_MPTP_2d,
                       rm_data.list$GSE17542_10day_VTA, 
                       rm_data.list$GSE60413_CE_Pink1_ko_6wk, rm_data.list$GSE60413_CE_Pink1_ko_12wk,
                       rm_data.list$GSE60413_CE_Pink1_ko_24wk, rm_data.list$GSE60413_CE_Pink1_ko_18mnth,
                       rm_data.list$GSE60413_MB_Pink1_ko_6wk, rm_data.list$GSE60413_MB_Pink1_ko_24wk,
                       rm_data.list$GSE60413_MB_Pink1_ko_18mnth, rm_data.list$GSE60413_ST_Pink1_ko_6wk,
                       rm_data.list$GSE60413_ST_Pink1_ko_12wk, rm_data.list$GSE60413_ST_Pink1_ko_24wk,
                       rm_data.list$GSE60413_ST_Pink1_ko_18mnth, rm_data.list$GSE60414_CE_Pink1_ko_6wk,
                       rm_data.list$GSE60414_CE_Pink1_ko_6mnth, rm_data.list$GSE60414_MB_Pink1_ko_6wk,
                       rm_data.list$GSE60414_MB_Pink1_ko_6mnth, rm_data.list$GSE60414_ST_Pink1_ko_6wk,
                       rm_data.list$GSE60414_ST_Pink1_ko_6mnth, rm_data.list$GSE66730_VMB_Tfr1_3wk,
                       rm_data.list$GSE66730_VMB_Tfr1_10wk, rm_data.list$GSE75000_VMA_Ercc1)

rr_data.list.1 <- list(rr_data.list$GSE93695_ST_LDOPA, rr_data.list$GSE74382_DST_lesion_saline,
                       rr_data.list$GSE58710_SN_6OHDA_wk1,
                       rr_data.list$GSE58710_SN_6OHDA_wk2, rr_data.list$GSE58710_SN_6OHDA_wk4,
                       rr_data.list$GSE58710_SN_6OHDA_wk6, rr_data.list$GSE58710_SN_6OHDA_wk16,
                       rr_data.list$GSE71968_HB_DJ1_KO, rr_data.list$GSE24233_ST_6OHDA)
library(plyr)
# Mouse
rm_data.list.df <- join_all(rm_data.list.1, by = 'Entrez', type = 'left', match = 'all')
rm_data.list.df <- data.frame(row.names = rm_data.list.df$Entrez, rm_data.list.df[,-1])
rm_data.list.df <- rm_data.list.df %>% dplyr::select(ends_with('logFC', ignore.case = T))
new_labs <- colnames(rm_data.list.df)
new_labs <- gsub('logFC', 'M', new_labs, ignore.case = F)
new_labs <- gsub('d', 'D', new_labs, ignore.case = F)
new_labs <- gsub('wk', 'W', new_labs, ignore.case = F)
new_labs <- gsub('mnth', 'M', new_labs, ignore.case = F)
new_labs <- gsub('Pink1', 'PNK', new_labs, ignore.case = F)
colnames(rm_data.list.df) <- new_labs
# rm_data.list.df_PV <- rm_data.list.df %>% dplyr::select(contains('Pval', ignore.case = T))

# Rat
rr_data.list.df <- join_all(rr_data.list.1, by = 'Entrez', type = 'left', match = 'all')
rr_data.list.df <- data.frame(row.names = rr_data.list.df$Entrez, rr_data.list.df[,-1])
rr_data.list.df <- rr_data.list.df %>% dplyr::select(ends_with('logFC', ignore.case = T))
new_labs <- colnames(rr_data.list.df)
new_labs <- gsub('logFC', 'R', new_labs, ignore.case = F)
new_labs <- gsub('wk', 'W', new_labs, ignore.case = F)
new_labs <- gsub('lesion_saline', 'SAL', new_labs, ignore.case = F)
colnames(rr_data.list.df) <- new_labs


# saveRDS(rm_data.list.df, file = 'G:\\ParkinsonsDatasets\\Rat & Mouse\\Rat & Mouse\\rm_data.list.df_FC.RDS')
# saveRDS(rr_data.list.df, file = 'G:\\ParkinsonsDatasets\\Rat & Mouse\\Rat & Mouse\\rr_data.list.df_FC.RDS')

###

# Mouse and Rat Human Gene Orthologs
ortholog_m <- read_xlsx(path = 'G:\\ParkinsonsDatasets\\Combined Datasets\\Orhan\\fcupdated\\mouse_orthologs.xlsx',
                        sheet = 1)
ortholog_m.1 <- merge(ortholog_m, modelIDs_ENTREZ, by.x = 'Human_Ensembl', by.y = 'ensembl_gene_id', all = T)
colnames(ortholog_m.1)[3] <- 'hEntrezID'


ortholog_r <- read_xlsx(path = 'G:\\ParkinsonsDatasets\\Combined Datasets\\Orhan\\fcupdated\\rat_orthologs.xlsx',
                        sheet = 1) 
ortholog_r.1 <- merge(ortholog_r, modelIDs_ENTREZ, by.x = 'Human_Ensembl', by.y = 'ensembl_gene_id', all = T)
colnames(ortholog_r.1)[3] <- 'hEntrezID'
###
# rat <- readRDS(file = 'Mouse_rat_FC.RDS')

rm_data_df <- data.frame(GeneID = ortholog_m.1$hEntrezID[match(row.names(rm_data.list.df), ortholog_m.1$Mouse_Entrez)], 
                         rm_data.list.df)

rr_data_df <- data.frame(GeneID = ortholog_r.1$hEntrezID[match(row.names(rr_data.list.df), ortholog_r.1$Rat_Entrez)], 
                         rr_data.list.df)
##

######
#Combining Human, Cell Line, Rat and Mouse Data 
library(plyr)
new_hFC_1 <- data.frame(GeneID = row.names(human_cell.ls$human), human_cell.ls$human)
new_cFC_1 <- data.frame(GeneID = row.names(human_cell.ls$cell), human_cell.ls$cell)

human_cell <- join(new_cFC_1, new_hFC_1, by = 'GeneID', type = 'left', match = 'all')
mouse_rat <- join(rm_data_df, rr_data_df, by = 'GeneID', type = 'left', match = 'all')
mouse_rat <- mouse_rat[!duplicated(mouse_rat$GeneID),]

saveRDS(new_hFC_1, file = 'G:\\ParkinsonsDatasets\\human.RDS')
saveRDS(new_cFC_1, file = 'G:\\ParkinsonsDatasets\\cell.RDS')
saveRDS(mouse_rat, file = 'G:\\ParkinsonsDatasets\\mouse_rat.RDS')


combined_hcrm <- join(human_cell,  mouse_rat, by = 'GeneID', type = 'left')
combined_hcrm <- data.frame(row.names = combined_hcrm$GeneID, combined_hcrm[, -1])

# combined_hcrm <- data.frame(row.names = combined_hcrm$GeneID, combined_hcrm[,-1])
samLabs <- colnames(combined_hcrm)
samLabs <- gsub('^FC_', '', samLabs, ignore.case = T)
# samLabs <- gsub('Pink1', 'PK', samLabs, ignore.case = T)
samLabs <- gsub('Rot', 'RT', samLabs, ignore.case = T)
samLabs <- gsub('lesion_saline', 'L-S', samLabs, ignore.case = T)
samLabs <- gsub('Tfr1.null', 'T1N', samLabs, ignore.case = T)
samLabs <- gsub('6OHDA', '6HDA', samLabs, ignore.case = T)
# samLabs <- gsub('mnth', 'MO', samLabs, ignore.case = T)
samLabs <- gsub('HtrA2', 'HA2', samLabs, ignore.case = T)
samLabs <- gsub('PM', 'MP', samLabs, ignore.case = T)
samLabs <- gsub('MPD', 'PD', samLabs, ignore.case = T)
# samLabs <- gsub('ILB', 'PD', samLabs, ignore.case = T)

colnames(combined_hcrm) <- samLabs
# combined_hcrm <- combined_hcrm[!duplicated(combined_hcrm$GeneID),]
# combined_hcrm1 <- data.frame(row.names = combined_hcrm$GeneID, combined_hcrm[,-1])
### Reference Gene Symbols
modelIDs_ref <- join(modelIDs_ENTREZ, modelIDs_HGNC, by = 'ensembl_gene_id')  

saveRDS(combined_hcrm, file = 'G:\\ParkinsonsDatasets\\Combined Datasets\\combined_old_hcrm_df.RDS')

###c_combined_hcrm.PV
set.seed(999)
cdf <- c_combined_hcrm.PV
na_count <- data.frame(rowSums(is.na(cdf))/ncol(cdf))
names(na_count) <- 'props'
to_keep <- data.frame(ifelse(na_count$props>.05, 0, 1))
names(to_keep) <- 'keep'
cdf <- data.frame(cdf[to_keep$keep ==1,]) %>% dplyr::select(-contains('GSE43490'))

## Genes expressed in very few samples
# na_countDF <- data.frame(Genes = modelIDs_ref$hgnc_symbol[match(row.names(na_count),modelIDs_ref$entrezgene_id)],
#                          na_count$props)
# colnames(na_countDF)[2] <- 'Percentage.of.NA'

# write.xlsx(na_countDF, file = 'G:\\ParkinsonsDatasets\\Combined Datasets\\na_countDF.xlsx')
# mydata <- combined_hcrm_df_clean
# na_prop <- sum(rowSums(is.na(mydata)))/(nrow(mydata)*ncol(mydata))*100
combined_hcrm_imp_cdf <- missForest(xmis = cdf, maxiter = 200, ntree = 15,
                                variablewise = FALSE,
                                decreasing = FALSE, verbose = FALSE,
                                mtry = floor(sqrt(ncol(cdf))), replace = TRUE,
                                classwt = NULL, cutoff = NULL, strata = NULL,
                                sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                                xtrue = NA, parallelize = "no")
combined_hcrm_imp_cdf <- data.frame(combined_hcrm_imp_cdf$ximp)
###
df_pathw_genes <- readxl::read_xlsx(path = 'G:\\ParkinsonsDatasets\\Combined Datasets\\human-GEM-MetPathGenes.xlsx')

dataDF <- combined_hcrm_imp_1
dataDF <- data.frame(GeneSym = modelIDs_ref$hgnc_symbol[match(row.names(dataDF), modelIDs_ref$entrezgene_id)],
                     dataDF)
dataDF_PV <- combined_hcrm_impPV
dataDF_PV <- data.frame(GeneSym = modelIDs_ref$hgnc_symbol[match(row.names(dataDF_PV), modelIDs_ref$entrezgene_id)],
                        dataDF_PV)
# row.names(dataDF) <- as.character(modelIDs_ref$hgnc_symbol[match(row.names(combined_hcrm), modelIDs_ref$entrezgene_id)])
path_Scores <- list()
path_i_score <- data.frame(matrix(data = NA, nrow = ncol(df_pathw_genes), ncol = ncol(dataDF_PV)))
for(n in 1:ncol(df_pathw_genes)){
  row.names(path_i_score)[n] <- colnames(df_pathw_genes[,n])
  path_n <- df_pathw_genes[,n]
  path_n <- path_n[!is.na(path_n)]
  dataDF_metFC <- dataDF[dataDF$GeneSym %in% path_n,]
  dataDF_metFC <- data.frame(GeneSym = dataDF_metFC$GeneSym, dataDF_metFC[,-1])
  dataDF_metPV <- dataDF_PV[dataDF_PV$GeneSym %in% path_n,]
  dataDF_metPV <- data.frame(GeneSym = dataDF_metPV$GeneSym, dataDF_metPV[,-1])
  dataDF_metRes <- data.frame(abs(as.matrix(dataDF_metFC[,-1])) * -log10(as.matrix(dataDF_metPV[,-1])))
  for(m in 1:ncol(dataDF_metRes)){
    path_i_score[n, m] <- sum(dataDF_metRes[, m])/sqrt(nrow(dataDF_metRes))
  }
}
colnames(path_i_score) <- colnames(dataDF_PV[2:ncol(dataDF_PV)])

saveRDS(path_i_score, file = 'G:\\ParkinsonsDatasets\\Combined Datasets\\Pathway Scores\\path_i_score.RDS')

#This function accepts three inputs as follows: df_FC - a data frame of imputed fold-changes
#df_PV - data frame of imputed p-values (both of these data frames must contain HUGO gene symbols
# as the first column) and df_Paths - data frame containing pathways with their 
#associated genes
pathScore <- function(df_FC, df_PV, df_Paths){
  path_i_score <- data.frame(matrix(data = NA, nrow = ncol(df_Paths), ncol = ncol(df_FC)-1))
  for(n in 1:ncol(df_Paths)){
    row.names(path_i_score)[n] <- colnames(df_Paths[,n])
    path_n <- df_Paths[,n]
    path_n <- path_n[!is.na(path_n)]
    dataDF_metFC <- df_FC[df_FC$GeneSym %in% path_n,]
    dataDF_metFC <- as.matrix(dataDF_metFC[,-1])
    # dataDF_metFC <- data.frame(GeneSym = dataDF_metFC$GeneSym, dataDF_metFC[,-1])
    dataDF_metPV <- df_PV[df_PV$GeneSym %in% path_n,]
    dataDF_metPV <- as.matrix(dataDF_metPV[,-1])
    # dataDF_metPV <- data.frame(GeneSym = dataDF_metPV$GeneSym, dataDF_metPV[,-1])
    dataDF_metRes <- data.frame(abs(dataDF_metFC) * -log10(dataDF_metPV))
    for(m in 1:ncol(dataDF_metRes)){
      path_i_score[n, m] <- sum(dataDF_metRes[, m])/sqrt(nrow(dataDF_metRes))
    }
  }
  colnames(path_i_score) <- colnames(df_FC[2:ncol(df_PV)])
  return(path_i_score)
}
# hsamlabs <- colnames(df_data1)
# c <- combined_hcrm_imp_cdf %>%dplyr::select(ends_with('PD')|ends_with('PC')) %>%
#   dplyr::filter(row.names(combined_hcrm_imp_cdf) %in% row.names(df_data1))
# cLab <- colnames(c)
# cLab <- gsub('cEMTAB812_PC','cEMTAB812_PC_PD',cLab)
# colnames(c) <- cLab
# c <- c[,str_order(colnames(c), decreasing = T)]
# c <- data.frame(GeneSym = modIDs$hgnc_symbol[match(row.names(c),modIDs$entrezgene_id)], c)
# d <- df_data1[,str_order(colnames(df_data1), decreasing = T)]
# d <- data.frame(GeneSym = modIDs$hgnc_symbol[match(row.names(d),modIDs$entrezgene_id)], d)
# humanPark_MetPath_cor <- pathScore(d, c, df_pathw_genes)
# pathScores <- data.frame(PathName = row.names(humanPark_MetPath_cor),
#                          PathScore = rowSums(humanPark_MetPath_cor))

hyperGeoScore <- function(df_FC, df_Paths){
  pathStatScore <- c(rep(NA,ncol(df_Paths)))
  for(n in 1:ncol(df_Paths)){
    # row.names(pathStatScore)[n] <- colnames(df_Paths[,n])
    path_n <- df_Paths[,n]
    path_n <- as.character(path_n[!is.na(path_n)])
    path_n1 <- length(path_n)
    total_genes <- 3628# ro should we use numeric(nrow(df_FC))?
    common_genes <- length(intersect(path_n, df_FC$GeneSym))
    observed_genes <- length(df_FC$GeneSym)
    pathStatScore[n] <- phyper(common_genes-1, path_n1, total_genes-path_n1, observed_genes,lower.tail= FALSE)
    # dataDF_metFC <- dataDF[dataDF$GeneSym %in% path_n,]
    # for(m in 1:ncol(df_Paths)){
    #   pathStatScore[m] <- phyper(common_genes-1, path_n, total_genes-path_n, observed_genes,lower.tail= FALSE)
    # }
  }
  # colnames(pathStatScore) <- colnames(dataDF_PV[2:ncol(df_FC)])
  FDR <-  p.adjust(pathStatScore, method='fdr', n=length(pathStatScore))
  pathStatScore_df <- data.frame(PathName = colnames(df_Paths),P.Val = pathStatScore, adjP.Val = FDR)
  return(pathStatScore_df)
}
humanPark_MetPath_pval <- hyperGeoScore(d,df_pathw_genes)
pathScoresDF <- join(pathScores, humanPark_MetPath_pval, by = 'PathName', type = 'left', match = 'all')
pathScoresDF.1 <- pathScoresDF[!is.na(pathScoresDF$PathScore),]

pathScoresDF.2 <- pathScoresDF.1[order(pathScoresDF.1$P.Val),] %>%
  head(10) 
pathScoresDF.2 <- pathScoresDF.2[,-4]
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots/MetPathPlots/top10.png', 
    width = 5, height = 3.5, units = 'in', res = 600)
pathScoresDF.2 %>% ggplot(aes(y=PathName,x=PathScore, size = -log10(P.Val), fill = P.Val)) +
  geom_point() +
  theme_classic()
dev.off()
###
#Pvalues: combined_hcrm
combined_hcrmPV <- combined_hcrm %>% dplyr::select(-contains('GSE43490'))
saveRDS(combined_hcrmPV, file = 'G:/ParkinsonsDatasets/Combined Datasets/combined_hcrmPV.RDS')
combined_hcrmPV <- readRDS(file = 'G:/ParkinsonsDatasets/Combined Datasets/combined_hcrm_pv.RDS')
set.seed(999)
df <- combined_hcrmPV
colnames(df) <- as.character(coln$...1)
na_count <- data.frame(rowSums(is.na(df))/ncol(df))
names(na_count) <- 'props'
to_keep <- data.frame(ifelse(na_count$props>.05, 0, 1))
names(to_keep) <- 'keep'
df_PV <- data.frame(df[to_keep$keep ==1,])

combined_hcrm_impPV <- missForest(xmis = df_PV, maxiter = 200, ntree = 15,
                                variablewise = FALSE,
                                decreasing = FALSE, verbose = FALSE,
                                mtry = floor(sqrt(ncol(df_PV))), replace = TRUE,
                                classwt = NULL, cutoff = NULL, strata = NULL,
                                sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                                xtrue = NA, parallelize = "no")
combined_hcrm_impPV <- data.frame(combined_hcrm_impPV$ximp)
colnames(combined_hcrm_impPV) <- as.character(coln$...1)
##
modIDs <- modelIDs_ref[!is.na(modelIDs_ref$entrezgene_id),]
df_PV1 <- combined_hcrm_impPV %>% dplyr::filter(as.character(row.names(combined_hcrm_impPV)) %in% as.character(p_10.15))
df_data1 <- combined_hcrm_imp_1 %>% dplyr::filter(as.character(row.names(combined_hcrm_imp_1)) %in% as.character(p_10.15))
### Substantia nigra
SN <- c('GSE49036_SN_SPD','GSE26927_SN_SPD','GSE8397_LSN_SPD',
        'GSE7621_SN_SPD','GSE8397_MSN_SPD','GSE20292_SN_SPD')
hSN <- df_data1 %>% dplyr::select(contains(SN))
hSN <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hSN),modIDs$entrezgene_id)],hSN)
hSN.PV <- df_PV1 %>% dplyr::select(contains(SN))
hSN.PV <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hSN.PV),modIDs$entrezgene_id)],hSN.PV)

PC <- c('GSE20168_PC_SPD','GSE68719_PC_SPD',
             'GSE8397_SFG_SPD','EMTAB812_PC_SPD')
hPC <- df_data1 %>% dplyr::select(contains(PC))
hPC <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hPC),modIDs$entrezgene_id)],hPC)
hPC.PV <- df_PV1 %>% dplyr::select(contains(PC))
hPC.PV <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hPC.PV),modIDs$entrezgene_id)],hPC.PV)

ME <- c('GSE19587_DNV_SPD','GSE19587_ION_SPD')
hME <- df_data1 %>% dplyr::select(contains(ME))
hME <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hME),modIDs$entrezgene_id)],hME)
hME.PV <- df_PV1 %>% dplyr::select(contains(ME))
hME.PV <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hME.PV),modIDs$entrezgene_id)],hME.PV)

CE <- c('GSE28894_CE_SPD','GSE20314_CE_SPD')
hCE <- df_data1 %>% dplyr::select(contains(CE))
hCE <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hCE),modIDs$entrezgene_id)],hCE)
hCE.PV <- df_PV1 %>% dplyr::select(contains(CE))
hCE.PV <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hCE.PV),modIDs$entrezgene_id)],hCE.PV)

SN.m <- c('GSE4788_SN_MPTP_MME_M','GSE4788_SN_MPTP_MML_M')
mSN <- df_data1 %>% dplyr::select(contains(SN.m))
mSN <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(mSN),modIDs$entrezgene_id)],mSN)
mSN.PV <- df_PV1 %>% dplyr::select(contains(SN.m))
mSN.PV <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(mSN.PV),modIDs$entrezgene_id)],mSN.PV)

SN.r <- c('GSE58710_SN_6HDA_W1_R','GSE58710_SN_6HDA_W2_R',
          'GSE58710_SN_6HDA_W6_R', 'GSE58710_SN_6HDA_W16_R')
rSN <- df_data1 %>% dplyr::select(contains(SN.r))
rSN <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(rSN),modIDs$entrezgene_id)],rSN)
rSN.PV <- df_PV1 %>% dplyr::select(contains(SN.r))
rSN.PV <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(rSN.PV),modIDs$entrezgene_id)],rSN.PV)
####
# hSN <- df_data1 %>% select(contains(c('GSE26927_SN_PD','GSE20292_SN_PD','GSE7621_SN_PD',
#                                       'GSE49036_SN_PD','GSE20164_SN_PD','GSE8397_LSN_PD',
#                                       'GSE8397_MSN_PD')))
# hSN <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hSN),modIDs$entrezgene_id)],hSN)
# hPC <- df_data1 %>% select(contains(c('GSE20168_PC_PD','GSE68719_PC_PD',
#                                       'GSE8397_SFG_PD','EMTAB812_PC_PD')))
# hPC <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hPC),modIDs$entrezgene_id)],hPC)
# ###
#hSN_mSN
humanPark_MetPath_hSN_mSN <- list(hSN = pathScore(hSN, hSN.PV, df_pathw_genes),
                                  rSN = pathScore(rSN, rSN.PV, df_pathw_genes),
                                  mSN = pathScore(mSN, mSN.PV, df_pathw_genes),
                                  hME = pathScore(hME, hME.PV, df_pathw_genes),
                                  hCE = pathScore(hCE, hCE.PV, df_pathw_genes),
                                  hPC = pathScore(hPC, hPC.PV, df_pathw_genes),
                                  pVal = hyperGeoScore(hSN,df_pathw_genes))

###
library(ggpubr)
theme_set(theme_pubr())
####
hSN_mSN <- data.frame(PathScore = apply(as.matrix(humanPark_MetPath_hSN_mSN[["hSN"]]),
                                        1, mean),humanPark_MetPath_hSN_mSN[["pVal"]])
hSN_mSN <- hSN_mSN[,c(2,1,3)]
hSN_mSN.h <- head(hSN_mSN[order(hSN_mSN$PathScore, decreasing = T),],10)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/hSN.png', 
    width = 5.5, height = 7.5, units = 'in', res = 600)
ggplot(data = hSN_mSN.h) + aes(PathScore, P.Val, size = -log10(P.Val), color = PathName, label = PathName) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf)) +
  labs(title = 'Human Substantia nigra')+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
dev.off()
###
hSN_mSN <- data.frame(PathScore = apply(as.matrix(humanPark_MetPath_hSN_mSN[["mSN"]]),
                                        1, mean),humanPark_MetPath_hSN_mSN[["pVal"]])
hSN_mSN <- hSN_mSN[,c(2,1,3)]
hSN_mSN.h <- head(hSN_mSN[order(hSN_mSN$PathScore, decreasing = T),],10)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/mSN.png', 
    width = 5.5, height = 7, units = 'in', res = 600)
ggplot(data = hSN_mSN.h) + aes(PathScore, P.Val, size = -log10(P.Val), color = PathName, label = PathName) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf))+
  labs(title = 'Mouse Substantia nigra')+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
dev.off()
###
hSN_mSN <- data.frame(PathScore = apply(as.matrix(humanPark_MetPath_hSN_mSN[["hPC"]]),
                                        1, mean),humanPark_MetPath_hSN_mSN[["pVal"]])
hSN_mSN <- hSN_mSN[,c(2,1,3)]
hSN_mSN.h <- head(hSN_mSN[order(hSN_mSN$PathScore, decreasing = T),],10)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/hPC.png', 
    width = 5.5, height = 7, units = 'in', res = 600)
ggplot(data = hSN_mSN.h) + aes(PathScore, P.Val, size = -log10(P.Val), color = PathName, label = PathName) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf))+
  labs(title = 'Human Prefrontal Cortex')+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
dev.off()
###
hSN_mSN <- data.frame(PathScore = apply(as.matrix(humanPark_MetPath_hSN_mSN[["hME"]]),
                                        1, mean),humanPark_MetPath_hSN_mSN[["pVal"]])
hSN_mSN <- hSN_mSN[,c(2,1,3)]
hSN_mSN.h <- head(hSN_mSN[order(hSN_mSN$PathScore, decreasing = T),],10)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/hME.png', 
    width = 5.5, height = 7, units = 'in', res = 600)
ggplot(data = hSN_mSN.h) + aes(PathScore, P.Val, size = -log10(P.Val), color = PathName, label = PathName) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf))+
  labs(title = 'Human Medulla')+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
dev.off()
###
hSN_mSN <- data.frame(PathScore = apply(as.matrix(humanPark_MetPath_hSN_mSN[["hCE"]]),
                                        1, mean),humanPark_MetPath_hSN_mSN[["pVal"]])
hSN_mSN <- hSN_mSN[,c(2,1,3)]
hSN_mSN.h <- head(hSN_mSN[order(hSN_mSN$PathScore, decreasing = T),],10)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/hCE.png', 
    width = 5.5, height = 7, units = 'in', res = 600)
ggplot(data = hSN_mSN.h) + aes(PathScore, P.Val, size = -log10(P.Val), color = PathName, label = PathName) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf))+
  labs(title = 'Human Cerebellum')+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
dev.off()
###
hSN_mSN <- data.frame(PathScore = apply(as.matrix(humanPark_MetPath_hSN_mSN[["rSN"]]),
                                        1, mean),humanPark_MetPath_hSN_mSN[["pVal"]])
hSN_mSN <- hSN_mSN[,c(2,1,3)]
hSN_mSN.h <- head(hSN_mSN[order(hSN_mSN$PathScore, decreasing = T),],10)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/rSN.png', 
    width = 5.5, height = 7, units = 'in', res = 600)
ggplot(data = hSN_mSN.h) + aes(PathScore, P.Val, size = -log10(P.Val), color = PathName, label = PathName) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf))+
  labs(title = 'Rat Substantia nigra')+
  theme(legend.position = 'none',plot.title = element_text(hjust = 0.5))
dev.off()
#hSN, hSN.PV
path_n <- df_pathw_genes[,1]
path_n <- path_n[!is.na(path_n)]
dataDF_PathnFC <- rSN[rSN$GeneSym %in% path_n,]
dataDF_PathnFC <- data.frame(GeneSym=modIDs$hgnc_symbol[match(dataDF_PathnFC$GeneSym,modIDs$ensembl_gene_id)], dataDF_PathnFC[,-1])

dataDF_PathnFC.2 <- dataDF_PathnFC %>%tidyr::pivot_longer(cols = starts_with('GSE', ignore.case = T),
                                                          values_to = 'FC', names_to = 'Samples')
dataDF_PathnFC.1 <- apply(dataDF_PathnFC[,-1], 1, mean)
dataDF_PathnPV <- hSN.PV[hSN.PV$GeneSym %in% path_n,]
dataDF_PathnPV.1 <- apply(dataDF_PathnPV[,-1], 1, mean)
dataDF_PathnPV.1.Gly <- data.frame(GeneSym=dataDF_PathnFC$GeneSym, FC = dataDF_PathnFC.1, PVal = dataDF_PathnPV.1)

png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/PathScore/gluco_glycon_pc_rNS.png', 
    width = 7.5, height = 7, units = 'in', res = 600)
ggplot(data = dataDF_PathnFC.2) + aes(Samples, FC, color = GeneSym, label = GeneSym) +
  geom_point() +
  geom_text_repel(xlim = c(NA, Inf),ylim = c(-Inf, Inf)) +
  theme(legend.position = "none") +
  labs(title = 'Glycolysis/Gluconeogenesis metabolism')+
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=6.5))
dev.off()
###

###
# Regression on Brain Regions

hSN.1 <- apply(abs(as.matrix(hSN[,-1])), 1, mean)
rSN.1 <- apply(abs(as.matrix(rSN[,-1])), 1, mean)
mSN.1 <- apply(abs(as.matrix(mSN[,-1])), 1, mean)
hME.1 <- apply(abs(as.matrix(hME[,-1])), 1, mean)
hCE.1 <- apply(abs(as.matrix(hCE[,-1])), 1, mean)
hPC.1 <- apply(abs(as.matrix(hPC[,-1])), 1, mean)
###
library(geodist)
# library(circlize)
#
df <- data.frame(y = hSN.1, x = hPC.1)
colnames(df) <- c('hSN', 'hPC')
df_fit <- lm(hPC ~ hSN, data = df)
# cor_stat <- Hmisc::rcorr(t(scale(as.matrix(t(df)), center = T, scale = T)), type = 'pearson')
anova(df_fit)
z <- summary(df_fit)
df <- df %>%
  rownames_to_column(var = 'GeneID')
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/Regression/pc_sn_human.png', 
    width = 7.0, height = 4.5, units = 'in', res = 550)
ggplot(df, aes(hSN, hPC)) +
  geom_point() +
  stat_smooth(method = lm)+
  labs(title = 'Human: S.nigra vs P.cortex') +
  annotate(geom = 'text', x = 2.5, y = 2.5, label = paste0('R^2 =', format(z$r.squared, digits = 2),'\np = ',format(z$cov.unscaled[2,2], digits = 2)),
           color = 'black')
dev.off()
###
df <- data.frame(y = hSN.1, x = hME.1)
colnames(df) <- c('hSN', 'hME')
df_fit <- lm(hME  ~ hSN, data = df)
# cor_stat <- Hmisc::rcorr(t(scale(as.matrix(t(df)), center = T, scale = T)), type = 'pearson')
anova(df_fit)
z <- summary(df_fit)
df <- df %>%
  rownames_to_column(var = 'GeneID')
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/Regression/sn_me_human.png', 
    width = 7.0, height = 4.5, units = 'in', res = 550)
ggplot(df, aes(hSN, hME)) +
  geom_point() +
  stat_smooth(method = lm)+
  labs(title = 'Human: S.nigra vs Medulla') +
  annotate(geom = 'text', x = 5, y = 6, label = paste0('R^2 =', format(z$r.squared, digits = 2),'\np = ',format(z$cov.unscaled[2,2], digits = 2)),
           color = 'black')
dev.off()
###
df <- data.frame(y = hSN.1, x = hCE.1)
colnames(df) <- c('hSN', 'hCE')
df_fit <- lm(hCE  ~ hSN, data = df)
# cor_stat <- Hmisc::rcorr(t(scale(as.matrix(t(df)), center = T, scale = T)), type = 'pearson')
anova(df_fit)
z <- summary(df_fit)
df <- df %>%
  rownames_to_column(var = 'GeneID')
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/Regression/sn_ce_human.png', 
    width = 7.0, height = 4.5, units = 'in', res = 550)
ggplot(df, aes(hSN, hCE)) +
  geom_point() +
  stat_smooth(method = lm)+
  labs(title = 'Human: S.nigra vs Cerebellum') +
  annotate(geom = 'text', x = 5, y = 4, label = paste0('R^2 =', format(z$r.squared, digits = 2),'\np = ',format(z$cov.unscaled[2,2], digits = 2)),
           color = 'black')
dev.off()
###
df <- data.frame(y = hSN.1, x = mSN.1)
colnames(df) <- c('hSN', 'mSN')
df_fit <- lm(hSN ~ mSN, data = df)
# cor_stat <- Hmisc::rcorr(t(scale(as.matrix(t(df)), center = T, scale = T)), type = 'pearson')
anova(df_fit)
z <- summary(df_fit)
df <- df %>%
  rownames_to_column(var = 'GeneID')
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/Regression/sn_mouse_human.png', 
    width = 7.0, height = 4.5, units = 'in', res = 550)
ggplot(df, aes(mSN, hSN)) +
  geom_point() +
  stat_smooth(method = lm)+
  labs(title = 'Substantia nigra: Mouse vs Human') +
  annotate(geom = 'text', x = 0.5, y = 5, label = paste0('R^2 =', format(z$r.squared, digits = 2),'\np = ',format(z$cov.unscaled[2,2], digits = 2)),
           color = 'black')
dev.off()
##
df <- data.frame(y = hSN.1, x = rSN.1)
colnames(df) <- c('hSN', 'rSN')
df_fit <- lm(hSN ~ rSN, data = df)
# cor_stat <- Hmisc::rcorr(t(scale(as.matrix(t(df)), center = T, scale = T)), type = 'pearson')
anova(df_fit)
z <- summary(df_fit)
df <- df %>%
  rownames_to_column(var = 'GeneID')
#df[df$GeneID %in% c(1644,489,55244,2027),] <- 0
#df <- df[!rowSums(df[,2:3])==0,]
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/Regression/sn_rat_human2.png', 
    width = 6.5, height = 4.5, units = 'in', res = 550)
ggplot(df, aes(rSN, hSN)) +
  geom_point() +
  stat_smooth(method = lm)+
  labs(title = 'Substantia nigra: Rat vs Human') +
  annotate(geom = 'text', x = 1.5, y = 5.5, label = paste0('R^2 =', format(z$r.squared, digits = 2),'\np = ',format(z$cov.unscaled[2,2], digits = 2)),
           color = 'black')
dev.off()
##
###
library("ggpubr")
hSN.1 <- apply(abs(as.matrix(hSN[,-1])), 1, mean)
rSN.1 <- apply(abs(as.matrix(rSN[,-1])), 1, mean)
mSN.1 <- apply(abs(as.matrix(mSN[,-1])), 1, mean)
hME.1 <- apply(abs(as.matrix(hME[,-1])), 1, mean)
hPC.1 <- apply(abs(as.matrix(hPC[,-1])), 1, mean)
#PC
df <- data.frame(y = hSN.1, x = hPC.1)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/SpCorrelation/sn_pc_human.png', 
    width = 6.5, height = 4.5, units = 'in', res = 550)
ggscatter(df, x = 'y', y = 'x', 
          add = 'reg.line', conf.int = TRUE, 
          cor.coef = TRUE, cor.method = 'spearman',
          xlab = 'Human SN', ylab = 'Human CX')
dev.off()
#ME
df <- data.frame(y = hSN.1, x = hME.1)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/SpCorrelation/sn_me_human.png', 
    width = 6.5, height = 4.5, units = 'in', res = 550)
ggscatter(df, x = 'y', y = 'x', 
          add = 'reg.line', conf.int = TRUE, 
          cor.coef = TRUE, cor.method = 'spearman',
          xlab = 'Human SN', ylab = 'Human ME')
dev.off()
#Rat
df <- data.frame(y = hSN.1, x = rSN.1)
df$GeneID <- row.names(df)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/SpCorrelation/sn_human_rat.png', 
    width = 6.5, height = 4.5, units = 'in', res = 550)
ggscatter(df, x = 'y', y = 'x', 
          add = 'reg.line', conf.int = TRUE, 
          cor.coef = TRUE, cor.method = 'spearman',
          xlab = 'Human SN', ylab = 'Rat SN')
dev.off()
###Cerebellum
df <- data.frame(y = hSN.1, x = hCE.1)
df$GeneID <- row.names(df)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/SpCorrelation/sn_ce.png', 
    width = 6.5, height = 4.5, units = 'in', res = 550)
ggscatter(df, x = 'y', y = 'x', 
          add = 'reg.line', conf.int = TRUE, 
          cor.coef = TRUE, cor.method = 'spearman',
          xlab = 'Human SN', ylab = 'Human CE')
dev.off()
#Mouse
df <- data.frame(y = hSN.1, x = mSN.1)
png(filename = 'G:/ParkinsonsDatasets/Combined Datasets/Plots2/SpCorrelation/sn_human_mouse.png', 
    width = 6.5, height = 4.5, units = 'in', res = 550)
ggscatter(df, x = 'y', y = 'x', 
          add = 'reg.line', conf.int = TRUE, 
          cor.coef = TRUE, cor.method = 'spearman',
          xlab = 'Human SN', ylab = 'Mouse SN') 
dev.off()
####
hSN <- df_data1 %>% select(contains(c('GSE26927_SN_PD','GSE20292_SN_PD','GSE7621_SN_PD',
                                      'GSE49036_SN_PD','GSE20164_SN_PD','GSE8397_LSN_PD',
                                      'GSE8397_MSN_PD')))
# hSN <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hSN),modIDs$entrezgene_id)],hSN)
hPC <- df_data1 %>% select(contains(c('GSE20168_PC_PD','GSE68719_PC_PD',
                                      'GSE8397_SFG_PD','EMTAB812_PC_PD')))
# hPC <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hPC),modIDs$entrezgene_id)],hPC)
hME <- df_data1 %>% select(contains(c('GSE19587_DNV_PD','GSE19587_ION_PD',
                                      'GSE28894_ME_PD')))
# hME <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hME),modIDs$entrezgene_id)],hME)
hCE <- df_data1 %>% select(contains(c('GSE28894_CE_PD','GSE20314_CE_PD')))
# hCE <- data.frame(GeneSym = modIDs$ensembl_gene_id[match(row.names(hCE),modIDs$entrezgene_id)],hCE)
###
human <- df_data1 %>% dplyr::select(ends_with('PD'))
models <- df_data1 %>% dplyr::select(ends_with('M')| ends_with('R')| ends_with('C'))
df <- data.frame(i = apply(as.matrix(hME),1,mean), j = apply(as.matrix(hCE),1,mean))
res <- c(rep(NA, nrow(df)))
for(i in 1:nrow(df)){
  gen_i <- df[i,]
  if(gen_i[1] <0 & gen_i[2] < 0){
    res[i] <- (-1)
  }
  else if(gen_i[1] >0 & gen_i[2] > 0){
    res[i] <- 1
  }
  else{
    res[i] <- 0
  }
}
##
#Modified Jaccard Index
#
df <- data.frame(i = apply(as.matrix(human),1,mean), j = apply(as.matrix(models),1,mean))
res <- matrix(data = NA, nrow = nrow(df), ncol = 2)
for(i in 1:nrow(df)){
  for(j in 1:ncol(df)){
    if(df[i,j] < (-1.2) ){
      res[i,j] <- (-1)
    }
    else if(df[i,j] >1.2 ){
      res[i,j] <- 1
    }
    else{
      res[i,j] <- 0
    }
  }
  # gen_i <- df[i,]
}
res <- data.frame(res)
#
for(j in 1:nrow(res)){
  gene_j <- res[j,]
  if((gene_j[1] > 0 & gene_j[2] > 0) | (gene_j[1] < 0 & gene_j[2] < 0)){
    num <- j+1
  }
  else if((gene_j[1] > 0 | gene_j[1] < 0) | (gene_j[2] > 0 | gene_j[2] < 0)){
    den <- j+1
  }
}
jac.idx <- num/den
###
colnames (hSN) <- colnames (hPC) <- c ("x", "y")
sn_pc <- geodist(as.matrix(hSN), as.matrix(hPC))

saveRDS(modIDs, file = 'G:/ParkinsonsDatasets/Combined Datasets/modIDs.RDS')
