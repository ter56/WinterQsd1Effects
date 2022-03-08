library(sommer);library(arm);library(lme4);library(Hmisc);library(readxl);
library(tibble);library(patchwork);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr)
library(BGLR)
library(tidyverse)
library(multidplyr)
library(dplyr);library(plyr)
load('GenotypeData/WinterGD.RData')
load('GenotypeData/WinterGM.RData')
WinterGM = WinterGM %>% arrange(Chromosome, Position)
WinterGD = WinterGD %>%
  dplyr::mutate(taxa = ifelse(taxa == "SY Tepee",'SY_Tepee',taxa),
                taxa = gsub(pattern = '-',replacement = '_',x = taxa),
                taxa = gsub(pattern = "[[:space:]]", replacement = '_',x = taxa),
                taxa1 = taxa) %>% remove_rownames()%>% 
  column_to_rownames('taxa1') %>%
  dplyr::mutate(Qsd1 = ifelse(taxa == 'SY_Tepee',0,Qsd1)) %>%arrange(taxa)
#SY_Tepee Family is segregating for Qsd1. Tepee must be 0 not 2. 

variousModelPAdetermination = function(yNA, myGDrr = DH_GD_GP_asMarkersrr){
  set.seed(1)
  Qsd1 = myGDrr$Qsd1
  inter = 20000
  bi = 5000
  # BRR
  set.seed(1)
  ETA_RR = list(list(X=myGDrr, model = 'BRR'))
  RR = BGLR(y = yNA,ETA = ETA_RR,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi ) 
  #BRR+fixed effect of Qsd1
  set.seed(1)
  ETA_RR_qsd1Fixed = list(list(X=Qsd1, model = 'FIXED'),list(X=myGDrr, model = 'BRR'))
  RR_qsd1_fixed = BGLR(y = yNA,ETA = ETA_RR_qsd1Fixed,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi) 
  #LASSO
  set.seed(1)
  ETA_LASSO = list(list(X = myGDrr, model = 'BL'))
  BL = BGLR(y = yNA,ETA = ETA_LASSO,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi) 
  # LASSE with Qsd1
  set.seed(1)
  ETA_LASSO_qsd1Fixed = list(list(X = Qsd1, model = "FIXED"), list(X = myGDrr, model = 'BL'))
  BL_qsd1_fixed = BGLR(y = yNA,ETA = ETA_LASSO_qsd1Fixed,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi) 
  #BayesB
  set.seed(1)
  ETA_BayesB = list(list(X = myGDrr, model = 'BayesB'))
  BayesB = BGLR(y = yNA, ETA = ETA_BayesB,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi)
  #BayesB with Qsd1 Fixed
  set.seed(1)
  ETA_BayesB_Qsd1Fixed = list(list(X = Qsd1, model = "FIXED"),list(X = myGDrr, model = 'BayesB'))
  BayesB_Qsd1Fixed = BGLR(y = yNA, ETA = ETA_BayesB_Qsd1Fixed,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi)
  #BayesC
  set.seed(1)
  ETA_BayesC = list(list(X = myGDrr, model = 'BayesC'))
  BayesC = BGLR(y = yNA, ETA = ETA_BayesC,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi)
  #BayesC with Qsd1 Fixed
  set.seed(1)
  ETA_BayesC_Qsd1Fixed = list(list(X = Qsd1, model = "FIXED"),list(X = myGDrr, model = 'BayesC'))
  BayesC_Qsd1Fixed = BGLR(y = yNA, ETA = ETA_BayesC_Qsd1Fixed,saveAt ='Analysis/BGLROutput/BGLRtemp/',nIter =inter,burnIn=bi )
  
  results = data.frame(Family= 'Overall',True = yPheno[tst], Qsd1Status = 'AllAlleles',
                       BRR = RR$yHat[tst], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[tst],
                       LASSO = BL$yHat[tst], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[tst],
                       BayesB = BayesB$yHat[tst], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[tst],
                       BayesC = BayesC$yHat[tst], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[tst]) %>%
    rbind(data.frame(Family= 'Overall',True = yPheno[tst], Qsd1Status = Qsd1[tst],
                     BRR = RR$yHat[tst], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[tst],
                     LASSO = BL$yHat[tst], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[tst],
                     BayesB = BayesB$yHat[tst], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[tst],
                     BayesC = BayesC$yHat[tst], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[tst])) %>%
    rbind(data.frame(Family= WinterGPdata$Family[tst],True = yPheno[tst],Qsd1Status = 'AllAlleles',
                     BRR = RR$yHat[tst], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[tst],
                     LASSO = BL$yHat[tst], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[tst],
                     BayesB = BayesB$yHat[tst], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[tst],
                     BayesC = BayesC$yHat[tst], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[tst]))%>%
    rbind(data.frame(Family= WinterGPdata$Family[tst],True = yPheno[tst],Qsd1Status = Qsd1[tst],
                     BRR = RR$yHat[tst], BRR_Qsd1_Fixed = RR_qsd1_fixed$yHat[tst],
                     LASSO = BL$yHat[tst], LASSO_Qsd1_Fixed = BL_qsd1_fixed$yHat[tst],
                     BayesB = BayesB$yHat[tst], BayesB_Qsd1_Fixed = BayesB_Qsd1Fixed$yHat[tst],
                     BayesC = BayesC$yHat[tst], BayesC_Qsd1_Fixed = BayesC_Qsd1Fixed$yHat[tst])) %>%
    pivot_longer(cols = !c('Family','True', Qsd1Status), names_to = 'Model', values_to = 'Predicted') %>%
    filter(Family!='DH130910') %>%
    group_by(Family, Model, Qsd1Status) %>%
    dplyr::summarise(correlation = cor(True,Predicted), n_test = n()) %>% ungroup() %>%
    join(.,as.data.frame(base::table(WinterGPdata$Family[tst])) %>% dplyr::rename(Family = Var1, n_train = Freq) %>% filter(Family!='DH130910') %>%
           add_row(Family = 'Overall',n_train = length(423 - tst)))
  
  return(results)
}
#####
load('Analysis/AllDHBluesPerYear.RData')
# Change this dataset and rerun all code below to get each yearly thing. 
WinterGPdata = AllDHBluesPerYear %>% filter(TP %in% c('TP1','TP2','TP3','TP4','TP5') & year == '2020/2021' | TP %in% c('TP1','TP2') & year == '2021') %>% 
  filter(type == 'BLUE',Family != 'Cha', Family != 'End') %>% 
  group_by(taxa) %>% add_tally() %>% ungroup() %>% filter(n == 14) %>%
  filter(taxa %in% WinterGD$taxa) %>% arrange(year, taxa) %>% join(WinterGD[,c('taxa','Qsd1')]) %>% 
  filter(Qsd1 != 1) %>% select(!c(PM_date,type,n)) %>% mutate(year =ifelse(year=='2020/2021','2020_2021', '2021')) %>%
  pivot_wider(names_from = c(TP,trait,year), values_from = value) %>% arrange(taxa)

DH_GD_GP_asMarkers = WinterGD %>% filter(taxa %in% WinterGPdata$taxa) %>% arrange(taxa)
sum(WinterGPdata$taxa==DH_GD_GP_asMarkers$taxa)
DH_GD_GP_asMarkersrr = DH_GD_GP_asMarkers[,-1]-1
######
set.seed(1)
trnNumbers = c(70,100,150,200,250,300,350)
CVsets = list()
for (i in 1:7){
  CVsets[[i]] = replicate(sort(sample(1:dim(WinterGPdata)[1],trnNumbers[i],replace = F)),simplify = F, n = 60)  
}
CVsets[[1]][[1]]
library(BGLR)
library(foreach)
library(parallel)
#TO Run on Mac:
cl = parallel::makeCluster(20)
doParallel::registerDoParallel(cl)
setwd(rprojroot::find_rstudio_root_file())

###### Finished with  Models PA framework #######
# Now to the structured GBLUP ######
#Full population for the structured GBLUP. 

N_full = WinterGPdata$taxa %>% unique() %>% length()
DH_GD_GP = WinterGD %>% filter(taxa %in% WinterGPdata$taxa) %>% arrange(taxa) %>% remove_rownames() %>%
  column_to_rownames(var = 'taxa')

# Genotype matrixes that are relevent: ######
library(rrBLUP) 
DHFullPop = A.mat(DH_GD_GP[,]-1) 
DHFullPop <- DHFullPop[order(row.names(DHFullPop)),]
Qsd1Groups = DH_GD_GP$Qsd1
FamilyGroups = WinterGPdata %>% mutate(Family = ifelse(Family =='DH130910','Scala/DH130910',Family)) %>% pull(Family)

# The Qsd1 Specific Matrices 
FullQsd1_0 = A.mat(DH_GD_GP[which(DH_GD_GP$Qsd1==0),]-1)
FullQsd1_2 = A.mat(DH_GD_GP[which(DH_GD_GP$Qsd1==2),]-1)

Full_byQsd1K = matrix(bdiag(FullQsd1_0,FullQsd1_2),nrow = N_full, ncol = N_full, 
                      dimnames = list(c(rownames(FullQsd1_0),rownames(FullQsd1_2)),
                                      c(rownames(FullQsd1_0),rownames(FullQsd1_2))))
Full_byQsd1K <- Full_byQsd1K[ order(row.names(Full_byQsd1K)), order(colnames(Full_byQsd1K)) ]
sum(rownames(DH_GD_GP)== rownames(Full_byQsd1K))
Full_byQsd1K_2 = Full_byQsd1K;
Full_byQsd1K_2[Qsd1Groups !=2,] <- 0
Full_byQsd1K_0 = Full_byQsd1K;
Full_byQsd1K_0[Qsd1Groups !=0,] <- 0

# now need the family specfic matrices ######
FullFamFlavia = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Flavia/DH130910'),]-1)
FullFamScala = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Scala/DH130910'|WinterGPdata$Family=='DH130910'),]-1)
FullFamSY_Tepee = A.mat(DH_GD_GP[which(WinterGPdata$Family=='SY_Tepee/DH130910'),]-1)
FullFamWintmalt = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Wintmalt/DH130910'),]-1)

Full_byFam <- matrix(bdiag(FullFamFlavia,FullFamScala,FullFamSY_Tepee,FullFamWintmalt), nrow = N_full, ncol =N_full,
                     dimnames = list(c(row.names(FullFamFlavia),row.names(FullFamScala),
                                       row.names(FullFamSY_Tepee),row.names(FullFamWintmalt)),
                                     c(row.names(FullFamFlavia),row.names(FullFamScala),
                                       row.names(FullFamSY_Tepee),row.names(FullFamWintmalt))))
Full_byFam <- Full_byFam[order(row.names(Full_byFam)),order(colnames(Full_byFam))]

Full_byFam_Flavia = Full_byFam ;Full_byFam_Flavia[FamilyGroups !='Flavia/DH130910',] <- 0 
Full_byFam_Scala = Full_byFam ;Full_byFam_Scala[FamilyGroups !='Scala/DH130910',] <- 0 
Full_byFam_SY_Tepee = Full_byFam;Full_byFam_SY_Tepee[FamilyGroups !='SY_Tepee/DH130910',] <- 0 
Full_byFam_Wintmalt = Full_byFam;Full_byFam_Wintmalt[(FamilyGroups !='Wintmalt/DH130910'),] <- 0 

# Lastly the Qsd1:Family specific ######
Flavia_2 = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Flavia/DH130910'),]-1) # this will be the same as flavia
Scala_0  = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Scala/DH130910'&WinterGPdata$Qsd1==0|
                                  WinterGPdata$Family=='DH130910'),]-1) 
Scala_2  = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Scala/DH130910'&WinterGPdata$Qsd1==2),]-1) 
Tepee_0 = A.mat(DH_GD_GP[which(WinterGPdata$Family=='SY_Tepee/DH130910'&WinterGPdata$Qsd1==0),]-1) 
Tepee_2 = A.mat(DH_GD_GP[which(WinterGPdata$Family=='SY_Tepee/DH130910'&WinterGPdata$Qsd1==2),]-1) 
Wintmalt_0 = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Wintmalt/DH130910'&WinterGPdata$Qsd1==0),]-1) 
Wintmalt_2 = A.mat(DH_GD_GP[which(WinterGPdata$Family=='Wintmalt/DH130910'&WinterGPdata$Qsd1==2),]-1) 

Full_byQsd1AndFam = matrix(bdiag(Flavia_2,Scala_0,Scala_2,Tepee_0,Tepee_2,Wintmalt_0,Wintmalt_2)+diag(0.001, nrow = N_full), 
                           nrow = N_full,ncol = N_full, 
                           dimnames = list(c(row.names(Flavia_2),row.names(Scala_0),row.names(Scala_2),row.names(Tepee_0),
                                             row.names(Tepee_2),row.names(Wintmalt_0),row.names(Wintmalt_2)),
                                           c(row.names(Flavia_2),row.names(Scala_0),row.names(Scala_2),row.names(Tepee_0),
                                             row.names(Tepee_2),row.names(Wintmalt_0),row.names(Wintmalt_2))))
Full_byQsd1AndFam  = Full_byQsd1AndFam[order(row.names(Full_byQsd1AndFam)),order(colnames(Full_byQsd1AndFam))]

Full_Flavia2 = Full_byQsd1AndFam; Full_Flavia2[FamilyGroups!= 'Flavia/DH130910' & Qsd1Groups != 2,] <- 0
Full_Scala2 = Full_byQsd1AndFam; Full_Scala2[FamilyGroups != 'Scala/DH130910' & Qsd1Groups != 2,] <- 0
Full_Scala0 = Full_byQsd1AndFam; Full_Scala0[FamilyGroups!= 'Scala/DH130910' & Qsd1Groups != 0,] <- 0
Full_SY_Tepee2 = Full_byQsd1AndFam; Full_SY_Tepee2[FamilyGroups!= 'SY_Tepee/DH130910' & Qsd1Groups != 2,] <- 0
Full_SY_Tepee0 = Full_byQsd1AndFam; Full_SY_Tepee0[FamilyGroups!= 'SY_Tepee/DH130910' & Qsd1Groups != 0,] <- 0
Full_Wintmalt2 = Full_byQsd1AndFam; Full_Wintmalt2[FamilyGroups!= 'Wintmalt/DH130910' & Qsd1Groups != 2,] <- 0
Full_Wintmalt0 = Full_byQsd1AndFam; Full_Wintmalt0[FamilyGroups!= 'Wintmalt/DH130910' & Qsd1Groups != 0,] <- 0

# Check Order for all G matrices ######
sum(row.names(DHFullPop)==row.names(Full_byQsd1K))
sum(row.names(DHFullPop)==row.names(Full_byFam))
sum(row.names(DHFullPop)==row.names(Full_byQsd1AndFam))
sum(WinterGPdata$taxa ==row.names(DHFullPop))

# Now put together all the FULL ETA for BGLR. 

Base.eta = list(list(K = DHFullPop, model = 'RKHS'))

Qsd1_full.eta =list(list(K = DHFullPop, model = 'RKHS'),
                    list(X = model.matrix(~as.factor(Qsd1Groups)), model = 'FIXED'),
                    list(K = Full_byQsd1K_0, model = 'RKHS'),
                    list(K = Full_byQsd1K_2, model = 'RKHS'))
Qsd1_BaseQsd1Fixeff.eta = list(list(X = model.matrix(~as.factor(Qsd1Groups)), model = 'FIXED'),
                               list(K = DHFullPop, model = 'RKHS'))
Qsd1_NoFixeff.eta =list(list(K = DHFullPop, model = 'RKHS'),
                        list(K = Full_byQsd1K_0, model = 'RKHS'),
                        list(K = Full_byQsd1K_2, model = 'RKHS'))
Qsd1_JustCov.eta = list(list(K = Full_byQsd1K_0, model = 'RKHS'),
                        list(K = Full_byQsd1K_2, model = 'RKHS'))
Qsd1_JustCovWithFixEff.eta = list(list(X = model.matrix(~as.factor(Qsd1Groups)), model = 'FIXED'),
                                  list(K = Full_byQsd1K_0, model = 'RKHS'),
                                  list(K = Full_byQsd1K_2, model = 'RKHS'))

Fam_full.eta =list(list(K = DHFullPop, model = 'RKHS'),
                   list(X = model.matrix(~factor(FamilyGroups)), model = 'FIXED'),
                   list(K = Full_byFam_Flavia, model = 'RKHS'),
                   list(K = Full_byFam_Scala, model = 'RKHS'),
                   list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                   list(K = Full_byFam_Wintmalt, model = 'RKHS'))
Fam_BaseFamFixeff.eta = list(list(X = model.matrix(~as.factor(FamilyGroups)), model = 'FIXED'),
                             list(K = DHFullPop, model = 'RKHS'))
Fam_NoFixeff.eta =list(list(K = DHFullPop, model = 'RKHS'),
                       list(K = Full_byFam_Flavia, model = 'RKHS'),
                       list(K = Full_byFam_Scala, model = 'RKHS'),
                       list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                       list(K = Full_byFam_Wintmalt, model = 'RKHS'))
Fam_JustCov.eta =list(list(K = Full_byFam_Flavia, model = 'RKHS'),
                      list(K = Full_byFam_Scala, model = 'RKHS'),
                      list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                      list(K = Full_byFam_Wintmalt, model = 'RKHS'))
Fam_JustCovWithFixEff.eta =list(list(X = model.matrix(~factor(FamilyGroups)), model = 'FIXED'),
                                list(K = Full_byFam_Flavia, model = 'RKHS'),
                                list(K = Full_byFam_Scala, model = 'RKHS'),
                                list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                                list(K = Full_byFam_Wintmalt, model = 'RKHS'))

Qsd1PlusFam_Full.eta =list(list(K = DHFullPop, model = 'RKHS'), 
                           list(X = model.matrix(~as.factor(Qsd1Groups)+FamilyGroups), model = 'FIXED'),
                           list(K = Full_byQsd1K_0, model = 'RKHS'),
                           list(K = Full_byQsd1K_2, model = 'RKHS'),
                           list(K = Full_byFam_Flavia, model = 'RKHS'),
                           list(K = Full_byFam_Scala, model = 'RKHS'),
                           list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                           list(K = Full_byFam_Wintmalt, model = 'RKHS'))
Qsd1PlusFam_BaseQPFFixeff.eta = list(list(X = model.matrix(~as.factor(FamilyGroups)), model = 'FIXED'),
                                     list(X = model.matrix(~as.factor(Qsd1Groups)), model = 'FIXED'),
                                     list(K = DHFullPop, model = 'RKHS'))
Qsd1PlusFam_NoFixeff.eta =list(list(K = DHFullPop, model = 'RKHS'), 
                               list(K = Full_byQsd1K_0, model = 'RKHS'),
                               list(K = Full_byQsd1K_2, model = 'RKHS'),
                               list(K = Full_byFam_Flavia, model = 'RKHS'),
                               list(K = Full_byFam_Scala, model = 'RKHS'),
                               list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                               list(K = Full_byFam_Wintmalt, model = 'RKHS'))
Qsd1PlusFam_JustCov.eta =list( list(K = Full_byQsd1K_0, model = 'RKHS'),
                               list(K = Full_byQsd1K_2, model = 'RKHS'),
                               list(K = Full_byFam_Flavia, model = 'RKHS'),
                               list(K = Full_byFam_Scala, model = 'RKHS'),
                               list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                               list(K = Full_byFam_Wintmalt, model = 'RKHS'))
Qsd1PlusFam_JustCovWithFixEff.eta =list(list(X = model.matrix(~as.factor(FamilyGroups)), model = 'FIXED'),
                                        list(X = model.matrix(~as.factor(Qsd1Groups)), model = 'FIXED'),
                                        list(K = Full_byQsd1K_0, model = 'RKHS'),
                                        list(K = Full_byQsd1K_2, model = 'RKHS'),
                                        list(K = Full_byFam_Flavia, model = 'RKHS'),
                                        list(K = Full_byFam_Scala, model = 'RKHS'),
                                        list(K = Full_byFam_SY_Tepee, model = 'RKHS'),
                                        list(K = Full_byFam_Wintmalt, model = 'RKHS'))


Qsd1byFam_full.eta =list(list(K = DHFullPop, model = 'RKHS'), 
                         list(X = model.matrix(~as.factor(Qsd1Groups):FamilyGroups), model = 'FIXED'),
                         list(K=Full_Scala2+diag(x = .001, nrow = 424), model = 'RKHS'),
                         list(K=Full_Scala0+diag(x = .001, nrow = 424), model = 'RKHS'),
                         list(K=Full_Wintmalt2+diag(x = .001, nrow = 424), model = 'RKHS'),
                         list(K=Full_Wintmalt0+diag(x = .001, nrow = 424), model = 'RKHS'),  
                         list(K=Full_SY_Tepee0+diag(x = .001, nrow = 424), model = 'RKHS'),
                         list(K=Full_SY_Tepee2+diag(x = .001, nrow = 424), model = 'RKHS'),
                         list(K=Full_Flavia2+diag(x = .001, nrow = 424), model = 'RKHS'))
Qsd1byFam_BaseQsd1byFamFixeff.eta = list(list(X = model.matrix(~as.factor(Qsd1Groups):FamilyGroups), model = 'FIXED'),
                                         list(K = DHFullPop, model = 'RKHS'))
Qsd1byFam_NoFixeff.eta =list(list(K = DHFullPop, model = 'RKHS'), 
                             list(K=Full_Scala2 +diag(x = .001, nrow = 424), model = 'RKHS'),
                             list(K=Full_Scala0 +diag(x = .001, nrow = 424), model = 'RKHS'),
                             list(K=Full_Wintmalt2+diag(x = .001, nrow = 424), model = 'RKHS'),
                             list(K=Full_Wintmalt0+diag(x = .001, nrow = 424), model = 'RKHS'),  
                             list(K=Full_SY_Tepee0+diag(x = .001, nrow = 424), model = 'RKHS'),
                             list(K=Full_SY_Tepee2+diag(x = .001, nrow = 424), model = 'RKHS'),
                             list(K=Full_Flavia2+diag(x = .001, nrow = 424), model = 'RKHS'))
Qsd1byFam_JustCov.eta =list(list(K=Full_Scala2+diag(x = .001, nrow = 424) , model = 'RKHS'),
                            list(K=Full_Scala0 +diag(x = .001, nrow = 424), model = 'RKHS'),
                            list(K=Full_Wintmalt2+diag(x = .001, nrow = 424), model = 'RKHS'),
                            list(K=Full_Wintmalt0+diag(x = .001, nrow = 424), model = 'RKHS'),  
                            list(K=Full_SY_Tepee0+diag(x = .001, nrow = 424), model = 'RKHS'),
                            list(K=Full_SY_Tepee2+diag(x = .001, nrow = 424), model = 'RKHS'),
                            list(K=Full_Flavia2+diag(x = .001, nrow = 424), model = 'RKHS'))
Qsd1byFam_JustCovWithFixEff.eta =list(list(X = model.matrix(~as.factor(Qsd1Groups):FamilyGroups), model = 'FIXED'),
                                      list(K=Full_Scala2+diag(x = .001, nrow = 424) , model = 'RKHS'),
                                      list(K=Full_Scala0+diag(x = .001, nrow = 424) , model = 'RKHS'),
                                      list(K=Full_Wintmalt2+diag(x = .001, nrow = 424), model = 'RKHS'),
                                      list(K=Full_Wintmalt0+diag(x = .001, nrow = 424), model = 'RKHS'),  
                                      list(K=Full_SY_Tepee0+diag(x = .001, nrow = 424), model = 'RKHS'),
                                      list(K=Full_SY_Tepee2+diag(x = .001, nrow = 424), model = 'RKHS'),
                                      list(K=Full_Flavia2+diag(x = .001, nrow = 424), model = 'RKHS'))



Full.eta = list(Base = Base.eta,
                Qsd1_full.eta = Qsd1_full.eta,
                Qsd1_BaseQsd1Fixeff.eta=Qsd1_BaseQsd1Fixeff.eta,
                Qsd1_NoFixeff.eta = Qsd1_NoFixeff.eta,
                Qsd1_JC.eta = Qsd1_JustCov.eta,
                Qsd1_JCWithFixEff.eta= Qsd1_JustCovWithFixEff.eta,
                
                Fam_full.eta = Fam_full.eta,
                Fam_BaseFamFixeff.eta= Fam_BaseFamFixeff.eta,
                Fam_NoFixeff.eta = Fam_NoFixeff.eta,
                Fam_JC.eta= Fam_JustCov.eta,
                Fam_JCWithFixEff.eta=Fam_JustCovWithFixEff.eta,
                
                Qsd1PlusFam_Full.eta = Qsd1PlusFam_Full.eta,
                Qsd1PlusFam_BaseQPFFixeff.eta = Qsd1PlusFam_BaseQPFFixeff.eta,
                Qsd1PlusFam_NoFixeff.eta=Qsd1PlusFam_NoFixeff.eta,
                Qsd1PlusFam_JC.eta = Qsd1PlusFam_JustCov.eta,
                Qsd1PlusFam_JCWithFixEff.eta = Qsd1PlusFam_JustCovWithFixEff.eta,
                
                Qsd1byFam_full.eta=Qsd1byFam_full.eta,
                Qsd1byFam_BaseQsd1byFamFixeff.eta = Qsd1byFam_BaseQsd1byFamFixeff.eta,
                Qsd1byFam_NoFixeff.eta = Qsd1byFam_NoFixeff.eta,
                Qsd1byFam_JC.eta = Qsd1byFam_JustCov.eta,
                Qsd1byFam_JCWithFixEff.eta = Qsd1byFam_JustCovWithFixEff.eta)


Full.etaList = c('Base','Qsd1_full.eta','Qsd1_BaseQsd1Fixeff.eta','Qsd1_NoFixeff.eta','Qsd1_JC.eta','Qsd1_JCWithFixEff.eta',
                 'Fam_full.eta','Fam_BaseFamFixeff.eta',  'Fam_NoFixeff.eta','Fam_JC.eta','Fam_JCWithFixEff.eta',
                 'Qsd1PlusFam_Full.eta','Qsd1PlusFam_BaseQPFFixeff.eta','Qsd1PlusFam_NoFixeff.eta','Qsd1PlusFam_JC.eta','Qsd1PlusFam_JCWithFixEff.eta',
                 'Qsd1byFam_full.eta','Qsd1byFam_BaseQsd1byFamFixeff.eta','Qsd1byFam_NoFixeff.eta','Qsd1byFam_JC.eta','Qsd1byFam_JCWithFixEff.eta')

'Qsd1byFam_full.eta'

c('Base','byFam','byQsd1','byQsd1PlusFam','byQsd1Fam')
n_train = 325
set.seed(1)
CVsetsFullPop = replicate(60,sort(sample(1:N_full,N_full-n_train,replace = F)), simplify = F)

niter = 20000
bi = 5000

for (j in 4:17)  {
  
  tempFileName <- paste0("Analysis/BGLROutput/WinterStrucGBLUP_", colnames(WinterGPdata)[j], ".csv")
  tempBGLR = "Analysis/BGLROutput/BGLRtemp/"
  
  tempYFull = pull(WinterGPdata[,j])
  Results = data.frame()
  for (ETA in Full.etaList){
    temp <- foreach (i = 1:60, .combine=rbind) %dopar% {
      import::from(magrittr, "%>%")
      import::from(dplyr, "group_by")
      y.na <- tempYFull
      tst <- CVsetsFullPop[[i]];
      y.na[tst] <- NA
      set.seed(1)
      ModelCall <-BGLR::BGLR(y= y.na,ETA=Full.eta[[ETA]],nIter=niter, burnIn=bi,saveAt=paste0(tempBGLR, i))
      data.frame(Pred = ModelCall$yHat[tst], True = tempYFull[tst], Family= 'Overall', Qsd1 = 'Overall') %>%
        rbind(data.frame(Pred = ModelCall$yHat[tst], True = tempYFull[tst], Family= WinterGPdata$Family[tst],Qsd1 = 'Overall'),
              data.frame(Pred = ModelCall$yHat[tst], True = tempYFull[tst], Family= WinterGPdata$Family[tst],Qsd1 = WinterGPdata$Qsd1[tst]),
              data.frame(Pred = ModelCall$yHat[tst], True = tempYFull[tst], Family= 'Overall',Qsd1 = WinterGPdata$Qsd1[tst])) %>%
        group_by(Family, Qsd1) %>% dplyr::summarize(correlation = cor(True,Pred))}
    Results = temp %>% mutate(ModelCalled = ETA, trait = colnames(WinterGPdata)[j], trainingPop = 'Full') %>% rbind(Results)
    print(colnames(WinterGPdata)[j])
    print(ETA)
    print('Finished with those above')
  }
  write.csv(Results,tempFileName)
}
