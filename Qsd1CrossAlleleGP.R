library(sommer);library(arm);library(lme4);library(Hmisc);library(readxl);
library(tibble);library(patchwork);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr)
library(BGLR)
library(tidyverse)
library(multidplyr)
library(dplyr);library(plyr)

load('Analysis/AllDHBluesPerYear.RData')

FoldCVGP_tidy_randomTrainS = function(df, groupvars= .y, mGD=WinterGD, numfolds = 100,
                                      OutGroup = F, OutgroupInfo = WinterGPdataNoFlavia){
  Phenotype = df %>% filter(taxa %in% mGD$taxa) %>% arrange(taxa) 
  myGD = mGD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  dim(Phenotype)
  dim(myGD)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  print(groupvars)
  TrueandPredicted = data.frame()
  Qsd1 = myGD$Qsd1
  start_time <- Sys.time() 
  set.seed(1)
  myGD = myGD[,-1]-1 
  phenotype = as.vector(Phenotype$value)
  num_entries = nrow(Phenotype)
  
  if ((OutGroup)){
    OutGroupdf = OutgroupInfo %>% filter(TP == groupvars$TP,
                                         trait ==groupvars$trait,
                                         Qsd1 != groupvars$Qsd1Train ) %>%
      arrange(taxa)
    outGroupGD = mGD %>% filter(taxa%in%OutGroupdf$taxa) %>% arrange(taxa)
    outGroupGD = outGroupGD[,-1]-1 
  }
  
  for (iii in c(40,70,100,150,200,250,300,350,0.5,0.75,0.9)){
    if(iii<1){
      Num_train = round(num_entries*iii,0)
    } else{
      Num_train = iii
    }
    
    if( num_entries < iii){
      next
    }
    
    
    for (i in 1:numfolds){
      
      trainingSet = sort(sample(1:num_entries,Num_train))
      testingSet = setdiff(1:num_entries,trainingSet)
      
      fam_test = Phenotype$Family[testingSet]
      fam_train = Phenotype$Family[trainingSet]
      
      y_train = phenotype[trainingSet] 
      y_test = phenotype[testingSet]
      
      marker_train = myGD[trainingSet,]
      marker_test = myGD[testingSet,]
      
      trained.Model = mixed.solve(y_train, Z = marker_train, K = NULL, SE =FALSE)
      
      PredictedPheno = as.matrix(marker_test) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
      print(cor(y_test,PredictedPheno))
      
      if((OutGroup)){
        OutGroupPrediction = as.matrix(outGroupGD) %*% as.matrix(trained.Model$u)+ as.numeric(trained.Model$beta)
        pred = data.frame(Family = OutGroupdf$Family, Predicted = OutGroupPrediction, TrueValue = OutGroupdf$value,Qsd1 = OutGroupdf$Qsd1) %>%
          rbind(data.frame(Family = "Overall", Predicted = OutGroupPrediction, TrueValue = OutGroupdf$value,Qsd1 = OutGroupdf$Qsd1))  %>%
          rbind(data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
          rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
          filter(Qsd1 != '1', Family != 'DH130910') %>%
          group_by(Family, Qsd1) %>%
          dplyr::summarize(correlation = cor(Predicted, TrueValue),Spearman_Cor = cor(Predicted,TrueValue, method = 'spearman'),n_test = n()) %>%
          ungroup() %>%
          join(as.data.frame(table(Phenotype$Family[trainingSet])) %>%
                 dplyr::rename(Family= Var1, n_train =Freq) %>% filter(Family!='DH130910') %>%
                 add_row(Family = 'Overall',n_train = Num_train)) %>% mutate(fold = i, trainingProportion = iii)
      }else{
        pred = data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = 'AllAlleles') %>% 
          rbind(data.frame(Family = fam_test, Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
          rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = Qsd1[testingSet])) %>%
          rbind(data.frame(Family = 'Overall', Predicted = PredictedPheno, TrueValue = y_test, Qsd1 = 'AllAlleles')) %>%
          filter(Qsd1 != '1', Family != 'DH130910') %>%
          group_by(Family, Qsd1) %>%
          dplyr::summarize(correlation = cor(Predicted, TrueValue), Spearman_Cor = cor(Predicted,TrueValue, method = 'spearman'),n_test = n()) %>%
          ungroup() %>%
          join(as.data.frame(table(Phenotype$Family[trainingSet])) %>%
                 dplyr::rename(Family= Var1, n_train =Freq) %>% filter(Family!='DH130910') %>%
                 add_row(Family = 'Overall',n_train = Num_train)) %>% mutate(fold = i, trainingProportion = iii)
      }
      
      TrueandPredicted = rbind(TrueandPredicted, pred)
    }}
  end_time <- Sys.time()
  print(end_time-start_time)
  return(TrueandPredicted)
}
# we have to compare apples to apples here - so train on the three families that have qsd1 segregating 
# as a whole then split into groups
WinterGPdataNoFlavia = AllDHBluesPerYear %>% filter(year == '2020/2021') %>% 
  filter(type == 'BLUE',Family != 'Cha') %>% filter(Family != 'Flavia/DH130910') %>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1)

NoFlavia_RandonTRN_RandomTST = WinterGPdataNoFlavia %>%  
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% 
  group_by(year, trait, TP) %>% group_modify(FoldCVGP_tidy_randomTrainS)

NoFlavia_inQsd1grp_RandTRN_RandTST_predictsOutgroup = WinterGPdataNoFlavia %>%  
  filter(!(trait =='GE' & TP %in% c('TP4','TP5'))) %>% 
  group_by(year, trait, TP, Qsd1) %>% dplyr::rename(Qsd1Train = Qsd1) %>% 
  group_modify(~FoldCVGP_tidy_randomTrainS(df = .x,groupvars = .y, OutGroup = T))

NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest = AllDHBluesPerYear %>% filter(year == '2021') %>%
  filter(type == 'BLUE',Family != 'Cha',Family != 'Flavia/DH130910') %>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>% filter(TP == 'TP1')%>% 
  group_by(year, trait, TP) %>% group_modify(FoldCVGP_tidy_randomTrainS)

NOFlaviaDH2021_TP1GE_GI_PredictOutgroup = AllDHBluesPerYear %>% filter(year == '2021') %>%
  filter(type == 'BLUE',Family != 'Cha',Family != 'Flavia/DH130910') %>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>% filter(TP == 'TP1')%>% 
  group_by(year, trait, TP, Qsd1) %>% dplyr::rename(Qsd1Train = Qsd1) %>% 
  group_modify(~FoldCVGP_tidy_randomTrainS(df = .x,groupvars = .y, OutGroup = T))

write.csv(NoFlavia_RandonTRN_RandomTST, file = 'Analysis/GP_NoFlavia_rndTest_rndtrain.csv')
write.csv(NoFlavia_inQsd1grp_RandTRN_RandTST_predictsOutgroup, file = 'Analysis/GP_NoFlavia_Perqsd1GroupPre_withOutGroupPredictions.csv')
write.csv(NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest, file = 'Analysis/NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest.csv')
write.csv(NOFlaviaDH2021_TP1GE_GI_PredictOutgroup, file = 'Analysis/NOFlaviaDH2021_TP1GE_GI_PredictOutgroup.csv')

