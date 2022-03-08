library(lme4);library(Hmisc);library(readxl);library(readr)
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x)library(plyr)
library(dplyr)
load('Analysis/AllDHBluesPerYear.RData')
load('GenotypeData/WinterGD.RData')
load('GenotypeData/WinterGM.RData')
WinterGM = WinterGM %>% arrange(Chromosome, Position)
WinterGD = WinterGD %>%
  mutate(taxa = gsub(pattern = '-',replacement = '_',taxa),
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         taxa1 = taxa) %>% remove_rownames()%>% 
  column_to_rownames('taxa1')

load('Analysis/OutputAsreml/FamilyMeansAndVarianceQsd1Asreml.RData')
load('Analysis/OutputAsreml/HeritabilitiesWinterGermTraits.RData')
load('Analysis/OutputAsreml/atFamQsd1TaxaVar.RData')

# Figure 1 Blues over time from combined dataset ####
png('picsPNGforQsd1Effects_paper/BluesByFamilyQsd12020_2021.png', 1400, 800, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Non-dormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020/2021'))  %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  mutate(Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  theme_bw()+
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free')+
  labs(y = 'trait value',x = 'Timepoint', fill = 'Qsd1 Status')
dev.off()

# Figure 2 Heritabilities of the traits over time  #####
HeritabilitiesWinterGermTraits$type %>% unique()

H2Names = c(expression(Cullis~H^2),
            expression(Cullis~H^2~Qsd1~Fixed),
            expression(h^2),
            expression(h^2 ~ Qsd1~Fixed),
            expression(h^2 ~ Qsd1[D]),
            expression(h^2 ~Qsd1[N]))
  
  # "H2_Cullis", "H2_Cullis_Qsd1Accounted", "h2_Narrow","h2_Narrow_Qsd1Accounted", "h2_NarrowU_Qsd1D","h2_NarrowU_Qsd1N",   

png(filename = '/picsPNGforQsd1Effects_paper/Heritabilites.png',
    width = 1000, height = 800, res = 120)

HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge')+theme_bw() +
  scale_fill_manual( values = c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'), labels = H2Names)+
  theme(legend.text.align = 0  )+
  facet_nested('Year'+year~'Trait'+trait) + labs(fill = 'Heritability type', y = 'Heritability',x = 'Timepoint')
dev.off()

HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>% 
  pivot_wider(names_from = 'type',values_from='value') %>% arrange(desc(trait)) %>%
  mutate(ratio = h2_NarrowU_Qsd1N/h2_NarrowU_Qsd1D)

#significance test: does a unique genetic variance per family improve model fit? #####
FamilyMeansAndVarianceQsd1 %>% filter(Coef == 'lrt1') %>% mutate(significant = value< 0.05) 

# Figure 3 Variance attributable to Qsd1  ######
FMV_Qsd1 = FamilyMeansAndVarianceQsd1 %>% ungroup() %>% filter(type =='Fixed') %>% group_by(year, TP,trait, model) %>% 
  group_modify(~{.x%>% mutate(value = .x$value + .x$value[length(.x$value)])}) %>%ungroup() %>%
  rbind(FamilyMeansAndVarianceQsd1 %>% filter(type != 'Fixed')) %>% arrange(year, TP, trait, model)

png('picsPNGforQsd1Effects_paper/VarAttributableToQsd1.png', 1400, 800, res =120)
FMV_Qsd1 %>% filter(model %in% c('base','difference')) %>% filter(trait =='GE') %>%
  filter(substr(Coef, 1,2) %in% c('Fa', 'at')) %>% 
  mutate(Coef = gsub(pattern = "at\\(Family, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef)),
         Coef = gsub('Family_','',Coef)) %>%
  mutate(type = ifelse(type=='Fixed',type,'Variance')) %>%
  filter(type == 'Variance') %>%filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  mutate(model2 = ifelse(model == 'base','Family Genetic Variance\nafter Qsd1 is accounted for',
                         'Variance attributable\nto Qsd1'),
         Coef = gsub(pattern = 'DH130910',replacement = 'Lightning',Coef)) %>%
  ggplot(aes(x = TP, y = value, group = Coef, fill = Coef, alpha = model2))+
  geom_col(position = 'dodge', color = 'black')+facet_nested(year~trait, scales= 'free')+scale_alpha_manual(values = c(.5,1))+
  labs(fill = 'Family', alpha = 'Variance attribution', y = 'Variance', x = 'Timepoint') + theme_bw(base_size = 13) +
  
  FMV_Qsd1 %>% filter(model %in% c('base','difference')) %>% filter(trait =='GI') %>%
  filter(substr(Coef, 1,2) %in% c('Fa', 'at')) %>% 
  mutate(Coef = gsub(pattern = "at\\(Family, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef)),
         Coef = gsub('Family_','',Coef)) %>%
  mutate(type = ifelse(type=='Fixed',type,'Variance')) %>%
  filter(type == 'Variance') %>%filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5')))%>%
  mutate(model2 = ifelse(model == 'base','Family Genetic Variance\nafter Qsd1 is accounted for',
                         'Variance attributable\nto Qsd1'),
         Coef = gsub(pattern = 'DH130910',replacement = 'Lightning',Coef)) %>%
  ggplot(aes(x = TP, y = value, group = Coef, fill = Coef, alpha = model2))+
  geom_col(position = 'dodge', color = 'black')+facet_nested(year~trait, scales= 'free')+scale_alpha_manual(values = c(.5,1))+
  labs(fill = 'Family', alpha = 'Variance attribution',y = 'Variance', x= 'Timepoint')+theme_bw(base_size = 13) +
  plot_layout(ncol = 2, guides = 'collect')
dev.off()
'Colors correspond to families, shading to attribution. The height of the bar indicates the Total genetic variance in the Family 
when Qsd1 is not accounted for, the darker shaded portion of the bar is the genetic variance attributable to Qsd1, while the 
lighter portion is the genetic variance remaining within each family after Qsd1 is accounted for as a fixed effect.'

# Figure 4 Different variances per Family:Qsd1 grouping######
png('picsPNGforQsd1Effects_paper/DifferentVarForFamilyQsd1.png', 1400, 800, res =120)
atFamQsd1TaxaVar %>% filter(substr(Coef,1,2)=='at' ) %>%
  mutate(Coef = gsub(pattern = "at\\(FamilyQsd1, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef))) %>%
  separate(col = Coef, into = c('Family','Qsd1'), sep = '130910',remove = F) %>% 
  mutate(Family = paste0(Family, '130910'), Qsd1lab = ifelse(Qsd1==2, 'Dormant','Non-dormant'),
         Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(TP, value, group = Coef, alpha = Qsd1lab,fill= Family))+
  geom_bar(position = 'dodge',stat= 'identity')+
  facet_nested(trait+year~Family, scales= 'free_y') +scale_alpha_manual(values = c(.5,1)) +theme_bw()+
  geom_text(inherit.aes = F, 
            data = atFamQsd1TaxaVar %>% filter(Coef == 'lrt'|substr(Coef,1,2)=='at') %>% 
              mutate(siglab = ifelse(value<0.01,'','ns'), ypos = max(value[1:7])*1.1) %>% filter(Coef == 'lrt',year != '2020/2021'),
            aes(x = TP, y = ypos, label = siglab))+
  labs(alpha = 'Qsd1 status', y = 'Variance', x = 'Timepoint')
dev.off()

# Figure 5 Correlation of Marker effects between Qsd1 groups #####
GetMarkerEffectsperGroup = function(df, groupvars, GD = WinterGD, skip){
  # I checked these results within the framework of the variable 'skip' to see if they were insensative to 
  # marker numbers. They are.
  Phenotype = df %>% filter(taxa %in% GD$taxa) %>% arrange(taxa) 
  myGD = GD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  Y = Phenotype$value %>% as.vector()
  #lets reduce marker count to make things more clear: maybe by 1/2?
  Msubset = seq(from = 1, to  = dim(myGD)[2]-1, by =skip) #5148 id qsd1.
  myGD = myGD[,-1]-1
  Model.effects = mixed.solve(y = Y, Z = myGD[,Msubset], K = NULL)
  MarkerEffects = Model.effects$u %>% data.frame() %>% rownames_to_column(var = 'Marker') %>% rename(Effect = '.')
  return(MarkerEffects)
}

GetMarkerEffectsperGroupRand = function(df, groupvars, GD = WinterGD, skip=20){
  # I checked these results within the framework of the variable 'skip' to see if they were insensative to 
  # marker numbers. They are.
  
  df.return = data.frame()
  Phenotype = df %>% filter(taxa %in% GD$taxa) %>% arrange(taxa) 
  myGD = GD %>% filter(taxa %in% Phenotype$taxa) %>% arrange(taxa)
  if(!length(unique(Phenotype$taxa))==sum(Phenotype$taxa == myGD$taxa)){
    stop('taxa lists are not correctly aligning')
  }
  Msubset = seq(from = 1, to  = dim(myGD)[2]-1, by =skip) #5148 id qsd1.
  subset.mark = WinterGM[Msubset,] %>% filter(Chromosome != '5') 
  
  Y = Phenotype$value %>% as.vector()
  n = length(Y)
  for (i in subset.mark$SNP) {
    print(sum(myGD[,i])/n/2)
    print(i)
    if (sum(myGD[,i])/n/2 < 0.3 | sum(GD[,i])/n/2 > 0.7){
      next
    }
    if (sum(round(myGD[,i])==1)>.1*n){
      next
    }
    
    group_1 = which(round(myGD[,i]) == 2)
    group_2 = which(round(myGD[,i]) != 2)
    
    y_1 = Y[group_1]
    y_2 = Y[group_2]
    
    M_1 = myGD[group_1,-1]-1
    M_2 = myGD[group_2,-1]-1
    
    model.effect_1 =  mixed.solve(y = y_1, Z = M_1[,Msubset], K = NULL)
    model.effect_2 = mixed.solve(y = y_2, Z = M_2[,Msubset], K = NULL)
    
    df.return = rbind(df.return, data.frame(cor = cor(model.effect_1$u,model.effect_2$u)))
    
  }
  return(df.return %>% tibble())
}

# RandomSNPsSorting = WinterMarkerEffectTesters %>% group_by(TP,trait,year) %>%
#   group_modify(~GetMarkerEffectsperGroupRand(df = .))
# RandomSNPsSorting %>% summarize(mean = mean(cor), sd = sd(cor))
# write.csv(RandomSNPsSorting,'WinterBarley/Analysis/OutputAsreml/RandomSNPASE.csv')

load('Analysis/AllDHBluesPerYear.RData')
WinterMarkerEffectTesters = AllDHBluesPerYear %>% filter(year == '2021') %>% 
  filter(type == 'BLUE',Family != 'Cha') %>% filter(Family != 'Flavia/DH130910') %>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) #%>% mutate(Qsd1 = ifelse(Family == 'Flavia/DH130910','Flavia',Qsd1))
df = WinterMarkerEffectTesters %>% filter(TP =='TP2', trait == 'GI')
WinterMarkerEffectTesters %>% filter(TP =='TP2',trait == 'GI') %>% select(Qsd1) %>% table()
#145 ND vs 208 Dormant
WinterMarkerEffectByQsd1 = WinterMarkerEffectTesters%>% 
  rbind(WinterMarkerEffectTesters %>% mutate(Qsd1 ='Overall')) %>%
  group_by(trait,TP,year,Qsd1)%>%
  group_modify(~GetMarkerEffectsperGroup(df = .,skip =50))

WinterMarkerEffectbyFamily = WinterMarkerEffectTesters %>% group_by(trait,TP, Family) %>%
  group_modify(~GetMarkerEffectsperGroup(df = .,skip =50))

RandomSNPsSorting = read_csv('Analysis/OutputAsreml/RandomSNPASE.csv')
RandomSNPsSorting %>% group_by(TP,trait, year) %>% summarize(mean = mean(cor), sd = sd(cor)) %>%
  filter(trait == 'GI') %>% rename(TPcommon = TP) 

png('picsPNGforQsd1Effects_paper/CorOfQsd1GroupMarkerEffects.png', 800, 500, res =120)
WinterMarkerEffectByQsd1 %>% filter(Qsd1 != 'Overall') %>% ungroup() %>% 
  mutate(Qsd1 = mapvalues(Qsd1,from = c('0','2','Overall'), to = c('ND','D','Overall'))) %>%
  pivot_wider(names_from = c(Qsd1,TP), values_from = Effect) %>%
  ungroup() %>% select(!c(year,Marker)) %>% group_by(trait) %>%
  group_modify(~{
    if (.y$trait == 'GE'){
    cort = cor(.x, use = 'complete.obs') %>% data.frame() 
    b = matrix( nrow = dim(cort)[1],ncol =dim(cort)[1] )
    b[upper.tri(b, diag = T)] = cort[upper.tri(cort,diag = T)]
    
    b %>% as.data.frame() 
    colnames(b) = colnames(cort)
    rownames(b) = rownames(cort)
    b %>% as.data.frame %>% rownames_to_column('TP_Qsd1') %>% pivot_longer(!TP_Qsd1) %>%
      filter(!is.na(value)| str_sub(TP_Qsd1,-3,-1)==str_sub(name,-3,-1))}
    else{
      cort = cor(.x, use = 'complete.obs') %>% data.frame() 
      b = matrix( nrow = dim(cort)[1],ncol =dim(cort)[1] )
      b[lower.tri(b, diag = F)] = cort[lower.tri(cort)]
      
      b %>% as.data.frame() 
      colnames(b) = colnames(cort)
      rownames(b) = rownames(cort)
      b %>% as.data.frame %>% rownames_to_column('TP_Qsd1') %>% pivot_longer(!TP_Qsd1) %>%
        filter(!is.na(value))} 
    })%>% 
  separate(TP_Qsd1, into = c('Qsd1','TP'), sep = '_') %>% separate(name, sep = '_',c('Qsd11', 'TPP')) %>%
  mutate(TPcommon = ifelse(TP==TPP, TP, NA))%>%
  join(RandomSNPsSorting %>% group_by(TP,trait, year) %>% summarize(mean = mean(cor), sd = sd(cor)) %>%
         rename(TPcommon = TP))%>%
  mutate(tval = pnorm(abs(value-mean)/sd, lower.tail = F)) %>%
  transmute(TP_Qsd1 = paste0(Qsd1, '_' ,TP), name = paste0(Qsd11,'_',TPP), 
            value, tval, label = paste0(round(value,2), ifelse(tval<0.05,'Sig','')), 
            color = ifelse(tval<0.05,'colors','no')) %>%
  mutate(label = gsub('NA','', label))%>% 
  ggplot(aes(x = TP_Qsd1, y = name, fill = value) )+geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5) +geom_text(aes(label = label), size = 2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + labs(fill = 'Correlation')+
  xlab('Qsd1 Allele_Timepoint')+ylab('Qsd1 Allele_Timepoint')
dev.off()

# Figure 6 GI model training and prediction by All Qsd1 or ind Qsd1 ######
NoFlavia_inQsd1grp_RandTRN_RandTST_predictsOutgroup = read_csv('Analysis/GP_NoFlavia_Perqsd1GroupPre_withOutGroupPredictions.csv')
NoFlavia_RandonTRN_RandomTST =read_csv('Analysis/GP_NoFlavia_rndTest_rndtrain.csv')
names = c(expression(h[a]),expression(h[aD]),expression(h[aN]))

png('picsPNGforQsd1Effects_paper/GP_GI_trainQsd1TestQsd1.png', 1400, 1000, res =120)
HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(year == '2020/2021', TP %nin% c('TP1.5','TP2.5','TP3.5'), trait == 'GI')  %>%
  filter(substr(type, 1,2)=='h2', type != 'h2_Narrow_Qsd1Accounted') %>%
  mutate(maxCor = sqrt(value),
         type =   mapvalues(x = type, from = c('h2_Narrow','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D' ), 
                            to = c('sigma[a]','sigma[aN]','sigma[aD]'))) %>% 
  ggplot(aes(type, maxCor, fill = type))+
  geom_bar(stat = 'identity',position = 'dodge')+facet_grid(TP~'Trait: GI',margins = F)+
  ylim(0,1)+ xlab('Term')+
  scale_x_discrete(labels = names)+
  theme_bw()+theme(legend.position = 'none', strip.background.y = element_blank(), strip.text.y = element_blank())+
  ylab('Max expected prediction accuracy')+

NoFlavia_inQsd1grp_RandTRN_RandTST_predictsOutgroup %>% 
  filter(Qsd1 != "AllAlleles", n_test >8) %>%
  group_by(trait,Family,Qsd1Train, TP, trainingProportion, Qsd1) %>% 
  dplyr::summarise(stdev = sd(correlation, na.rm = T),
                   correlation = mean(correlation,na.rm = T), 
                   n_train = mean(n_train), 
                   n_test = mean(n_test)) %>% mutate(type = 'Training_Qsd1') %>%
  dplyr::rename(Qsd1Pred = Qsd1) %>% 
  mutate(Qsd1Train = as.character(Qsd1Train),
         Qsd1Pred = as.character(Qsd1Pred)) %>%
  rbind(NoFlavia_RandonTRN_RandomTST %>% mutate(Qsd1Train = 'AllAlleles') %>%
          dplyr::rename(Qsd1Pred = Qsd1) %>%
          group_by(trait, Family, Qsd1Train, TP, trainingProportion, Qsd1Pred) %>%
          filter(n_test>5) %>%
          dplyr::summarize(stdev = sd(correlation, na.rm = T),
                           correlation = mean(correlation,na.rm = T), 
                           n_train = mean(n_train), 
                           n_test = mean(n_test)) %>% 
          mutate(type = 'Training_Overall'))  %>%
  filter(TP %nin% c('TP1.5','TP2.5','TP3.5'), Family == 'Overall', trait == 'GI') %>%
  mutate(Qsd1Train = mapvalues(x = Qsd1Train, from = c('0','2','AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles')),
         Qsd1Pred = mapvalues(x = Qsd1Pred, from = c('0','2', 'AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles')),
         Trained_Predicted = paste0(Qsd1Train, ' \U2192 ',Qsd1Pred)) %>%
ggplot(aes(x = n_train, y = correlation, color = Trained_Predicted, linetype = Trained_Predicted))+
  # geom_point(aes(shape = Trained_Predicted))+
  facet_grid(TP~ 'Trait: GI', scales = 'free') + theme_bw()+
  xlab('Number of training taxa') +
  labs(color = 'Qsd1 Trained \U2192 Predicted', linetype = 'Qsd1 Trained \U2192 Predicted',shape = 'Qsd1 Trained \U2192 Predicted') +
  geom_smooth(method = 'loess', se =T, size = 1.2)+labs(y = 'Realized prediction accuracy')+
  scale_color_manual(values = c('salmon','salmon','salmon','mediumseagreen','mediumseagreen','cornflowerblue','cornflowerblue'))+
  # scale_linetype_manual(values = c('3131','dotted','solid','dotted','solid','dotted','solid'))+
  scale_linetype_manual(values = c('solid','dotted','3131','solid','dotted','dotted','solid'))+
  ylim(0,1)+
  plot_layout(ncol = 2, design = c('ABBBBB'))
# scale_shape_manual(values = c(0,0,0,2,2,3,3))
dev.off()

# Figure 7 Structured GBLUP results.  ######
filesModel = dir("Analysis/BGLROutput/", full.names = T, pattern = "WinterStrucGBLUP_")
FilesStructuredGBLUP = lapply(filesModel,FUN = read.csv) %>% bind_rows() %>% dplyr::select(!X)%>% 
  separate(col = trait, into = c('TP','trait', 'year1','year2')) %>%
  mutate(year = paste0(year1,'_',year2), year = gsub('NA','',year)) %>%
  separate(col = "ModelCalled", sep = "_", into =c('Model', 'Effects')) %>%
  dplyr::mutate(Effects = ifelse(Effects =='Full.eta','full.eta',Effects),
                Effects = substr(Effects,1,4),
                Effects = ifelse(Effects =='Base','Base_fixef',Effects),
                Effects = ifelse(is.na(Effects), 'Base',Effects)) %>%
  filter(!is.na(correlation), Family == 'Overall') %>%
  mutate(Qsd1 = ifelse(Qsd1 == 'Overall','AllAlleles',Qsd1),
         Effects = mapvalues(x = Effects, from = c('JC.e','JCWi','NoFi','full'), to = c('JC','JC_Fix','NoFix','Full'))) 
FilesStructuredGBLUP %>% head()

FilesStructuredGBLUP %>% filter(trait == 'GI') %>% mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2'), to = c('Non-Dormant','Dormant'))) %>%
  filter(year == '2020_2021') %>%
  ggplot(aes(x = Model, y = correlation,fill = Effects, color = Effects))+
  geom_boxplot()+facet_nested(Qsd1~trait+TP, space = 'free_x')+
  geom_hline(yintercept = 0)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))

Significant_differents = FilesStructuredGBLUP %>%
  group_by(year,trait, TP, Qsd1) %>%
  group_modify(~{
    Base = .x %>% filter(Model == 'Base')
    .x %>% filter(Model != 'Base') %>% group_by(Model, Effects) %>%
      group_modify(~{
        ttest = t.test(x = .$correlation, y = Base$correlation)
        return(data.frame(P.value = ttest$p.value, MeanDiffer = mean(.$correlation)-mean(Base$correlation)))
      })
  })

Significant_differents %>% filter(P.value<0.01/300 & MeanDiffer>0) %>% ungroup() %>% arrange(trait,TP,Model, Effects) %>%
  view()

png('picsPNGforQsd1Effects_paper/StrucutredGBLUPGI_Results.png', 1400, 800, res =120)
FilesStructuredGBLUP %>% filter(trait == 'GI') %>%  filter(year == '2020_2021') %>% filter(Qsd1 == 'AllAlleles') %>%
  mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
         Variation = Effects) %>%
  ggplot(aes(x = Model, y = correlation, fill = Variation, color = Variation))+
  geom_boxplot()+facet_wrap(~TP)+
  geom_hline(yintercept = 0)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  geom_text(data = Significant_differents %>% filter(P.value<0.01/300 & MeanDiffer>0) %>% ungroup() %>% 
              filter(trait == 'GI', year == '2020_2021')%>%filter(Qsd1 == 'AllAlleles') %>%
              mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2','AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
                     label = mapvalues(Effects, from = c('Base_fixef','Full','JC','JC_Fix','NoFix'),
                                       to = c('*         ','  *       ','   *   ','      * ','          *')),
                     Variation = Effects),
            aes(y = 1, label =label, group = Variation))+
  labs(y = 'Prediction Accuracy')
dev.off()


# Supp Figure 1 Various model and the PA of each tests #####
filesModel = dir("Analysis/BGLROutput/", full.names = T, pattern = "Various")
VariousModels = lapply(filesModel,FUN = read.csv) %>% bind_rows() %>% select(!X)
VariousModels %>%head()

png('/picsPNGforQsd1Effects_paper/Sup_VariousModelsPA.png', 1400, 1200, res =120)
VariousModels %>% filter(Family == 'Overall') %>% 
  mutate(Qsd1Status = mapvalues(Qsd1Status, from = c('-1','1','AllAlleles'), to = c('NonDormant','Dormant','Both Alleles'))) %>%
  ggplot(aes(x = Model, y = correlation, fill = Model, color = Qsd1Status))+
  geom_boxplot()+ labs(y = 'Prediction Accuracy')+
  facet_grid(trait+n_train~Family, scales ='free_x')
dev.off()

# Supp Figure 2 All the blues over time from each datatset #####
png('picsPNGforQsd1Effects_paper/Sup_BluesByFamilyQsd1AllYears.png', 2400, 1500, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Non-dormant'),
         Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  theme_bw()+
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free', space = 'free_x')+
  labs(x = 'Timepoint',fill = ' Qsd1 status')
dev.off()

# Supp Figure 3 Family:Qsd1 unique genetic variances as heritability ###### 
png('picsPNGforQsd1Effects_paper/Sup_H2_DifferentVarForFamilyQsd1.png', 1400, 800, res =120)
atFamQsd1TaxaVar %>% filter(substr(Coef,1,2)=='at' | Coef == 'units!units') %>%
  group_modify(~{
    error = .x %>% filter(Coef == 'units!units') %>% select(value) %>% as.numeric() 
    .x = .x %>% filter(substr(Coef,1,2)=='at') %>% mutate(value = value/(value+error))
    return(.x)
  }) %>%
  mutate(Coef = gsub(pattern = "at\\(FamilyQsd1, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef))) %>%
  separate(col = Coef, into = c('Family','Qsd1'), sep = '130910',remove = F) %>% 
  mutate(Family = paste0(Family, '130910'), Qsd1lab = ifelse(Qsd1==2, 'Dormant','Non-dormant'),
          Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(TP, value, group = Coef, alpha = Qsd1lab,fill= Family))+
  geom_bar(position = 'dodge',stat= 'identity')+
  facet_nested(trait+year~Family, scales= 'free_y') +scale_alpha_manual(values = c(.5,1)) +theme_bw()+
  geom_text(inherit.aes = F, 
            data = atFamQsd1TaxaVar %>% filter(Coef == 'lrt'|substr(Coef,1,2)=='at') %>% 
              mutate(siglab = ifelse(value<0.01,'','ns'), ypos = max(value[1:7])*1.1) %>% filter(Coef == 'lrt',year != '2020/2021'),
            aes(x = TP, y = ypos, label = siglab))+
  labs(alpha = 'Qsd1 status', y = 'Heritability', x = 'Timepoint')+
  ylim(0,1)
dev.off()


# Supp Figure 4 Heritability of the remining variance within families after Qsd1 is accounted for #####
# What is the ratio of error variance to the at(Family):taxa geentic variance? ie a heritability after Qsd1 is accounted for?
png(filename = 'picsPNGforQsd1Effects_paper/RemainingFamilyGenVarH2.png',
    width = 1000, height = 800, res = 120)
FamilyMeansAndVarianceQsd1 %>% filter(model =='Family:Qsd1_Fixed', type == 'Variance') %>%
  group_modify(~{
    error = .x %>% filter(Coef == 'units!units') %>% select(value) %>% as.numeric()
    FamilyH2 = .x %>% filter(substr(Coef,1,2)=='at') %>% mutate(H2 = value/(value+error))
    return(FamilyH2)
  }) %>% 
  mutate(Coef = gsub(pattern = "at\\(Family, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef)),
         Coef = gsub('Family_','',Coef),
         Coef = gsub(pattern = 'DH130910',replacement = 'Lightning',Coef)) %>%
  filter(!(year == '2020/2021'& TP %in% c('TP1.5','TP2.5','TP3.5')))%>%
  ggplot(aes(x = TP, y = H2,group = Coef, fill = Coef))+theme_bw()+
  geom_col(position = 'dodge')+facet_nested('Year'+ year~'Trait'+trait, scales= 'free')+
  labs(group = 'Family',fill = 'Family', y = 'Broad Sense H2', x = 'Timepoint')
dev.off()


# Supp Figure 5 GE prediction by Qsd1 overall of ind Qsd1 ######

png('picsPNGforQsd1Effects_paper/Sup_GE_Qsd1trainQsd1Test.png', 1400, 800, res =120)
HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(year == '2020/2021', TP %nin% c('TP1.5','TP2.5','TP3.5','TP4','TP5'), trait == 'GE')  %>%
  filter(substr(type, 1,2)=='h2', type != 'h2_Narrow_Qsd1Accounted') %>%
  mutate(maxCor = sqrt(value),
         type =   mapvalues(x = type, from = c('h2_Narrow','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D' ), 
                            to = c('sigma[a]','sigma[aN]','sigma[aD]'))) %>% 
  ggplot(aes(type, maxCor, fill = type))+
  geom_bar(stat = 'identity',position = 'dodge')+facet_grid(TP~'Trait: GE',margins = F)+
  ylim(0,1)+ xlab('Term')+
  scale_x_discrete(labels = names)+
  theme_bw()+theme(legend.position = 'none', strip.background.y = element_blank(), strip.text.y = element_blank())+
  ylab('Max expected prediction accuracy')+
  
  NoFlavia_inQsd1grp_RandTRN_RandTST_predictsOutgroup %>% 
  filter(Qsd1 != "AllAlleles", n_test >8) %>% filter(Family == 'Overall') %>%
  group_by(trait,Family,Qsd1Train, TP, trainingProportion, Qsd1) %>% 
  dplyr::summarise(stdev = sd(correlation, na.rm = T),
                   correlation = mean(correlation,na.rm = T),
                   n_train = mean(n_train),
                   n_test = mean(n_test)) %>% mutate(type = 'Training_Qsd1') %>%
  dplyr::rename(Qsd1Pred = Qsd1) %>% 
  mutate(Qsd1Train = as.character(Qsd1Train),
         Qsd1Pred = as.character(Qsd1Pred)) %>%
  rbind(NoFlavia_RandonTRN_RandomTST %>% mutate(Qsd1Train = 'AllAlleles') %>%
          dplyr::rename(Qsd1Pred = Qsd1) %>%
          group_by(trait, Family, Qsd1Train, TP, trainingProportion, Qsd1Pred) %>%
          filter(n_test>5) %>%
          dplyr::summarize(stdev = sd(correlation, na.rm = T),
                           correlation = mean(correlation,na.rm = T), 
                           n_train = mean(n_train), 
                           n_test = mean(n_test)) %>% 
          mutate(type = 'Training_Overall'))  %>%
  filter(TP %nin% c('TP1.5','TP2.5','TP3.5'), Family == 'Overall', trait == 'GE') %>%
  mutate(Qsd1Train = mapvalues(x = Qsd1Train, from = c('0','2','AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles')),
         Qsd1Pred = mapvalues(x = Qsd1Pred, from = c('0','2', 'AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles')),
         Trained_Predicted = paste0(Qsd1Train, ' \U2192 ',Qsd1Pred)) %>%
  ggplot(aes(x = n_train, y = correlation, color = Trained_Predicted, linetype = Trained_Predicted))+
  # geom_point(aes(shape = Trained_Predicted))+
  facet_grid(TP~ 'Trait: GE', scales = 'free') + theme_bw()+
  xlab('Number of training taxa') +
  labs(color = 'Qsd1 Trained \U2192 Predicted', linetype = 'Qsd1 Trained \U2192 Predicted',shape = 'Qsd1 Trained \U2192 Predicted') +
  geom_smooth(method = 'loess', se =T, size = 1.2,)+labs(y = 'Realized prediction accuracy')+
  scale_color_manual(values = c('salmon','salmon','salmon','mediumseagreen','mediumseagreen','cornflowerblue','cornflowerblue'))+
  # scale_linetype_manual(values = c('3131','dotted','solid','dotted','solid','dotted','solid'))+
  scale_linetype_manual(values = c('solid','dotted','3131','solid','dotted','dotted','solid'))+
  ylim(0,1)+
  plot_layout(ncol = 2, design = c('ABBBBB'))
dev.off()
# Supp Figure 6 GE Structured GBLUP #####
png('picsPNGforQsd1Effects_paper/Sup_StrucutredGBLUPGE_Results.png', 1400, 800, res =120)
FilesStructuredGBLUP %>% filter(trait == 'GE') %>% mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2','AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
                                                          Variation = Effects) %>%
  ggplot(aes(x = Model, y = correlation, fill = Variation, color = Variation))+
  geom_boxplot()+facet_nested(Qsd1~trait+TP, space = 'free_x')+
  geom_hline(yintercept = 0)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  geom_text(data = Significant_differents %>% filter(P.value<0.01/300 & MeanDiffer>0) %>% ungroup() %>%  
              filter(trait == 'GE', year == '2020_2021')%>%
              mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2','AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
                     label = mapvalues(Effects, from = c('Base_fixef','Full','JC','JC_Fix','NoFix'),
                                       to = c('*    ',' *   ','  *  ','   * ','    *')),
                     Variation = Effects),
            aes(y = 1, label =label, group = Variation))+
  labs(y = 'Prediction Accuracy')
dev.off()

# Supp Figure 7 GI Structured GBLUP ######
png('picsPNGforQsd1Effects_paper/StrucutredGBLUPGI_Results_sup.png', 1400, 800, res =120)
FilesStructuredGBLUP %>% filter(trait == 'GI') %>%  filter(year == '2020_2021') %>%
  mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
         Variation = Effects) %>%
  ggplot(aes(x = Model, y = correlation, fill = Variation, color = Variation))+
  geom_boxplot()+facet_nested(Qsd1~trait+TP, space = 'free_x')+
  geom_hline(yintercept = 0)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  geom_text(data = Significant_differents %>% filter(P.value<0.01/300 & MeanDiffer>0) %>% ungroup() %>% 
              filter(trait == 'GI', year == '2020_2021')%>%
              mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2','AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
                     label = mapvalues(Effects, from = c('Base_fixef','Full','JC','JC_Fix','NoFix'),
                                       to = c('*    ',' *   ','  *  ','   * ','    *')),
                     Variation = Effects),
            aes(y = 1, label =label, group = Variation))+
  labs(y = 'Prediction Accuracy')
dev.off()


# Supp Figure 8 Segregation distortion tests #######
taxa.Family = AllDHBluesPerYear %>% select(Family, taxa) %>% unique()
FamWinterGD = taxa.Family %>% join(WinterGD)  

Seg_distortionTest = data.frame()
for (i in c('BS7','BS8','BS9')){
  tempWinterGD = WinterGD %>% filter(substr(taxa,1,3) == i)
  Qsd1_0 = tempWinterGD %>% filter(Qsd1 == 0)
  Qsd1_2 = tempWinterGD %>% filter(Qsd1 == 2)
  for (colnam in colnames(tempWinterGD[,-1])){
    Ttest.results = tryCatch(expr = t.test(Qsd1_0[[colnam]],Qsd1_2[[colnam]]),
                             error = function(e){data.frame(p.value = NA)})
    
    Seg_distortionTest = rbind(Seg_distortionTest, 
                               data.frame(Family = i, SNP = colnam, p.value = Ttest.results$p.value))
    
  }
}

ChrBreaks = WinterGM %>% dplyr::mutate(ordinal = 1:n()) %>% group_by(Chromosome) %>% slice_head(n = 1) %>% select(Chromosome, ordinal)

WinterGM %>% dplyr::mutate(ordinal = 1:n()) %>% filter(SNP == 'Qsd1')

png(filename = 'picsPNGforQsd1Effects_paper/SegDistortionTest.png', width = 1000, height = 600, res = 120)
Seg_distortionTest %>% join(WinterGM %>% dplyr::mutate(ordinal = 1:n())) %>% 
  mutate(Family = mapvalues(Family, from = c('BS7','BS8','BS9'),
                            to = c('Scala/Lightning','SY_Tepee/Lightning','Wintmalt/Lightning')),
         logPval = -log10(p.value)) %>%
  arrange(Chromosome, Position) %>%
  ggplot(aes(x = ordinal, y = logPval))+geom_point()+facet_grid(rows =vars(Family))+xlab('Chromosome')+
  geom_vline(xintercept = ChrBreaks$ordinal)+
  scale_x_continuous(label = c("1H","2H", "3H", "4H", "5H", "6H", "7H", "UN"),
                     breaks = ChrBreaks$ordinal+300)+
  geom_vline(xintercept = 5148, color = 'red')+geom_text(x = 5148, y = 50, label = 'Qsd1', color = 'red')
dev.off()



# Supp Figure 9 Markers responding to Qsd1 #####
png(filename = 'picsPNGforQsd1Effects_paper/Loci1toLoci2.png', width = 700, height = 600, res = 120)
data.frame(x = 1:100) %>% mutate(Dormant_density = pnorm(mean = 20, sd = 5, q = x), 
                                 NonDormant_density = pnorm(mean = 10, sd = 5, q = x),
                                 Dormant_cdf = pnorm(mean = 70, sd = 5, q = x), 
                                 NonDormant_cdf = pnorm(mean = 50, sd = 5, q = x)) %>%
  pivot_longer(cols = !x) %>% separate(name, into = c('name', 'facet'))%>% 
  mutate(facet = ifelse(facet=='density',"Non-dormant Qsd1 \nEffect of Locus2 on dormancy",
                        "Dormant Qsd1 \nEffect of Locus2 on dormancy"))%>%
  ggplot(aes(x =x, y= value, group = name, color = name))+geom_line(size = 2) +
  facet_grid(rows = vars(facet), scales ='free', switch = 'y')+
  labs(color = "Status of second locus",x = 'Time post maturity', y =NULL)+theme_bw()+
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), )
dev.off()


# Supp Figure 10 GE and GI TP1 model training and prediction by Qsd1 grouping.  ######
NOFlaviaDH2021_TP1GE_GI_PredictOutgroup = read_csv('Analysis/NOFlaviaDH2021_TP1GE_GI_PredictOutgroup.csv')
NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest = read_csv('Analysis/NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest.csv')

png(filename = 'picsPNGforQsd1Effects_paper/Sup_GEGI2021TP1Predictions.png', width = 1000, height = 700, res = 120)
HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(year == '2021', TP == 'TP1')  %>%
  filter(substr(type, 1,2)=='h2', type != 'h2_Narrow_Qsd1Accounted') %>%
  mutate(maxCor = sqrt(value),
         type =   mapvalues(x = type, from = c('h2_Narrow','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D' ), 
                            to = c('sigma[a]','sigma[aN]','sigma[aD]'))) %>% 
  ggplot(aes(type, maxCor, fill = type))+
  geom_bar(stat = 'identity',position = 'dodge')+facet_grid(TP+trait~'h^2',margins = F)+
  ylim(0,1)+ xlab('Term')+scale_x_discrete(labels = names)+
  theme_bw()+ theme(legend.position = 'none', strip.background.y = element_blank(), strip.text.y = element_blank())+
  ylab('Max expected prediction accuracy')+
  
NOFlaviaDH2021_TP1GE_GI_PredictOutgroup %>% 
  filter(Qsd1 != "AllAlleles", n_test >8) %>%
  group_by(trait,Family,Qsd1Train, TP, trainingProportion, Qsd1) %>% 
  dplyr::summarise(stdev = sd(correlation, na.rm = T),
                   correlation = mean(correlation,na.rm = T), 
                   n_train = mean(n_train), 
                   n_test = mean(n_test)) %>% mutate(type = 'Training_Qsd1') %>%
  dplyr::rename(Qsd1Pred = Qsd1) %>% 
  mutate(Qsd1Train = as.character(Qsd1Train),
         Qsd1Pred = as.character(Qsd1Pred)) %>%
  rbind(NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest %>% mutate(Qsd1Train = 'AllAlleles') %>%
          dplyr::rename(Qsd1Pred = Qsd1) %>%
          group_by(trait, Family, Qsd1Train, TP, trainingProportion, Qsd1Pred) %>%
          filter(n_test>5) %>%
          dplyr::summarize(stdev = sd(correlation, na.rm = T),
                           correlation = mean(correlation,na.rm = T), 
                           n_train = mean(n_train), 
                           n_test = mean(n_test)) %>% 
          mutate(type = 'Training_Overall'))  %>%
  filter(Family == 'Overall') %>%
  mutate(Qsd1Train = mapvalues(x = Qsd1Train, from = c('0','2','AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles')),
         Qsd1Pred = mapvalues(x = Qsd1Pred, from = c('0','2', 'AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles')),
         Trained_Predicted = paste0(Qsd1Train, ' \U2192 ',Qsd1Pred)) %>%
  ggplot(aes(x = n_train, y = correlation, color = Trained_Predicted, linetype = Trained_Predicted))+
  # geom_point(aes(shape = Trained_Predicted))+
  facet_grid(trait+TP~'Predictions', scales = 'free') + theme_bw()+
  xlab('Number of training taxa') +
  labs(color = 'Qsd1 Trained \U2192 Predicted', linetype = 'Qsd1 Trained \U2192 Predicted',shape = 'Qsd1 Trained \U2192 Predicted') +
  geom_smooth(method = 'loess', se =T, size = 1.2)+labs(y = 'Realized prediction accuracy')+
  scale_color_manual(values = c('salmon','salmon','salmon','mediumseagreen','mediumseagreen','cornflowerblue','cornflowerblue'))+
  # scale_linetype_manual(values = c('3131','dotted','solid','dotted','solid','dotted','solid'))+
  scale_linetype_manual(values = c('solid','dotted','3131','solid','dotted','dotted','solid'))+
  ylim(0,1)+  plot_layout(ncol = 2, design = c('ABBBBB'))
dev.off()

# Supp figure ? Extra 2021 Structured GBLUP ######

png('picsPNGforQsd1Effects_paper/Supp_2021StrucutredGBLUP_Results.png', 1400, 800, res =120)
FilesStructuredGBLUP %>%  filter(year == '2021_') %>%
  mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant', 'Both Alleles')),
         Variation = Effects) %>%
  ggplot(aes(x = Model, y = correlation, fill = Variation, color = Variation))+
  geom_boxplot()+facet_nested(Qsd1~trait+TP, space = 'free_x')+
  geom_hline(yintercept = 0)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))+
  geom_text(data = Significant_differents %>% filter(P.value<0.01/300 & MeanDiffer>0) %>% ungroup() %>%  filter(year == '2021_')%>%
              mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles')),
                     label = mapvalues(Effects, from = c('Base_fixef','Full','JC','JC_Fix','NoFix'),
                                       to = c('*    ',' *   ','  *  ','   * ','    *')),
                     Variation = Effects),
            aes(y = 1, label =label, group = Variation))+
  ylab('Prediction Accuracy')
dev.off()


# Supp Figure ? Correlation of marker effects by family.  ######
png('picsPNGforQsd1Effects_paper/Sup_CorOfFamGroupMarkerEffects.png', 800, 600, res =120)
WinterMarkerEffectbyFamily %>% filter(trait == 'GI', Family != 'DH130910')%>% ungroup() %>%
  mutate(Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  pivot_wider(names_from = c(Family,TP), values_from = Effect) %>%
  ungroup() %>% select(!c(trait,Marker)) %>% cor(use = 'complete.obs') %>% data.frame() %>% 
  rownames_to_column(var='TP_Family') %>%
  pivot_longer(cols = !TP_Family) %>% ggplot(aes(x = TP_Family, y = name, fill = value) )+geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5) +
  geom_text(aes(label = round(value,2)), size = 2) +
  xlab('Family_Timepoint')+ylab('Family_Timepoint')
dev.off()
