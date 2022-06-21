library(lme4);library(Hmisc);library(plyr);library(readxl);
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x)
library(tidyverse)
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
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/BluesByFamilyQsd12020_2021.png', 1400, 800, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Non-dormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020/2021'))  %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  mutate(Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family),
         trait = ifelse(trait == 'GI','Germiantion Rate','Germination Percent')) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  theme_bw()+
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free')+
  labs(y = 'trait value',x = 'Days of after-ripening (PM+days)', fill = 'Qsd1 Status')+
  scale_x_discrete(label = c(5,19,47,96,152)) 

dev.off()

# Figure 1 JPEG ######
jpeg('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/figure1.jpg', 1400, 800, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Non-dormant')) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% filter(year %in% c('2020/2021'))  %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(!(trait =='GE' &value>1.05)) %>%
  mutate(Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family),
         trait = ifelse(trait == 'GI','Germiantion Rate','Germination Percent')) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  theme_bw()+
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free')+
  labs(y = 'trait value',x = 'Days of after-ripening (PM+days)', fill = 'Qsd1 Status')+
  scale_x_discrete(label = c(5,19,47,96,152)) 
dev.off()

Qsd1Key = WinterGD[,c('taxa','Qsd1')]

AllDHBluesPerYear %>% filter(TP == 'TP1',trait == 'GE',year == '2020/2021') %>%
  join(Qsd1Key) %>% select(Family, Qsd1) %>% table()



FMV_Qsd1 = FamilyMeansAndVarianceQsd1 %>% ungroup() %>% filter(type =='Fixed') %>% group_by(year, TP,trait, model) %>% 
  group_modify(~{.x%>% mutate(value = .x$value + .x$value[length(.x$value)])}) %>%ungroup() %>%
  rbind(FamilyMeansAndVarianceQsd1 %>% filter(type != 'Fixed')) %>% arrange(year, TP, trait, model)

# Figure 2 Variance partitioning and heritability  #######
HeritabilityPlot =
  HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(year == '2020/2021',trait == 'GI') %>%
  mutate(peralle = ifelse(substr(type,1,2)=='H2', "overall",'Per allele')) %>%
  ggplot(aes(x = TP, y = value, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge')+theme_bw() +
  scale_fill_manual( values = c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'), labels = H2Names)+
  theme(legend.text.align = 0)+
  labs(fill = 'Heritability type', y = 'Heritability',x = NULL)+
  scale_x_discrete(label = c(expression(PM[5]),expression(PM[19]),expression(PM[47]),expression(PM[96]), expression(PM[152])))

GIVarianceAttribution = FMV_Qsd1 %>% filter(model %in% c('base','difference')) %>% filter(trait =='GI',year =='2020/2021') %>% 
  filter(substr(Coef, 1,2) %in% c('Fa', 'at')) %>% 
  mutate(Coef = gsub(pattern = "at\\(Family, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef)),
         Coef = gsub('Family_','',Coef)) %>%
  mutate(type = ifelse(type=='Fixed',type,'Variance')) %>%
  filter(type == 'Variance') %>%filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5')))%>%
  mutate(model2 = ifelse(model == 'base','Family Genetic Variance\nafter Qsd1 is accounted for',
                         'Variance attributable\nto Qsd1'),
         Coef = gsub(pattern = 'DH130910',replacement = 'Lightning',Coef)) %>%
  ggplot(aes(x = TP, y = value, group = Coef, fill = Coef, alpha = model2))+
  geom_col(position = 'dodge', color = 'black')+scale_alpha_manual(values = c(.5,1), guide = 'none')+
  labs(fill = 'Family', alpha = 'Variance attribution',y = 'Genetic Variance', x= NULL)+theme_bw(base_size = 13) +
  scale_x_discrete(label = c(expression(PM[5]), expression(PM[19]),expression(PM[47]), expression(PM[96]),  expression(PM[152])))

Qsd1Variance = atFamQsd1TaxaVar %>% filter(substr(Coef,1,2)=='at' ) %>%
  mutate(Coef = gsub(pattern = "at\\(FamilyQsd1, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef))) %>%
  separate(col = Coef, into = c('Family','Qsd1'), sep = '130910',remove = F) %>% 
  mutate(Family = paste0(Family, '130910'), Qsd1lab = ifelse(Qsd1==2, 'Dormant','Non-dormant'),
         Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  filter(trait == 'GI',year == '2020/2021') %>%
  ggplot(aes(TP, value,fill= Qsd1lab))+
  geom_bar(position = 'dodge',stat= 'identity', color = 'black') +
  facet_wrap(~Family, ncol = 1) +
  scale_alpha_manual(values = c(.3,1)) + theme_bw()+
  labs(fill = 'Qsd1 status', y = 'Genetic variance', x = NULL)+
  scale_x_discrete(label = c(expression(PM[5]),expression(PM[19]),expression(PM[47]),expression(PM[96]),  expression(PM[152]))) 

png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/ComboFigure.png', 1200, 800, res =120)
HeritabilityPlot+#
  GIVarianceAttribution+
  Qsd1Variance+
  plot_layout(design = 'AAC\nBBC')+plot_annotation(tag_levels = 'A')
dev.off()

jpeg('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/figure2.jpg', 1200, 800, res =120)
HeritabilityPlot+#
  GIVarianceAttribution+
  Qsd1Variance+
  plot_layout(design = 'AAC\nBBC')+plot_annotation(tag_levels = 'A')
dev.off()

# Figure 3 Correlation of Marker effects between Qsd1 groups #####
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

png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/CorOfQsd1GroupMarkerEffects.png', 800, 500, res =120)
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
            PM = mapvalues(TP, from = c('TP1','TP1.5','TP2','TP2.5','TP3','TP3.5','TP4','TP5'),
                           to = c(5,12,19,33,47,68,96,152)),
            QSD1_PM = paste0(Qsd1, '_' ,PM),
            name2 = paste0(Qsd11,'_',mapvalues(TPP, from = c('TP1','TP1.5','TP2','TP2.5','TP3','TP3.5','TP4','TP5'),
                                               to = c(5,12,19,33,47,68,96,152))),
            value, tval, label = paste0(round(value,2), ifelse(tval<0.05,'Sig','')), 
            color = ifelse(tval<0.05,'colors','no')) %>%
  mutate(label = gsub('NA','', label),
         QSD1_PM = factor(QSD1_PM, levels = c('D_5','D_12','D_19','D_33','D_47','D_68','D_96','D_152',
                                              'ND_5','ND_12','ND_19','ND_33','ND_47','ND_68','ND_96','ND_152')),
         name2 = factor(name2, levels = c('D_5','D_12','D_19','D_33','D_47','D_68','D_96','D_152',
                                          'ND_5','ND_12','ND_19','ND_33','ND_47','ND_68','ND_96','ND_152')),
         label = ifelse(label == '1Sig','1',label))%>% 
  ggplot(aes(x = QSD1_PM, y = name2, fill = value) )+geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5) +
  geom_text(aes(label = label), size = 2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + labs(fill = 'Correlation')+
  xlab('Qsd1 Allele_Timepoint')+ylab('Qsd1 Allele_Timepoint') +
  geom_vline(xintercept = 8.5)+geom_hline(yintercept = 8.5)
dev.off()

jpeg('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/figure3.jpg', 800, 500, res =120)
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
            PM = mapvalues(TP, from = c('TP1','TP1.5','TP2','TP2.5','TP3','TP3.5','TP4','TP5'),
                           to = c(5,12,19,33,47,68,96,152)),
            QSD1_PM = paste0(Qsd1, '_' ,PM),
            name2 = paste0(Qsd11,'_',mapvalues(TPP, from = c('TP1','TP1.5','TP2','TP2.5','TP3','TP3.5','TP4','TP5'),
                                               to = c(5,12,19,33,47,68,96,152))),
            value, tval, label = paste0(round(value,2), ifelse(tval<0.05,'Sig','')), 
            color = ifelse(tval<0.05,'colors','no')) %>%
  mutate(label = gsub('NA','', label),
         QSD1_PM = factor(QSD1_PM, levels = c('D_5','D_12','D_19','D_33','D_47','D_68','D_96','D_152',
                                              'ND_5','ND_12','ND_19','ND_33','ND_47','ND_68','ND_96','ND_152')),
         name2 = factor(name2, levels = c('D_5','D_12','D_19','D_33','D_47','D_68','D_96','D_152',
                                          'ND_5','ND_12','ND_19','ND_33','ND_47','ND_68','ND_96','ND_152')),
         label = ifelse(label == '1Sig','1',label))%>% 
  ggplot(aes(x = QSD1_PM, y = name2, fill = value) )+geom_tile()+
  scale_fill_gradient2(limits=c(0,1),low = 'blue',high = 'red', midpoint =.5) +
  geom_text(aes(label = label), size = 2)+
  theme_bw()+theme(axis.text.x = element_text(angle = 90)) + labs(fill = 'Correlation')+
  xlab('Qsd1 Allele_Timepoint')+ylab('Qsd1 Allele_Timepoint') +
  geom_vline(xintercept = 8.5)+geom_hline(yintercept = 8.5)
dev.off()


# Figure 4 cross allele predictions #######
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Bettern_100CrossAllele2.png', 800, 400, res =120)
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
  mutate(Qsd1Train = mapvalues(x = Qsd1Train, from = c('0','2','AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles'),warn_missing = F),
         Qsd1Pred = mapvalues(x = Qsd1Pred, from = c('0','2', 'AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles'),warn_missing = F),
         Trained_Predicted = paste0(Qsd1Train, ' \U2192 ',Qsd1Pred)) %>%
  filter(trainingProportion == 100) %>%
  # mutate(TP = mapvalues(TP, from = c('TP1','TP2','TP3','TP4','TP5'),
  #                       to = c(5,19,47,96, 152), warn_missing = F),
  #        TP = factor(TP, levels = c(5,19,47,96, 152))) %>%
  ggplot(aes( x = TP, y = correlation, color = Trained_Predicted, shape = Trained_Predicted,group = Trained_Predicted))+
  geom_point(size = 4, stroke = 5)+theme_bw()+geom_line()+
  scale_shape_manual(values = c('1','2','3','4','5','6','7'))+
  scale_color_manual(values = c('salmon','salmon','salmon','mediumseagreen','mediumseagreen','cornflowerblue','cornflowerblue'))+
  labs(y = 'Prediction accuracy', x = 'Days of after-ripening')+ ylim(0,1)+geom_hline(yintercept = 0)+
  scale_x_discrete(label = c(5,19,47,96,152))
dev.off()  

# Graphical abstract structured gblup results #####
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/StructGIforPrettyPlot.png', 800, 500, res =120)
FilesStructuredGBLUP %>% filter(year == '2020_2021') %>% 
  filter(Qsd1 == 'AllAlleles') %>%
  mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles'),warn_missing =F ),
         Variation = Effects) %>%filter(Qsd1=='Both Alleles') %>%
  filter(Model == 'Base', trait == 'GI')%>% group_by(Family, Qsd1, Model, Effects,TP, trait, year) %>%
  summarize(Correl = mean(correlation), sd = sd(correlation)) %>% ungroup() %>%
  
  ggplot(aes(x = TP, y = Model))+geom_text(aes(label =round(Correl, 3)))+theme_bw()+
  labs(x = 'Days of after-ripening' ,y = NULL)+  geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+
  scale_y_discrete(label = expression(GBLUP: ~y %==% mu + bold(Z[0]) * bold(mu[0])))+
  # annotate('text', label='Prediction ability within GBLUP', x=-Inf, y=Inf, hjust=-0.1, vjust=0.90, size =2)+
  # scale_x_discrete(label = c(5,19,47,96, 152), position = 'top')+
  scale_x_discrete(position = 'top', label = c(5,19,47,96,152))+
  theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 12))+
  
  FilesStructuredGBLUP %>% filter(year == '2020_2021') %>% 
  filter(Qsd1 == 'AllAlleles') %>%
  mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles'),warn_missing =F ),
         Variation = Effects) %>%filter(Qsd1=='Both Alleles') %>%
  filter(Model == 'Qsd1', trait == 'GI')%>% group_by(Family, Qsd1, Model, Effects,TP, trait, year) %>%
  summarize(Correl = mean(correlation), sd = sd(correlation)) %>% join(Significant_differents) %>%
  mutate(Correl = ifelse(P.value<0.001, Correl, NA), #PercentDif = ifelse(P.value < 0.001, PercentDif, NA),
         Variation = Effects) %>%
  filter(Variation != 'JC') %>%
  mutate(Variation = factor(Variation, levels = c('Full','NoFix','JC_Fix','Base_fixef')))%>%
  ggplot(aes(x = TP,y = Variation, color = PercentDif, size = PercentDif))+geom_point() +theme_bw()+
  labs(color = 'Percent\nincrease\nfrom GBLUP', size = 'Percent\nincrease\nfrom GBLUP',
       #x = 'Days post maturity',
       x = NULL,
       y = NULL)+
  scale_color_continuous(limits = c(-3, 14),breaks = c(0,4.,8.,12.), low = 'white', high = 'green')+
  guides(color=guide_legend(), size = guide_legend())+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.text.y = element_text(size = 14))+
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5))+
  scale_y_discrete(label =
                     c(expression(y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
                       expression(y%==% mu+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
                       expression(y%==% mu+bold(X)*bold(beta)+bold(Z[1])*bold(mu[1])),
                       expression(y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0]))
                     )
  )+
  
  plot_layout(nrow = 2, design = c('A\nB\nB\nB\nB\nB'))
dev.off()

# Figure 5 GI Structure gblups ######
BetterPLotsSTGBLUP= FilesStructuredGBLUP %>%
  group_by(year,trait, TP, Qsd1) %>%
  group_modify(~{
    Base = .x %>% filter(Model == 'Base')
    .x %>% group_by(Model, Effects) %>%
      group_modify(~{
        ttest = t.test(x = .$correlation, y = Base$correlation)
        data.frame(P.value = ttest$p.value, MeanDiffer = mean(.$correlation)-mean(Base$correlation), 
                   PercentDif = 100*(mean(.$correlation)-mean(Base$correlation))/mean(Base$correlation)) 
        
      }) %>% mutate(label = ifelse(Model=='Base',as.character(round(mean(Base$correlation),2)),NA),
                    label = ifelse(PercentDif< -15, 'D',label),
                    PercentDif = ifelse(PercentDif< -15, -15, PercentDif),
                    PercentDif = ifelse(PercentDif>65, 65, PercentDif))
  }) %>% mutate(Qsd1 = mapvalues(Qsd1, from = c('0','2', 'AllAlleles'), to = c('Non-dormant','Dormant','Both Alleles'),warn_missing =F ))

png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/StructGI.png', 1200, 700, res =120)
BetterPLotsSTGBLUP %>% filter(year  == '2020_2021') %>% ungroup() %>%
  filter(Qsd1 == 'Both Alleles',trait == 'GI', Model != 'Base') %>%
  rbind(BetterPLotsSTGBLUP %>% filter(year == '2020_2021') %>% 
          filter(Qsd1 == 'Both Alleles',trait == 'GI', Model == 'Base')%>% ungroup() %>%
          rbind(.,.,.,.) %>% mutate(Model = rep(c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),each =5))) %>%
  mutate(Variation = factor(Effects, levels = c('Base','Base_fixef','JC','JC_Fix','NoFix','Full')),
         Model = mapvalues(Model, from = c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),
                           to = c('Family','Qsd1','Qsd1:Family interaction','Qsd1 + Family'))) %>% 
  filter(Variation != 'JC') %>%
  ggplot(aes(x = TP, y = Variation, size = abs(PercentDif), color = PercentDif))+geom_point()+
  geom_text(aes(label =label), color = 'black', size = 3.88)+theme_bw()+
  labs(x = 'Days of after-ripening', y = NULL)+ facet_wrap(~Model)+
  scale_y_discrete(label = c(
    expression(Full:~ y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(NoFix:~y%==% mu+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(JC_Fix:~ y%==% mu+bold(X)*bold(beta)+bold(Z[1])*bold(mu[1])),
    # expression(JC:~y%==% mu+bold(Z[1])*bold(mu[1])),
    expression(Base_fixef:~y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])),
    expression(Base~GBLUP: ~y %==% mu + bold(Z[0]) * bold(mu[0]))
  ), limits = rev)+
  scale_colour_gradient2(name = 'Percent Difference\nfrom GBLUP',low = 'brown2', mid = 'white',high = 'green')+
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))+
  scale_size_continuous(name = NULL, range = c(0,15), guide='none')+
  geom_vline(xintercept =  c(1.5,2.5,3.5,4.5), alpha = 0.3)+
  scale_x_discrete(label = c(5,19,47,96,152))

dev.off()

jpeg('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/figure5.jpg', 1200, 700, res =120)
BetterPLotsSTGBLUP %>% filter(year  == '2020_2021') %>% ungroup() %>%
  filter(Qsd1 == 'Both Alleles',trait == 'GI', Model != 'Base') %>%
  rbind(BetterPLotsSTGBLUP %>% filter(year == '2020_2021') %>% 
          filter(Qsd1 == 'Both Alleles',trait == 'GI', Model == 'Base')%>% ungroup() %>%
          rbind(.,.,.,.) %>% mutate(Model = rep(c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),each =5))) %>%
  mutate(Variation = factor(Effects, levels = c('Base','Base_fixef','JC','JC_Fix','NoFix','Full')),
         Model = mapvalues(Model, from = c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),
                           to = c('Family','Qsd1','Qsd1:Family interaction','Qsd1 + Family'))) %>% 
  filter(Variation != 'JC') %>%
  ggplot(aes(x = TP, y = Variation, size = abs(PercentDif), color = PercentDif))+geom_point()+
  geom_text(aes(label =label), color = 'black', size = 3.88)+theme_bw()+
  labs(x = 'Days of after-ripening', y = NULL)+ facet_wrap(~Model)+
  scale_y_discrete(label = c(
    expression(Full:~ y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(NoFix:~y%==% mu+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(JC_Fix:~ y%==% mu+bold(X)*bold(beta)+bold(Z[1])*bold(mu[1])),
    # expression(JC:~y%==% mu+bold(Z[1])*bold(mu[1])),
    expression(Base_fixef:~y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])),
    expression(Base~GBLUP: ~y %==% mu + bold(Z[0]) * bold(mu[0]))
  ), limits = rev)+
  scale_colour_gradient2(name = 'Percent Difference\nfrom GBLUP',low = 'brown2', mid = 'white',high = 'green')+
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))+
  scale_size_continuous(name = NULL, range = c(0,15), guide='none')+
  geom_vline(xintercept =  c(1.5,2.5,3.5,4.5), alpha = 0.3)+
  scale_x_discrete(label = c(5,19,47,96,152))

dev.off()


# Supplemental Figure 1 Various model and the PA of each tests #####
filesModel = dir("Analysis/BGLROutput/", full.names = T, pattern = "Various")
VariousModels = lapply(filesModel,FUN = read.csv) %>% bind_rows() %>% select(!X)
VariousModels %>%head()
library(stringr)
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Sup_VariousModelsPA.png', 1400, 1200, res =120)
VariousModels %>% filter(Family == 'Overall') %>% mutate(fixed = ifelse(str_sub(Model,-5,-1)=='Fixed','Qsd1_Fixed','No_fixed'),
                                                         Model = gsub('_Qsd1_Fixed','',Model)) %>%
  group_by(trait, n_train, Qsd1Status, fixed) %>%
  group_modify(~{
    base = .x %>% filter(Model == 'BRR'|Model == 'BRR_Qsd1_Fixed')
    meanbase = mean(base$correlation)
    .x %>% group_by(Model) %>% group_modify(~{
      .x %>% summarise(percentDif = (mean(correlation)-meanbase)/meanbase*100)
    })%>% mutate(label = ifelse(Model == 'BRR', as.character(round(meanbase,2)),''))
  }) %>% 
  mutate(Qsd1Status = mapvalues(Qsd1Status, from = c('-1','1','AllAlleles'), to = c('NonDormant','Dormant','Both Alleles'), warn_missing = F),
         Model = factor(Model, levels = c('BRR','BayesB','BayesC','LASSO'))) %>%
  ggplot(aes(x = Model, Qsd1Status, color = percentDif, size = abs(percentDif)))+
  geom_point()+ labs(y = 'Prediction Accuracy')+ geom_text(aes(Model, Qsd1Status, label = label), color = 'black', size = 3) +
  scale_color_gradient2(low = 'blue',high = 'red',limits = c(-5,5))+ labs(color = 'Percent Difference')+
  scale_size_continuous(name = NULL, guide='none') + geom_hline(yintercept = c(1.5,2.5))+
  facet_grid(n_train~trait+fixed, scales ='free_x')

dev.off()

VariousModels %>% filter(Family == 'Overall') %>%
  group_by(trait, n_train, Qsd1Status) %>%
  group_modify(~{
    base = .x %>% filter(Model == 'BRR')
    meanbase = mean(base$correlation)
    .x %>% group_by(Model) %>% group_modify(~{
      .x %>% summarise(percentDif = (mean(correlation)-meanbase)/meanbase, label = '')
    }) %>% mutate(label = ifelse(Model == 'BRR'))
  }) %>%
  mutate(Qsd1Status = mapvalues(Qsd1Status, from = c('-1','1','AllAlleles'), to = c('NonDormant','Dormant','Both Alleles'), warn_missing = F)) %>%
  ggplot(aes(x = Model, Qsd1Status, color = percentDif, size = abs(percentDif)))+
  geom_point()+ labs(y = 'Prediction Accuracy')+
  scale_color_gradient2()+
  facet_grid(n_train~trait, scales ='free_x')




VariousModels %>% filter(Family == 'Overall', n_train == 70) %>%
  ggplot(aes(x = Model, y = correlation, fill = Model, color = Qsd1Status))+
  geom_boxplot()+ labs(y = 'Prediction Accuracy')+
  facet_grid(trait+n_train~Family, scales ='free_x')


# Supplemental Figure 2 All the blues over time from each datatset #####
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Sup_BluesByFamilyQsd1AllYears.png', 2400, 1500, res =120)
AllDHBluesPerYear %>%  filter(type == 'BLUE') %>% mutate(year = factor(year, levels = c('2020','2021','2020/2021')))%>%
  join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1!= 1) %>% mutate(Qsd1= as.factor(Qsd1)) %>%
  filter(Family %nin% c('Cha','End','DH130910')) %>% 
  mutate(Qsd1= ifelse(Qsd1==2,'Dormant','Non-dormant'),
         Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  filter(!(trait =='GE' &value>1.05)) %>% mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  ggplot(aes(x = TP, y = value, fill = Qsd1))+  theme_bw()+
  geom_boxplot()+facet_nested(trait~year+Family, scales = 'free_y', space = 'free_x')+
  labs(x = 'Timepoint',fill = ' Qsd1 status')+
  scale_x_discrete(label = c(5,12,19,33,47,68,96,152))
dev.off()

# Supplemental figure 3 Heritabilities of the traits over time  #####
HeritabilitiesWinterGermTraits$type %>% unique()

H2Names = c(expression(Cullis~H^2),
            expression(Cullis~H^2~Qsd1~Fixed),
            expression(h^2),
            expression(h^2 ~ Qsd1~Fixed),
            expression(h^2 ~ Qsd1[D]),
            expression(h^2 ~ Qsd1[N]))

# "H2_Cullis", "H2_Cullis_Qsd1Accounted", "h2_Narrow","h2_Narrow_Qsd1Accounted", "h2_NarrowU_Qsd1D","h2_NarrowU_Qsd1N",   

png(filename = 'WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Heritabilites.png',
    width = 1000, height = 800, res = 120)
HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%  ggplot(aes(x = TP, y = value, fill = type)) +
  geom_bar(stat = 'identity', position = 'dodge')+theme_bw() +
  scale_fill_manual( values = c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd'), labels = H2Names)+
  theme(legend.text.align = 0  )+
  facet_nested('Year'+year~trait) + labs(fill = 'Heritability type', y = 'Heritability',x = 'Timepoint')+
  scale_x_discrete(label = c(expression(PM[5]),  expression(PM[12]),expression(PM[19]),expression(PM[33]),expression(PM[47]),
                             expression(PM[68]),expression(PM[96]),  expression(PM[152])))

dev.off()


HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>% 
  pivot_wider(names_from = 'type',values_from='value') %>% arrange(desc(trait)) %>%
  mutate(ratio = h2_NarrowU_Qsd1N/h2_NarrowU_Qsd1D) %>% view()



# Supplemental figure 4 Variance attributable to Qsd1  ######
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/VarAttributableToQsd1.png', 1400, 800, res =120)
FMV_Qsd1 %>% filter(model %in% c('base','difference')) %>% filter(trait =='GE') %>%
  filter(substr(Coef, 1,2) %in% c('Fa', 'at')) %>% 
  mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
  mutate(Coef = gsub(pattern = "at\\(Family, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef)),
         Coef = gsub('Family_','',Coef)) %>%
  mutate(type = ifelse(type=='Fixed',type,'Variance')) %>%
  filter(type == 'Variance') %>%filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  mutate(model2 = ifelse(model == 'base','Family Genetic Variance\nafter Qsd1 is accounted for',
                         'Variance attributable\nto Qsd1'),
         Coef = gsub(pattern = 'DH130910',replacement = 'Lightning',Coef)) %>%
  ggplot(aes(x = TP, y = value, group = Coef, fill = Coef, alpha = model2))+
  geom_col(position = 'dodge', color = 'black')+facet_nested(year~trait, scales= 'free')+scale_alpha_manual(values = c(.5,1))+
  labs(fill = 'Family', alpha = 'Variance attribution', y = 'Genetic Variance', x = 'Timepoint') + theme_bw(base_size = 13) +
  scale_x_discrete(label = c(expression(PM[5]),  expression(PM[12]),expression(PM[19]),expression(PM[33]),expression(PM[47]),
                             expression(PM[68]),expression(PM[96]),  expression(PM[152])))+
  
  FMV_Qsd1 %>% filter(model %in% c('base','difference')) %>% filter(trait =='GI') %>%
  mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
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
  labs(fill = 'Family', alpha = 'Variance attribution',y = 'Genetic Variance', x= 'Timepoint')+theme_bw(base_size = 13) +
  scale_x_discrete(label = c(expression(PM[5]),  expression(PM[12]),expression(PM[19]),expression(PM[33]),expression(PM[47]),
                             expression(PM[68]),expression(PM[96]),  expression(PM[152]))) +
  plot_layout(ncol = 2, guides = 'collect')
dev.off()

# Supplemental Figure 5 Heritability of the remining variance within families after Qsd1 is accounted for #####
# What is the ratio of error variance to the at(Family):taxa geentic variance? ie a heritability after Qsd1 is accounted for?
png(filename = 'WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/RemainingFamilyGenVarH2.png',
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
  mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
  ggplot(aes(x = TP, y = H2,group = Coef, fill = Coef))+theme_bw()+
  geom_col(position = 'dodge')+facet_nested('Year'+ year~'Trait'+trait, scales= 'free')+
  labs(group = 'Family',fill = 'Family', y = 'Broad Sense H2', x = 'Timepoint')+
  scale_x_discrete(label = c(expression(PM[5]),  expression(PM[12]),expression(PM[19]),expression(PM[33]),expression(PM[47]),
                             expression(PM[68]),expression(PM[96]),  expression(PM[152])))

dev.off()

# Supplemental Figure 6 Different variances per Family:Qsd1 grouping######
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/DifferentVarForFamilyQsd1.png', 1400, 800, res =120)
atFamQsd1TaxaVar %>% filter(substr(Coef,1,2)=='at' ) %>%
  mutate(Coef = gsub(pattern = "at\\(FamilyQsd1, ", replacement ="",x = gsub(pattern = '\\):taxa',replacement = '', x = Coef))) %>%
  separate(col = Coef, into = c('Family','Qsd1'), sep = '130910',remove = F) %>% 
  mutate(Family = paste0(Family, '130910'), Qsd1lab = ifelse(Qsd1==2, 'Dormant','Non-dormant'),
         Family = gsub(pattern = 'DH130910',replacement = 'Lightning',Family)) %>%
  filter(!(year == '2020/2021' & TP %in% c('TP1.5','TP2.5','TP3.5'))) %>%
  mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
  ggplot(aes(TP, value,fill= Qsd1lab))+
  geom_bar(position = 'dodge',stat= 'identity')+
  facet_nested(trait+year~Family, scales= 'free_y') +scale_alpha_manual(values = c(.5,1)) +theme_bw()+
  geom_text(inherit.aes = F, 
            data = atFamQsd1TaxaVar %>% filter(Coef == 'lrt'|substr(Coef,1,2)=='at') %>% 
              mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
              mutate(siglab = ifelse(value<0.01,'','ns'), ypos = max(value[1:7])*1.1) %>% filter(Coef == 'lrt',year != '2020/2021'),
            aes(x = TP, y = ypos, label = siglab))+
  labs(alpha = 'Qsd1 status', y = 'Genetic variance', x = 'Timepoint')+
  scale_x_discrete(label = c(expression(PM[5]),  expression(PM[12]),expression(PM[19]),expression(PM[33]),expression(PM[47]),
                             expression(PM[68]),expression(PM[96]),  expression(PM[152]))) 
dev.off()



# Supplemental Figure 7 Family:Qsd1 unique genetic variances as heritability ###### 
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Sup_H2_DifferentVarForFamilyQsd1.png', 1400, 800, res =120)
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
  mutate(trait = ifelse(trait == 'GE','Germiantion Percent','Germination Rate')) %>%
  ggplot(aes(TP, value, fill= Qsd1lab))+
  geom_bar(position = 'dodge',stat= 'identity')+
  facet_nested(trait+year~Family, scales= 'free_y') +scale_alpha_manual(values = c(.5,1)) +theme_bw()+
  labs(fill = 'Qsd1 status', y = 'Heritability', x = 'Timepoint')+
  ylim(0,1)
dev.off()




# Supplemental Figure 8 GI model training and prediction by All Qsd1 or ind Qsd1 ######
NoFlavia_inQsd1grp_RandTRN_RandTST_predictsOutgroup = read_csv('Analysis/GP_NoFlavia_Perqsd1GroupPre_withOutGroupPredictions.csv')
NoFlavia_RandonTRN_RandomTST =read_csv('Analysis/GP_NoFlavia_rndTest_rndtrain.csv')
names = c(expression(h[a]),expression(h[aD]),expression(h[aN]))

png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/GP_GI_trainQsd1TestQsd1.png', 1400, 1000, res =120)
HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(year == '2020/2021', TP %nin% c('TP1.5','TP2.5','TP3.5'), trait == 'GI')  %>%
  filter(substr(type, 1,2)=='h2', type != 'h2_Narrow_Qsd1Accounted') %>%
  mutate(maxCor = sqrt(value),
         type =   mapvalues(x = type, from = c('h2_Narrow','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D' ), 
                            to = c('sigma[a]','sigma[aN]','sigma[aD]'))) %>% 
  ggplot(aes(type, maxCor, fill = type))+
  geom_bar(stat = 'identity',position = 'dodge')+facet_grid(TP~'Germination Rate',margins = F)+
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
  mutate(Qsd1Train = mapvalues(x = Qsd1Train, from = c('0','2','AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles'),warn_missing = F),
         Qsd1Pred = mapvalues(x = Qsd1Pred, from = c('0','2', 'AllAlleles'), to = c('Non-dormant', 'Dormant','Both Alleles'),warn_missing = F),
         Trained_Predicted = paste0(Qsd1Train, ' \U2192 ',Qsd1Pred)) %>%
  ggplot(aes(x = n_train, y = correlation, color = Trained_Predicted, linetype = Trained_Predicted))+
  # geom_point(aes(shape = Trained_Predicted))+
  facet_grid(TP~ 'Germination Rate', scales = 'free') + theme_bw()+
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




# Supplemental Figure 9 GE prediction by Qsd1 overall of ind Qsd1 ######

png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Sup_GE_Qsd1trainQsd1Test.png', 1400, 800, res =120)
HeritabilitiesWinterGermTraits %>% filter(type != "H2_VarG_over_VarP") %>% 
  mutate(type = ifelse(type == 'H2_Broad_Qsd1Acounted','H2_Cullis_Qsd1Accounted', type)) %>%
  filter(year == '2020/2021', TP %nin% c('TP1.5','TP2.5','TP3.5','TP4','TP5'), trait == 'GE')  %>%
  filter(substr(type, 1,2)=='h2', type != 'h2_Narrow_Qsd1Accounted') %>%
  mutate(maxCor = sqrt(value),
         type =   mapvalues(x = type, from = c('h2_Narrow','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D' ), 
                            to = c('sigma[a]','sigma[aN]','sigma[aD]'))) %>% 
  ggplot(aes(type, maxCor, fill = type))+
  geom_bar(stat = 'identity',position = 'dodge')+facet_grid(TP~'Germination Percent',margins = F)+
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
  facet_grid(TP~ 'Germination Percent', scales = 'free') + theme_bw()+
  xlab('Number of training taxa') +
  labs(color = 'Qsd1 Trained \U2192 Predicted', linetype = 'Qsd1 Trained \U2192 Predicted',shape = 'Qsd1 Trained \U2192 Predicted') +
  geom_smooth(method = 'loess', se =T, size = 1.2,)+labs(y = 'Realized prediction accuracy')+
  scale_color_manual(values = c('salmon','salmon','salmon','mediumseagreen','mediumseagreen','cornflowerblue','cornflowerblue'))+
  # scale_linetype_manual(values = c('3131','dotted','solid','dotted','solid','dotted','solid'))+
  scale_linetype_manual(values = c('solid','dotted','3131','solid','dotted','dotted','solid'))+
  ylim(0,1)+
  plot_layout(ncol = 2, design = c('ABBBBB'))
dev.off()



# Supplemental Figure 10 Struct GE preds within Qsd1 #######
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/StructGESup.png', 1400, 900, res =120)
BetterPLotsSTGBLUP %>% filter(year == '2020_2021') %>% ungroup() %>%
  filter(trait == 'GE', Model != 'Base') %>%
  rbind(BetterPLotsSTGBLUP %>% filter(year == '2020_2021') %>% 
          filter(trait == 'GI', Model == 'Base')%>% ungroup() %>%
          rbind(.,.,.,.) %>% mutate(Model = rep(c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),each =15))) %>%
  mutate(Variation = factor(Effects, levels = c('Base','Base_fixef','JC','JC_Fix','NoFix','Full')),
         Model = mapvalues(Model, from = c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),
                           to = c('Family','Qsd1','Qsd1:Family interaction','Qsd1 + Family'))) %>% 
  # filter(Variation != 'JC') %>%
  ggplot(aes(x = TP, y = Variation, size = abs(PercentDif), color = PercentDif))+geom_point()+
  geom_text(aes(label =label), color = 'black', size = 3.88)+theme_bw()+
  labs(x = NULL, y = NULL)+ facet_grid(Qsd1~Model)+
  scale_y_discrete(label = c(
    expression(Full:~ y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(NoFix:~y%==% mu+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(JC_Fix:~ y%==% mu+bold(X)*bold(beta)+bold(Z[1])*bold(mu[1])),
    expression(JC:~y%==% mu+bold(Z[1])*bold(mu[1])),
    expression(Base_fixef:~y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])),
    expression(Base~GBLUP: ~y %==% mu + bold(Z[0]) * bold(mu[0]))
  ), limits = rev)+
  scale_colour_gradientn(name = 'Percent Difference\nfrom GBLUP', breaks = c(-15,0,15,65), 
                         colors = c('brown2', 'white', 'green', 'blue'),values = scales::rescale(c(-15,0,15,65)))+
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))+
  scale_size_continuous(name = NULL, range = c(0,15), guide='none')+
  geom_vline(xintercept =  c(1.5,2.5,3.5,4.5), alpha = 0.3)
dev.off()




# Supplemental Figure 11 Struct GI preds within Qsd1 #######
png('WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/StructGISup.png', 1400, 900, res =120)
BetterPLotsSTGBLUP %>% filter(year == '2020_2021') %>% ungroup() %>%
  filter(trait == 'GI', Model != 'Base') %>%
  rbind(BetterPLotsSTGBLUP %>% filter(year == '2020_2021') %>% 
          filter(trait == 'GI', Model == 'Base')%>% ungroup() %>%
          rbind(.,.,.,.) %>% mutate(Model = rep(c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),each =15))) %>%
  mutate(Variation = factor(Effects, levels = c('Base','Base_fixef','JC','JC_Fix','NoFix','Full')),
         Model = mapvalues(Model, from = c('Fam','Qsd1','Qsd1byFam','Qsd1PlusFam'),
                           to = c('Family','Qsd1','Qsd1:Family interaction','Qsd1 + Family'))) %>% 
  # filter(Variation != 'JC') %>%
  ggplot(aes(x = TP, y = Variation, size = abs(PercentDif), color = PercentDif))+geom_point()+
  geom_text(aes(label =label), color = 'black', size = 3.88)+theme_bw()+
  labs(x = NULL, y = NULL)+ facet_grid(Qsd1~Model)+
  scale_y_discrete(label = c(
    expression(Full:~ y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(NoFix:~y%==% mu+bold(Z[0])*bold(mu[0])+bold(Z[1])*bold(mu[1])),
    expression(JC_Fix:~ y%==% mu+bold(X)*bold(beta)+bold(Z[1])*bold(mu[1])),
    expression(JC:~y%==% mu+bold(Z[1])*bold(mu[1])),
    expression(Base_fixef:~y%==% mu+bold(X)*bold(beta)+bold(Z[0])*bold(mu[0])),
    expression(Base~GBLUP: ~y %==% mu + bold(Z[0]) * bold(mu[0]))
  ), limits = rev)+
  # scale_x_discrete(label = c(5,19,47,96, 152))+
  scale_colour_gradientn(name = 'Percent Difference\nfrom GBLUP', breaks = c(-15,0,15,65), 
                         colors = c('brown2', 'white', 'green', 'blue'),values = scales::rescale(c(-15,0,15,65)))+
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12))+
  scale_size_continuous(name = NULL, range = c(0,15), guide='none')+
  geom_vline(xintercept =  c(1.5,2.5,3.5,4.5), alpha = 0.3)
dev.off()


# Supplemental Figure 12 Segregation distortion tests #######
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

png(filename = 'WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/SegDistortionTest.png', width = 1000, height = 600, res = 120)
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








# Supplemental Figure 13 Markers responding to Qsd1 #####
png(filename = 'WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Loci1toLoci2.png', width = 700, height = 600, res = 120)
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



# Supplemental Figure 14 GE and GI TP1 model training and prediction by Qsd1 grouping.  ######
NOFlaviaDH2021_TP1GE_GI_PredictOutgroup = read_csv('Analysis/NOFlaviaDH2021_TP1GE_GI_PredictOutgroup.csv')
NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest = read_csv('Analysis/NOFlaviaDH2021_TP1GE_GI_RandTrain_RandTest.csv')

png(filename = 'WinterDHGerminationPaper/picsPNGforQsd1Effects_paper/Sup_GEGI2021TP1Predictions.png', width = 1000, height = 700, res = 120)
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
  mutate(trait = ifelse(trait == 'GE','Germination Percent','Germination Rate')) %>%
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














# Linkage decay across families: ##########
library(sommer)
WinterGD[1:5,1:5]
WinterGM[1:5,]
FamilyKey = AllDHBluesPerYear %>% select(taxa, Family) %>% unique()
LD.Overall = LD.decay(markers = WinterGD[,-1]-1, map = WinterGM %>% rename(Locus = SNP, LG = Chromosome),
                      silent = F, unlinked = FALSE)
LD.Flavia = LD.decay(markers = (WinterGD %>% filter(taxa %in% (FamilyKey %>% filter(Family %in% 'Flavia/DH130910'))$taxa))[,-1]-1, 
                     map = WinterGM %>% rename(Locus = SNP, LG = Chromosome),
                     silent = F, unlinked = FALSE)
LD.Wintmalt = LD.decay(markers = (WinterGD %>% filter(taxa %in% (FamilyKey %>% filter(Family %in% 'Wintmalt/DH130910'))$taxa))[,-1]-1, 
                       map = WinterGM %>% rename(Locus = SNP, LG = Chromosome),
                       silent = F, unlinked = FALSE)
LD.Scala = LD.decay(markers = (WinterGD %>% filter(taxa %in% (FamilyKey %>% filter(Family %in% 'Scala/DH130910'))$taxa))[,-1]-1, 
                    map = WinterGM %>% rename(Locus = SNP, LG = Chromosome),
                    silent = F, unlinked = FALSE)

LD.SY_Tepee= LD.decay(markers = (WinterGD %>% filter(taxa %in% (FamilyKey %>% filter(Family %in% 'SY_Tepee/DH130910'))$taxa))[,-1]-1, 
                      map = WinterGM %>% rename(Locus = SNP, LG = Chromosome),
                      silent = F, unlinked = FALSE)

drm(d ~ r2, data = LD.Overall$all.LG, fct = L.4())
LDmatrix = LD.Overall$all.LG
x$d = x$d/1000000

get_LD_distance = function(LDmatrix) {
  HW.st<-c(C=0.1)
  LDmatrix
  LDmatrix$d = LDmatrix$d/1000000
  HW.nonlinear<-nls(r2~((10+C*d)/((2+C*d)*(11+C*d)))*(1+((3+C*d)*(12+12*C*d+(C*d)^2))/(n*(2+C*d)*(11+C*d))),
                    data =  LDmatrix,
                    start=HW.st,control=nls.control(maxiter=100))
  tt<-summary(HW.nonlinear)
  new.rho<-tt$parameters[1]
  distance = seq(1:150)
  fpoints<-((10+new.rho*distance)/((2+new.rho*distance)*(11+new.rho*distance)))*(1+((3+new.rho*distance)*(12+12*new.rho*distance+(new.rho*distance)^2))/(n*(2+new.rho*distance)*(11+new.rho*distance)))
  length(fpoints[fpoints>0.2])
}

lapply(names(LD.Overall$by.LG), FUN = get_LD_distance)
get_LD_distance(LD.Overall$all.LG)
get_LD_distance(LD.Flavia$all.LG)
get_LD_distance(LD.Wintmalt$all.LG)
get_LD_distance(LD.Scala$all.LG)
get_LD_distance(LD.SY_Tepee$all.LG)


get_LD_distance(LD.Overall$by.LG$`5`)















