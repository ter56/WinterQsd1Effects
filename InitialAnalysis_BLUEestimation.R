# Script for the analysis, time series, genomic prediction, etc... Of the
# Cornell Winter Maling Barley DH families. 

library(sommer);library(arm);library(lme4);library(Hmisc);library(plyr);library(readxl);
library(tibble);library(patchwork);library(ggplot2);library(fda) ; library(magic); 
library(drc);library(rrBLUP);library(tidyr);library(ggh4x);library(dplyr);
setwd(rprojroot::find_rstudio_root_file())
load('WinterBarley/GenotypeData/WinterGD.RData')
load('WinterBarley/GenotypeData/WinterGM.RData')
WinterGM = WinterGM %>% arrange(Chromosome, Position)
WinterGD = WinterGD %>%
  mutate(taxa = gsub(pattern = '-',replacement = '_',taxa),
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         taxa1 = taxa) %>% remove_rownames()%>% 
  column_to_rownames('taxa1')
sum(WinterGM$SNP == colnames(WinterGD[,-1]))
# Make sure things are in the right order
# Sum should = 8384
dim(WinterGD)

# PCA plot of pop structure #####
WinterRelationship = rrBLUP::A.mat(WinterGD[,-1]-1, impute.method = 'EM', return.imputed = F)
WinterPCA = eigen(WinterRelationship)
WinterPVEPCA = WinterPCA$values/sum(WinterPCA$values)
data.frame(ordinal = 1:10, PVE = WinterPVEPCA[1:10]) %>%plot(., xlab = 'PC', col = 'red') 
winterlinePCAvalues = WinterPCA$vectors %>% data.frame()%>% 
  mutate(family = mapvalues(substr(WinterGD$taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','Sca','Win'), 
                            to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                   'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910')),
         taxa = WinterGD$taxa,
         shapes = ifelse(taxa %in% c('DH130910', 'Flavia','SY_Tepee','Wintmalt','Scala'), taxa, 'Lines'),
         size = ifelse(taxa %in% c('DH130910', 'Flavia','SY_Tepee','Wintmalt','Scala'), 3, 2))

winterlinePCAvalues %>% ggplot(aes(x = X1, y = X2, color = family)) + geom_point()+
  winterlinePCAvalues%>% ggplot(aes(x = X1, y = X3,color = family)) + geom_point()+
  winterlinePCAvalues%>% ggplot(aes(x = X2, y = X3,color = family)) + geom_point()

winterlinePCAvalues %>% filter(family != 'Cha') %>%
  ggplot(aes(x = X1, y = X2, color = family, shape = shapes)) + geom_point(aes(size = size))+theme_bw() +guides(size = "none")+
  xlab('PC1')+ylab('PC2')

# max 2 PCs accounts for population structure, though 1 would likely work
# Functions H2 get_blues #####
BLUPH2 = function(trait.lm) {
  ses<- se.ranef(trait.lm)$'taxa' #where 'm' is your model object from 'lmer' (replace 'genotypes' with whatever you call your individuals in the data)
  v_BLUP<- ses^2
  sigma2_g=VarCorr(trait.lm, comp="Variance")$'taxa'[1]
  Reliability<- 1- v_BLUP/ (2*sigma2_g)  #where sigma2_g is the genetic variance estimated
  H2<- round(mean(Reliability),3)
  return(H2)
}
##### Load raw data for the things #######
DHs2020 = rbind(read_excel("PhenotypeData/2020Harvest/DHtp1_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP1',PM_date =5),
                read_excel("PhenotypeData/2020Harvest/DHtp2_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP2',PM_date =19),
                read_excel("PhenotypeData/2020Harvest/DHtp3_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP3',PM_date =47),
                read_excel("PhenotypeData/2020Harvest/DHtp4_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP4',PM_date =96),
                read_excel("PhenotypeData/2020Harvest/DHtp5_all.xlsx", guess_max = 1723)%>%mutate(TP = 'TP5',PM_date =152)) %>%
  dplyr::mutate(GE =(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
                GI = 10*GE*(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+2*Day2Germ+3*Day3Germ),
                GI = ifelse(is.nan(GI),0,GI),
                GE5 = (Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ+KernelsLeft),
                GI5 = 10*GE5*(Day1Germ+Day2Germ+Day3Germ+Day4Germ+Day5Germ)/(Day1Germ+Day2Germ*2+Day3Germ*3+Day4Germ*4+Day5Germ*5),
                Location = ifelse(substr(PLOT,1,2)=='11' | substr(PLOT,1,2)=='a1','2020Snyder','2020Caldwell'),
                rep = as.factor(replication),
                taxa = mapvalues(Entry, from = c('Check 1','Check 2','Check 3','Check 4','Check 5','Check 6'),
                                 to = c('Flavia', 'Scala','DH130910','SY_Tepee','Wintmalt','Charles')),
                taxa =  mapvalues(taxa, from = c("Check 1-Flavia","Check 2-Scala","Check 3-DH130910","Check 3-SY Tepee",
                                                 "Check 4-SY Tepee","Check 5-Wintmalt","Check 6-Charles"),
                                  to = c('Flavia', 'Scala','DH130910','DH130910','SY_Tepee','Wintmalt','Charles')),
                taxa = gsub(pattern = '-',replacement = '_',taxa),
                taxa = gsub(pattern = ' ', replacement = '_',taxa),
                year = '2020',
                Family = mapvalues(substr(taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','Sca','Win','KWS'), 
                                   to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                          'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910','Scala/DH130910'))) %>%
  filter(Entry %nin% c('BBBdormSNLine', 'BBBdormSRLine'))

DHs2021 = rbind(read_excel('PhenotypeData/2021Phenotyping/Data/DHs_GGS_TP1_PM5_full.xlsx') %>% mutate(notes = NA, TP = 'TP1'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_12_full.xlsx')%>% mutate(notes = NA,TP = 'TP1.5'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_19_GGS_TP2_full.xlsx') %>% mutate(notes = NA, TP = 'TP2'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_33_GGS_TP3_full.xlsx') %>% mutate(TP = 'TP2.5'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_47_GGS_TP4_full.xlsx') %>% mutate(TP = 'TP3'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_68_GGS_TP5_full.xlsx')%>% mutate(TP='TP3.5'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_96_full.xlsx') %>% mutate(TP = 'TP4'),
                read_excel('PhenotypeData/2021Phenotyping/Data/DHs_PM_150_GGS_TP7_full.xlsx') %>% mutate(Pmdate = 152,notes = NA,TP = 'TP5')) %>%
  filter(!is.na(Day2Germ)) %>% filter(AssayNumber != 10 & AssayNumber != 318 & AssayNumber != 190) %>% #These look to be contamination
  filter(Location == 'McGowan' | Location == 'Ketola') %>%
  dplyr::select(!c('seed_mass', 'GerminationAssay', 'MaltingQualitySampling')) %>%
  dplyr::mutate(GE =(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+Day2Germ+Day3Germ+Day4Germ+KernelsLeft),
                GI = 10*GE*(Day1Germ+Day2Germ+Day3Germ)/(Day1Germ+2*Day2Germ+3*Day3Germ), Rep = as.factor(Rep)) %>%
  dplyr::rename('taxa' = Entry, 'PM_date' = Pmdate, 'rep' = Rep)  %>%
  mutate(Family = mapvalues(substr(taxa,1,3), from = c('BS6','BS7','BS8','BS9','DH1','Fla','SY_','KWS','Win'), 
                            to = c('Flavia/DH130910','Scala/DH130910','SY_Tepee/DH130910','Wintmalt/DH130910',
                                   'DH130910','Flavia/DH130910','SY_Tepee/DH130910','Scala/DH130910','Wintmalt/DH130910')),
         GI = ifelse(is.nan(GI),0,GI), 
         taxa = ifelse(AssayNumber == 408, 'Charles',taxa), #AssayNumber 408 is 'Charles' not the entry it claims to be!
         taxa = gsub(pattern = '-',replacement = '_',taxa), 
         taxa = gsub(pattern = ' ', replacement = '_',taxa),
         taxa = ifelse(taxa == 'KWS_Scala','Scala',taxa), year = '2021')

DHs2021 %>% filter(substr(taxa,1,2) !='BS') %>% select(taxa) %>% unique()

# 2020 GE and GI anova #####
anova(lm(GI~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==5))) #rep insignificant
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==19)))
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==47)))
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==96)))
anova(lm(GI ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==152)))

anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==5)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==19)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==47)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==96)))
anova(lm(GE ~ taxa +Location+rep, DHs2020 %>% filter(PM_date ==152)))

# 2021 GE and  GI anova #####
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==5)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==12)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==19)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==33)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==47)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==68)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==96)))
anova(lm(GE ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==152)))

anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==5)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==12)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==19)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==33)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==47)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==68)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==96)))
anova(lm(GI ~ taxa +Location+rep, DHs2021 %>% filter(PM_date ==152)))

# Location is pretty much always significant, rep bounces around. 
# 2020,2021 combined, only for TPs that are phenotypes at both years #####
DHs20202021 = DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date, year) %>%
  rbind(.,DHs2020 %>% select(taxa, rep, Location,TP, GE, GI, PM_date, year)) %>%
  mutate(year = factor(year,levels = c('2020','2021')))

anova(lm(GE ~ taxa +year+year %in% Location+rep, DHs20202021 %>% filter(PM_date ==5)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==19)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==47)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==96)))
anova(lm(GE ~ taxa +Location+rep+year, DHs20202021 %>% filter(PM_date ==152)))

anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==5)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==19)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==47)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==96)))
anova(lm(GI ~ taxa +Location+rep, DHs20202021 %>% filter(PM_date ==152)))
# Location very significant!
# Lets plot at the raw phentypes for the sake of the story: ######
DHs20202021 %>% pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% filter(trait =='GE') %>%
  ggplot(aes(x = value, fill = TP))+geom_density()+facet_nested(TP ~trait+year, scales = 'free')+
  theme_bw()+guides(fill = "none")+
  DHs20202021 %>% pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% filter(trait =='GI') %>%
  ggplot(aes(x = value, fill = TP))+geom_density()+facet_nested(TP ~trait+year, scales = 'free')+
  theme_bw()+ guides(fill = "none")+
  plot_layout(ncol = 2)

DHs20202021 %>%pivot_longer(cols = c(GE,GI), names_to = 'trait') %>% 
  ggplot(aes(x = PM_date, y = value, group = taxa))+geom_line()+facet_grid(trait~year, scales = 'free')

# Lets extract out 2020, 2021, and 2020/2021 values for GWA: Always including rep and location #####

BlueBlupsH2_Location_rep_taxa <- function(d, groupvars) {
  trait.lmer <- lmer(formula = value ~(1|taxa)+Location+rep, 
                     data = d)
  lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
    mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
  trait.lm = broom::tidy(lm(value~ taxa+Location+rep, data=d))
  
  first_taxa = d %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
  Intercept = trait.lm %>% filter(term == '(Intercept)') %>% select(estimate) %>% as.numeric()
  lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
    add_row(term = paste0('taxa',first_taxa),
            estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
    transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
              value = BLUE, type = 'BLUE')
  H2 = BLUPH2(trait.lmer)
  return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa) %>%
           join(d %>% select(taxa, Family) %>% unique()))
}

d2 = DHs2020 %>% select(taxa, rep, Location,TP, GE, GI,PM_date,year) %>%
  rbind(., DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year)) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%
  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>% filter(TP == 'TP2', trait == 'GE')

BlueBlupsH2_Year_rep_taxa <- function(d2, groupvars) {
  if (length(unique(d2$year))==2) {
    trait.lmer <- lmer(formula = value ~(1|taxa)+Location + rep, data = d2)
    # fixef(trait.lmer)[2]/need to think about this more. 
    
    Cept = (fixef(trait.lmer)[1]*4+sum(fixef(trait.lmer)[2:4]))/4
    lineEf = (ranef(trait.lmer)$taxa + Cept) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
      mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
    
    trait.lm = broom::tidy(lm(value~ taxa+Location+ rep, data=d2))
    first_taxa = d2 %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
    Intercept_list = trait.lm %>% filter(term == '(Intercept)'|substr(term,1,8)=='Location')
    Intercept = (Intercept_list$estimate[1]*4+sum(Intercept_list$estimate[2:4]))/4
    
    lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
      add_row(term = paste0('taxa',first_taxa),
              estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
      transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
                value = BLUE, type = 'BLUE')
    H2 = BLUPH2(trait.lmer)
    return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2')%>% join(d2 %>% select(taxa, Family)%>% unique()))
  }
  else {
    trait.lmer <- lmer(formula = value ~(1|taxa)+Location+rep, 
                       data = d2)
    lineEf = (ranef(trait.lmer)$taxa + fixef(trait.lmer)[1]) %>% as.data.frame() %>% rownames_to_column('taxa') %>% 
      mutate(type = 'BLUP') %>% rename(value = '(Intercept)')
    trait.lm = broom::tidy(lm(value~ taxa+Location+rep, data=d2))
    
    first_taxa = d2 %>% arrange(taxa) %>% slice_head( n = 1) %>% select(taxa) %>% as.character()
    Intercept = trait.lm %>% filter(term == '(Intercept)') %>% select(estimate) %>% as.numeric()
    lineBLUE = trait.lm %>% filter(substr(term,1,4)=='taxa') %>% 
      add_row(term = paste0('taxa',first_taxa),
              estimate = 0) %>% mutate(BLUE = estimate + Intercept) %>%
      transmute(taxa = gsub(pattern = 'taxa',replacement = '',x = term),
                value = BLUE,type = 'BLUE')
    H2 = BLUPH2(trait.lmer)
    return(rbind(lineEf, lineBLUE) %>% add_row(value = H2, type = 'H2') %>% arrange(type, taxa) %>%
             join(d2 %>% select(taxa, Family) %>% unique()))
  }
}

DH2020Estimates = DHs2020 %>% select(taxa, rep, Location,TP, GE, GI, PM_date,year, Family) %>% 
  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP, PM_date, trait, year) %>% group_modify(BlueBlupsH2_Location_rep_taxa) %>% ungroup()

DH2021Estimates = DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family) %>% pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP,PM_date, trait, year) %>% group_modify(BlueBlupsH2_Location_rep_taxa) %>% ungroup()

DHCombined = DHs2020 %>% select(taxa, rep, Location,TP, GE, GI,PM_date,year,Family) %>%
  rbind(., DHs2021 %>% select(taxa, rep, Location, TP, GE,GI,PM_date,year,Family)) %>% 
  mutate(year = factor(year, levels = c('2021','2020'))) %>%  pivot_longer(cols = c(GE, GI), names_to = 'trait') %>%
  group_by(TP,PM_date, trait) %>% group_modify(BlueBlupsH2_Year_rep_taxa)  %>% mutate(year = '2020/2021')

AllDHBluesPerYear = rbind(DH2020Estimates, DH2021Estimates,DHCombined) %>% filter(type =='BLUE') %>% ungroup()
save(AllDHBluesPerYear, file = 'WinterBarley/Analysis/AllDHBluesPerYear.RData')