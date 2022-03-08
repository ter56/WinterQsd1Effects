
library(Hmisc);library(plyr);library(readxl);
library(tibble);
library(rrBLUP);library(tidyr);
library(lme4);
library(asreml)
library(dplyr)

load('GenotypeData/WinterGD.RData')
load('GenotypeData/WinterGM.RData')
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

DHsRaw = DHs2021 %>% dplyr::select(year, Family, taxa, rep, Location, TP, GE,GI,PM_date) %>%
  rbind(.,DHs2020 %>% dplyr::select(year, Family,taxa, rep, Location,TP, GE, GI, PM_date)) %>%
  rbind(DHs2021 %>% dplyr::select(year, Family, taxa, rep, Location, TP, GE,GI,PM_date) %>%
          rbind(.,DHs2020 %>% dplyr::select(year, Family,taxa, rep, Location,TP, GE, GI, PM_date))%>%
          mutate(year = '2020/2021'))%>%
  mutate(year = factor(year,levels = c('2020','2021','2020/2021')),
         Location = as.factor(Location),
         taxa = as.factor(taxa),
         Family = factor(Family)) %>% join(WinterGD[,c('taxa','Qsd1')]) %>%
  filter(Family %nin% c('Cha','End','DH130910'))

DHsRaw_Qsd1Filt =  DHs2021 %>% dplyr::select(year, Family, taxa, rep, Location, TP, GE,GI,PM_date) %>%
  rbind(.,DHs2020 %>% dplyr::select(year, Family,taxa, rep, Location,TP, GE, GI, PM_date)) %>%
  rbind(DHs2021 %>% dplyr::select(year, Family, taxa, rep, Location, TP, GE,GI,PM_date) %>%
          rbind(.,DHs2020 %>% dplyr::select(year, Family,taxa, rep, Location,TP, GE, GI, PM_date))%>%
          mutate(year = '2020/2021')) %>%
  mutate(year = factor(year,levels = c('2020','2021','2020/2021')),
         Location = as.factor(Location),
         taxa = as.factor(taxa),
         Family = factor(Family)) %>% join(WinterGD[,c('taxa','Qsd1')]) %>% filter(Qsd1 != 1) %>%
  mutate(Qsd1 = factor(Qsd1),FamilyQsd1 = factor(paste0(Family,Qsd1)))%>%
  filter(Family %nin% c('Cha','End','DH130910'))

DHs20202021 = DHs2021 %>% select(taxa,Family, rep, Location, TP, GE,GI,PM_date, year) %>%
  rbind(.,DHs2020 %>% select(taxa,Family, rep, Location,TP, GE, GI, PM_date, year)) %>%
  mutate(year = factor(year,levels = c('2020','2021'))) %>% filter(Family %nin% c('Cha','End','DH130910')) %>%
  mutate(Location = as.factor(Location),
         taxa = as.factor(taxa),
         Family = factor(Family))
### Narrow and broad sense eritability ######
df = DHsRaw_Qsd1Filt %>% pivot_longer(cols = c(GE,GI), values_to = 'value',names_to = 'trait') %>% 
  group_by(year, TP, trait) %>% filter(year =='2020/2021',TP =='TP1',trait =='GI') %>% arrange(Qsd1)
NarrowAndBroadSenseH2 = function(df, groupvars, GeneticData = WinterGD){
  GDsub = GeneticData %>% filter(taxa %in% df$taxa) %>% arrange(taxa)
  GMat = rrBLUP::A.mat(GDsub[,-1]-1, impute.method = 'EM', return.imputed = F)
  # GMatInv =solve(GMat + diag(.001, nrow(GMat))) #adding smaell
  # attr(GMatInv, "INVERSE") <- TRUE
  # attr(GMatInv, "colNames") <- as.character(df$taxa)
  df = df %>% filter(taxa %in% GDsub$taxa) %>% arrange(taxa)
  print(groupvars)
  additiveModel = tryCatch(expr = suppressWarnings(asreml(fixed = value~Location+rep, 
                                                          random = ~vm(taxa,GMat, singG = "PSD"),
                                                          residual = ~idv(units),data = df)),
                           error = function(e){NA})
  h2Narrow = tryCatch(expr = suppressWarnings(summary(additiveModel)$varcomp[1,1]/(summary(additiveModel)$varcomp[1,1]+summary(additiveModel)$varcomp[2,1])),
                      error = function(e){NA})
  
  BroadSenseModel = asreml(fixed = value ~ Location+rep, random = ~idv(taxa), residual = ~idv(units), data = df)
  H2BroadVargOverVarP = summary(BroadSenseModel)$varcomp[1,1]/(summary(BroadSenseModel)$varcomp[1,1]+summary(BroadSenseModel)$varcomp[2,1])
  H2Cullis = mean(1 - (summary(BroadSenseModel, coef = TRUE)$coef.random[,2])^2/(2*summary(BroadSenseModel)$varcomp[1,1]))
  
  BroadSenseAferQsd1 = asreml(fixed = value~Location+rep+Qsd1, random = ~idv(taxa), residual = ~idv(units), data = df)
  BroadSenseAferQsd1H2 =  summary(BroadSenseAferQsd1)$varcomp[1,1]/
    (summary(BroadSenseAferQsd1)$varcomp[1,1]+summary(BroadSenseAferQsd1)$varcomp[2,1])
  
  additiveModel_qsd1 = tryCatch(expr = suppressWarnings(asreml(fixed = value~Location+rep+Qsd1,
                                                               random = ~vm(taxa,GMat, singG = "PSD"), residual = ~idv(units),data = df)),
                                error = function(e){NA})
  
  h2Narrow_qsd1Accounted = tryCatch(expr = suppressWarnings(summary(additiveModel_qsd1)$varcomp[1,1]/
                                                              (summary(additiveModel_qsd1)$varcomp[1,1]+summary(additiveModel_qsd1)$varcomp[2,1])),
                                    error = function(e){NA})
  H2NarrowSensePerQSD1 = tryCatch(expr = suppressWarnings(
    asreml(fixed = value~Location+rep+Qsd1, random = ~at(Qsd1):vm(taxa,GMat, singG = "PSD"), residual = ~idv(units),data = df)),
    error = function(e){NA})
  summary(H2NarrowSensePerQSD1)
  
  H2NarrowPerQsd1N = 
    tryCatch(expr = suppressWarnings(summary(H2NarrowSensePerQSD1)$varcomp[1,1]/
                                       (summary(H2NarrowSensePerQSD1)$varcomp[1,1]+summary(H2NarrowSensePerQSD1)$varcomp[3,1])),
             error = function(e){NA})
  
  H2NarrowPerQsd1D = 
    tryCatch(expr = suppressWarnings(summary(H2NarrowSensePerQSD1)$varcomp[2,1]/
                                       (summary(H2NarrowSensePerQSD1)$varcomp[2,1]+summary(H2NarrowSensePerQSD1)$varcomp[3,1])),
             error = function(e){NA})
  print(H2NarrowPerQsd1D)
  return(data.frame(type = c('h2_Narrow', 'H2_Cullis', 'H2_VarG_over_VarP', 'H2_Broad_Qsd1Acounted', 'h2_Narrow_Qsd1Accounted','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D'), 
                    value = c(h2Narrow,H2Cullis, H2BroadVargOverVarP,BroadSenseAferQsd1H2, h2Narrow_qsd1Accounted,H2NarrowPerQsd1N,H2NarrowPerQsd1D)))
}
NarrowAndBroadSenseH2v2 = function(df, groupvars, GeneticData = WinterGD){
  GDsub = GeneticData %>% filter(taxa %in% df$taxa) %>% arrange(taxa)
  GMat = rrBLUP::A.mat(GDsub[,-1]-1, impute.method = 'EM', return.imputed = F)
  # GMatInv =solve(GMat + diag(.001, nrow(GMat))) #adding smaell
  # attr(GMatInv, "INVERSE") <- TRUE
  # attr(GMatInv, "colNames") <- as.character(df$taxa)
  df = df %>% filter(taxa %in% GDsub$taxa) %>% arrange(Qsd1, taxa)
  print(groupvars)
  additiveModel = tryCatch(expr = suppressWarnings(asreml(fixed = value~Location+rep, 
                                                          random = ~vm(taxa,GMat, singG = "PSD"),
                                                          residual = ~idv(units),data = df)),
                           error = function(e){NA})
  h2Narrow = tryCatch(expr = suppressWarnings(summary(additiveModel)$varcomp[1,1]/(summary(additiveModel)$varcomp[1,1]+summary(additiveModel)$varcomp[2,1])),
                      error = function(e){NA})
  
  BroadSenseModel = asreml(fixed = value ~ Location+rep, random = ~idv(taxa), residual = ~idv(units), data = df)
  H2BroadVargOverVarP = summary(BroadSenseModel)$varcomp[1,1]/(summary(BroadSenseModel)$varcomp[1,1]+summary(BroadSenseModel)$varcomp[2,1])
  H2Cullis = mean(1 - (summary(BroadSenseModel, coef = TRUE)$coef.random[,2])^2/(2*summary(BroadSenseModel)$varcomp[1,1]))
  
  BroadSenseAferQsd1 = asreml(fixed = value~Location+rep+Qsd1, random = ~idv(taxa), residual = ~idv(units), data = df)
  BroadSenseAferQsd1H2 =  summary(BroadSenseAferQsd1)$varcomp[1,1]/
    (summary(BroadSenseAferQsd1)$varcomp[1,1]+summary(BroadSenseAferQsd1)$varcomp[2,1])
  
  additiveModel_qsd1 = tryCatch(expr = suppressWarnings(asreml(fixed = value~Location+rep+Qsd1,
                                                               random = ~vm(taxa,GMat, singG = "PSD"), residual = ~idv(units),data = df)),
                                error = function(e){NA})
  
  h2Narrow_qsd1Accounted = tryCatch(expr = suppressWarnings(summary(additiveModel_qsd1)$varcomp[1,1]/
                                                              (summary(additiveModel_qsd1)$varcomp[1,1]+summary(additiveModel_qsd1)$varcomp[2,1])),
                                    error = function(e){NA})
  H2NarrowSensePerQSD1 = tryCatch(expr = suppressWarnings(
    asreml(fixed = value~Location+rep+Qsd1, random = ~at(Qsd1):vm(taxa,GMat, singG = "PSD"), residual = ~dsum(~idv(units)|Qsd1),data = df)),
    error = function(e){NA})
  summary(H2NarrowSensePerQSD1)
  
  H2NarrowPerQsd1N = 
    tryCatch(expr = suppressWarnings(summary(H2NarrowSensePerQSD1)$varcomp[1,1]/
                                       (summary(H2NarrowSensePerQSD1)$varcomp[1,1]+summary(H2NarrowSensePerQSD1)$varcomp[4,1])),
             error = function(e){NA})
  
  H2NarrowPerQsd1D = 
    tryCatch(expr = suppressWarnings(summary(H2NarrowSensePerQSD1)$varcomp[2,1]/
                                       (summary(H2NarrowSensePerQSD1)$varcomp[2,1]+summary(H2NarrowSensePerQSD1)$varcomp[6,1])),
             error = function(e){NA})
  print(H2NarrowPerQsd1D)
  return(data.frame(type = c('h2_Narrow', 'H2_Cullis', 'H2_VarG_over_VarP', 'H2_Broad_Qsd1Acounted', 'h2_Narrow_Qsd1Accounted','h2_NarrowU_Qsd1N','h2_NarrowU_Qsd1D'), 
                    value = c(h2Narrow,H2Cullis, H2BroadVargOverVarP,BroadSenseAferQsd1H2, h2Narrow_qsd1Accounted,H2NarrowPerQsd1N,H2NarrowPerQsd1D)))
}

# Function to get genetic variance per family, Account for Family effects, get genetic variance per family:Qsd1 #####
GetFamilyMeansAndVarianceWinter = function(df, groupvars) {
  idvTaxa.asm = asreml(fixed=value~Location+rep+Family, random = ~idv(taxa), residual = ~idv(units), data = df)
  FMV.asr = asreml(fixed =   value~Location+rep+Family, random = ~at(Family):taxa, residual = ~idv(units), data = df)
  Family.fixed = summary(FMV.asr, coef = TRUE)$coef.fixed %>%as.data.frame() %>% rownames_to_column(var = 'Coef')
  Family.random = summary(FMV.asr)$varcomp %>% as.data.frame() %>% rownames_to_column(var = 'Coef')
  
  df.returnModel1 = Family.fixed %>% transmute(type = 'Fixed', value = solution, Coef = Coef, std.error = `std error`, z.ratio = z.ratio, bound = NA, `%ch` = NA) %>%
    rbind(Family.random %>% transmute(type = 'Variance', value = component, Coef = Coef, std.error = std.error, z.ratio = z.ratio, bound = bound, `%ch` = `%ch`))%>%
    mutate(model ='base')
  FMV.asrQsd1 = asreml(fixed = value~Location+rep+Family+Family:Qsd1, random = ~at(Family):taxa, residual = ~idv(units), data = df)
  Family.fixed2 = summary(FMV.asrQsd1, coef = TRUE)$coef.fixed %>%as.data.frame() %>% rownames_to_column(var = 'Coef')
  
  Family.random2 = summary(FMV.asrQsd1)$varcomp %>% as.data.frame() %>% rownames_to_column(var = 'Coef')
  
  df.returnModel2 = Family.fixed2 %>% transmute(type = 'Fixed', value = solution, Coef = Coef, std.error = `std error`, z.ratio = z.ratio, bound = NA, `%ch` = NA) %>%
    rbind(Family.random2 %>% transmute(type = 'Variance', value = component, Coef = Coef, std.error = std.error, z.ratio = z.ratio, bound = bound, `%ch` = `%ch`)) %>%
    mutate(model = 'Family:Qsd1_Fixed')
  
  VarAccounted = Family.random$component - Family.random2$component
  to_return = rbind(df.returnModel1,df.returnModel2) %>% 
    add_row(type = 'VarianceExplained', Coef = Family.random$Coef[1:4], value = round(VarAccounted,3)[1:4],model = 'difference')%>%
    add_row(model = 'lrt idv(taxa) to at(Family):taxa', value = lrt(idvTaxa.asm,FMV.asr)$`Pr(Chisq)`, Coef = 'lrt1')
  
  return(to_return)
}
# Model with: fixed = value~Location+rep+Family+Family:Qsd1, random = ~at(FamilyQsd1):taxa, residual = ~idv(units) and lrt #####
At_FamilyQsd1_taxa = function(df, groupvars){
  
  FMV.asrQsd1 =   asreml(fixed = value~Location+rep+Family+Family:Qsd1, random = ~at(Family):taxa, residual = ~idv(units), data = df)
  
  FamQsd1Taxa.asr=asreml(fixed = value~Location+rep+Family+Family:Qsd1, random = ~at(FamilyQsd1):taxa, residual = ~idv(units), data = df)
  
  Family.fixed3 = summary(FamQsd1Taxa.asr, coef = TRUE)$coef.fixed %>%as.data.frame() %>% rownames_to_column(var = 'Coef')
  Family.random3 = summary(FamQsd1Taxa.asr)$varcomp %>% as.data.frame() %>% rownames_to_column(var = 'Coef')
  
  df.returnModel3 = Family.fixed3 %>% transmute(type = 'Fixed', value = solution, Coef = Coef, std.error = `std error`, z.ratio = z.ratio, bound = NA, `%ch` = NA) %>%
    rbind(Family.random3 %>% transmute(type = 'Variance', value = component, Coef = Coef, std.error = std.error, z.ratio = z.ratio, bound = bound, `%ch` = `%ch`))%>%
    mutate(model ='atFamilyQsd1_taxa')%>%
    add_row(model = 'lrt at(Family):taxa to at(FamilyQsd1):taxa', value = lrt(FMV.asrQsd1,FamQsd1Taxa.asr)$`Pr(Chisq)`,
            Coef = 'lrt')
  
  lrt(FMV.asrQsd1,FamQsd1Taxa.asr)$`Pr(Chisq)`
  return(df.returnModel3)
}

# Winter running of things and save output Save OUtput #####
HeritabilitiesWinterGermTraits = DHsRaw_Qsd1Filt %>% pivot_longer(cols = c(GE,GI), values_to = 'value',names_to = 'trait') %>% 
  group_by(year, TP, trait) %>% arrange(Qsd1)%>%
  group_modify(NarrowAndBroadSenseH2)
HeritabilitiesWinterGermTraitsv2 = DHsRaw_Qsd1Filt %>% pivot_longer(cols = c(GE,GI), values_to = 'value',names_to = 'trait') %>% 
  group_by(year, TP, trait) %>% arrange(Qsd1)%>%
  group_modify(NarrowAndBroadSenseH2v2)

FamilyMeansAndVarianceQsd1 = DHsRaw_Qsd1Filt %>% pivot_longer(cols = c(GE,GI), values_to = 'value',names_to = 'trait') %>% group_by(year, TP, trait) %>%
  group_modify(GetFamilyMeansAndVarianceWinter)
atFamQsd1TaxaVar = DHsRaw_Qsd1Filt %>% pivot_longer(cols = c(GE,GI), values_to = 'value',names_to = 'trait') %>% group_by(year, TP, trait) %>%
  group_modify(At_FamilyQsd1_taxa)
save(HeritabilitiesWinterGermTraits, file = 'Analysis/OutputAsreml/HeritabilitiesWinterGermTraits.RData')
save(FamilyMeansAndVarianceQsd1, file = 'Analysis/OutputAsreml/FamilyMeansAndVarianceQsd1Asreml.RData')
save(atFamQsd1TaxaVar, file = 'Analysis/OutputAsreml/atFamQsd1TaxaVar.RData')
save(HeritabilitiesWinterGermTraitsv2, file = 'Analysis/OutputAsreml/HeritabilitiesWinterGermTraitsv2.RData')
