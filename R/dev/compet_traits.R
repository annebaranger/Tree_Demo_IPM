library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(lme4)
library(mice)
#### get basic data ####
#%%%%%%%%%%%%%%%%%%%%%%%
tar_load(species.list.ipm)
tar_load(species.list.disturbance)
sp_list=gsub("_"," ",species.list.ipm)
sp_list[sp_list=="Betula"]<-"Betula pendula"
sp_list=c(sp_list,"Betula pubescens")

shade=read.csv2("data/data_Niinemets&Valladares_2006.csv") |> 
  rename(species=Species) |> 
  filter(species %in% sp_list) |> 
  mutate(shade_tolerance.mean=as.numeric(shade_tolerance.mean)) |> 
  select(species,shade_tolerance.mean)

try_data<-read.csv2("data/try/try_query.csv")

try_format<-try_data %>% 
  filter(!is.na(TraitID)) %>% 
  filter(TraitID!=18) %>% 
  select(AccSpeciesName,ObservationID,ObsDataID,TraitID,TraitName,StdValue) %>% 
  unique() %>% 
  mutate(StdValue=as.numeric(StdValue))
# table(try_format[try_format$TraitID==3116,]$UnitName) # checked if units were consistent across traits

leaf_trait<-try_format %>%
  mutate(species=case_when(AccSpeciesName=="Betula pendula"~"Betula",
                           AccSpeciesName=="Betula pubescens"~"Betula",
                           TRUE~AccSpeciesName)) %>% 
  group_by(species,TraitName) %>% 
  filter(!(TraitName=="Seed dry mass"&StdValue>quantile(StdValue,probs=0.95,na.rm=TRUE))) %>% 
  filter(!(TraitName=="Seed dry mass"&StdValue<quantile(StdValue,probs=0.05,na.rm=TRUE))) %>% 
  summarise(trait_mean=median(StdValue,na.rm=TRUE),
            trait_sd=sd(StdValue,na.rm=TRUE)) %>% 
  mutate(trait_mean=case_when(is.nan(trait_mean)~NA,
                              TRUE~trait_mean),
         trait_sd=case_when(is.nan(trait_sd)~NA,
                              TRUE~trait_sd),
         cv=trait_sd/trait_mean) %>% 
  ungroup() %>% 
  rename(trait=TraitName)

# get max height from FUNDIV data
tar_load(FUNDIV_data)
species_meta<-tar_read(species_meta) %>% 
  mutate(species=gsub("_"," ",species)) %>% 
  select(species,taxa)

max_height<- FUNDIV_data %>% unique() %>% 
  filter(species%in%c(sp_list,"Betula")) %>% 
  filter(height1>0) %>% 
  group_by(species) %>% 
  filter(!is.na(height1)) %>% 
  filter(!is.na(weight1)) %>% 
  summarise(trait_mean=weighted.quantile(height1,w=weight1,prob=0.99),
            trait="max_height") 

mean_trait<-rbind(leaf_trait[,c("species","trait","trait_mean")],
                  max_height) %>% 
  left_join(species_meta) %>% 
  left_join(shade,by="species") %>% 
  pivot_wider(names_from = trait,
              values_from = trait_mean)
colnames(mean_trait)<-c("species","taxa","shade","SLA","LN","LT","SDM","WD","HM")
mean_trait<-mean_trait[,c("species","taxa","shade","SLA","LN","LT","WD","HM")]
wooddensity<-read.csv2("data/try/GlobalWoodDensityDatabase.csv") %>% 
  filter(Binomial%in%sp_list) %>% 
  mutate(species=case_when(Binomial=="Betula pendula"~"Betula",
                           Binomial=="Betula pubescens"~"Betula",
                           TRUE~Binomial)) %>% 
  group_by(species) %>% 
  summarise(WD=mean(Wood.density))

wooddensity<-data.frame(species=gsub("_"," ",species.list.ipm)) %>% 
  left_join(wooddensity,by='species') %>%
  left_join(mean_trait[,c("species","WD")],by="species") %>% 
  mutate(WD=case_when(is.na(WD.x)~WD.y,
                      TRUE~WD.x))

mean_trait$WD<-wooddensity$WD
imputed_data <- mice(mean_trait, method = 'pmm', m = 5, maxit = 50, seed = 500)
trait_complete <- complete(imputed_data, 1)  # '1' refers to the first imputed dataset

imputed_data_cf <- mice(mean_trait[mean_trait$taxa=="conifer",], method = 'pmm', m = 5, maxit = 50, seed = 500)
trait_complete_cf <- complete(imputed_data_cf, 1)  # '1' refers to the first imputed dataset

trait_complete<-rbind(trait_complete_br,trait_complete_cf)
pca_trait<- prcomp(trait_complete[,c("WD","HM","shade")], center = TRUE, scale = TRUE)
library(factoextra)
fviz_pca_var(pca_trait, repel = TRUE)


load("data/try/seedTraits.rdata")


trait_complete$pca1=pca_trait$x[,1]
trait_complete$pca2=pca_trait$x[,2]

write.csv(trait_complete,"data/traits_complete.csv")
