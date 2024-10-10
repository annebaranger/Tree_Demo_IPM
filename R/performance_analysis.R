library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(lme4)
#### get basic data ####
#%%%%%%%%%%%%%%%%%%%%%%%
tar_load(species.list.ipm)
tar_load(species.list.disturbance)

## get demo traits
# tar_load(mean_demo)

# pca=prcomp(mean_demo[,c("inv_50","ba_equil")],scale=TRUE,center=TRUE)
# factoextra::fviz_pca_var(pca)

# mean_demo$pca1<-pca$x[,1]
# 
# 
traits<-read.csv("data/traits_complete.csv") %>%
  select(species,shade) %>%
  mutate(species=gsub(" ","_",species)) %>%
  rename(shade=shade)

# traits<-mean_demo %>% 
#   select(species,inv_50) %>% 
#   mutate(species=gsub(" ","_",species)) %>%
#   rename(shade=pca1)
  
# compute mean climate var by categories
mean_pca<-tar_read(climate.cat)$FUNDIV_plotcat %>% 
  group_by(species,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2)) %>% 
  ungroup()

# clean list of combinations

species.combination.excl<-tar_read(sim_forest_excl)$list.forests %>% 
  filter(forest.real) %>% 
  dplyr::select(species,clim_id,species_combination,excluded,competexcluded) %>% unique()

species.combination.select<- tar_read(species.combination.select) %>% 
  # mutate(pca1=(clim_up-clim_low)/2) %>% 
  left_join(mean_pca,by=c("species","clim_id")) %>% 
  left_join(species.combination.excl,by=c("species","clim_id","species_combination")) %>% 
  dplyr::select(species,species_combination,n_species,excluded,competexcluded,
                clim_id,pca1,pca2,wai,sgdd,clim_up,clim_low) 
  


#### analyse species exclusion ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species.combination.excl$smallcombi<-NA
for(i in 1:dim(species.combination.excl)[1]){
  if(species.combination.excl$competexcluded[i]!=""){
    s_p=gsub(" ","_",species.combination.excl$species[i])
    sp_combi=unlist(strsplit(species.combination.excl$species_combination[i], split = "\\."))
    # sp_combi=sp_combi[sp_combi!=s_p]
    sp_ex=unlist(strsplit(species.combination.excl$competexcluded[i],"\\."))[-1]
    
    new_combi=paste(sort(sp_combi[!sp_combi%in%sp_ex]),collapse=".")
    ex<-species.combination.excl %>% 
      filter(species==species.combination.excl$species[i]&
               clim_id==species.combination.excl$clim_id[i]&
               species_combination==new_combi)
    if(dim(ex)[1]==1){
      species.combination.excl$smallcombi[i]="present"
    }else{
      species.combination.excl$smallcombi[i]="absent"
    }
  }
  
}
# for checking the proportion of combination represented after species exclusion
# species.combination.excl %>% 
#   filter(competexcluded!="") %>% 
#   filter(smallcombi=="absent") 

#### get invasion metrics ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get disturbance metrics from mean 
disturbance_metric<-tar_read(disturbance_ba) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>%
  dplyr::select(species,elast,species_combination,clim_id,vr,resilience,resistance,recovery,matches(species.list.ipm))
invasion_metric<-tar_read(invasion_ba) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>% ## attention erruer
  dplyr::select(species,elast,species_combination,clim_id,vr,inv_mean,inv_max,inv_50,BA_100,BA_500,BA_1000,matches(species.list.ipm)) 


# compute mean elasticity per vital rates, for each species combi/clim
disturbance <- disturbance_metric %>% 
  arrange(species,clim_id,species_combination,elast) %>% 
  pivot_longer(cols=matches(species.list.ipm)) %>% 
  left_join(traits, by = c("name" = "species")) %>% 
  group_by(species,elast,species_combination,clim_id,vr) %>% 
  mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                               TRUE~0)) %>% 
  mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
         ba_target=sum(value*ba_multiple,na.rm = TRUE), 
         rel_ba_partner=ba_partner/(ba_target+ba_partner),
         nih=sum((sum(shade*ba_multiple)-shade)*value,na.rm=TRUE)/ba_partner,
         nid=sum(abs(sum(shade*ba_multiple)-shade)*value,na.rm=TRUE)/ba_partner) %>% 
  ungroup() %>% 
  dplyr::select(-c("ba_multiple","shade")) %>% 
  pivot_wider(names_from = name,
              values_from = value) %>% 
  dplyr::select(!matches(species.list.ipm)) 
invasion <- invasion_metric %>% ungroup() %>% 
  arrange(species,species_combination,clim_id) %>% 
  pivot_longer(cols=matches(species.list.ipm)) %>% 
  left_join(traits, by = c("name" = "species")) %>% 
  group_by(species,elast,species_combination,clim_id,vr) %>% 
  mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                               TRUE~0),
         value=case_when(!is.na(value)~1,
                         gsub("_"," ",name)==species~1,
                         TRUE~value)) %>%  # *1 because invasion are made at constant basal area of compet
  mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
         ba_target=sum(value*ba_multiple,na.rm = TRUE),
         rel_ba_partner=if_else(ba_partner!=0,
                                ba_partner/(ba_target+ba_partner),
                                0),
         nih=sum((sum(shade*ba_multiple)-shade)*value,na.rm=TRUE)/ba_partner, 
         nid=sum(abs(sum(shade*ba_multiple)-shade)*value,na.rm=TRUE)/ba_partner) %>% 
  
  ungroup() %>% 
  dplyr::select(-c("ba_multiple","shade")) %>% 
  pivot_wider(names_from = name,
              values_from = value) %>% 
  dplyr::select(!matches(species.list.ipm))%>% 
  left_join(species.combination.excl,by=c("species","clim_id","species_combination"))
ba_dif<-tar_read(ba_dif) %>% 
  select(!matches(c('ba_target',"ba_partner","n_species"))) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") 


# gather all data
performance <-invasion %>% 
  left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>% 
  left_join(ba_dif,by=c("species","elast","species_combination","clim_id","vr"))%>% 
  arrange(species,clim_id,species_combination,elast,vr) %>% 
  pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max",
                      "inv_50","BA_100","BA_500","BA_1000","ba_dif"),
               names_to="metric",
               values_to="metric_val") %>% 
  arrange(species,clim_id,elast,species_combination)%>% 
  mutate(ba_partner=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~ba_partner.x,
                              TRUE~ba_partner.y),
         ba_target=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~ba_target.x,
                             TRUE~ba_target.y),
         rel_ba_partner=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~rel_ba_partner.x,
                                  TRUE~rel_ba_partner.y),
         nih=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~nih.x,
                       TRUE~nih.y),
         nid=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~nid.x,
                       TRUE~nid.y)) %>% 
  dplyr::select(!matches("\\.x")) %>% dplyr::select(!matches(".y")) 
  # filter(excluded!="excluded") %>% 
  # filter(is.na(smallcombi)|smallcombi!="present")

#### plots ####
#%%%%%%%%%%%%%%

# variation of ba_dif with relative competition of partners
performance %>% 
  filter(metric=="ba_dif") %>% 
  ggplot(aes(rel_ba_partner,metric_val,color=excluded))+
  geom_point()

# covariations in inv performances
invasion %>% 
  ggplot(aes(BA_500,inv_50,color=species))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)
confint(lmer(inv_50~(1|species)+BA_100,data=invasion))

# variation of performances with climate and compet

performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="resilience") %>%
  filter(species%in% c("Abies alba","Carpinus betulus","Fagus sylvatica","Fraxinus excelsior",
                       "Picea abies","Quercus ilex","Quercus petraea","Quercus robur")) %>%
  filter(vr=="mean") %>%
  filter(is.na(smallcombi)|smallcombi=="absent") %>% 
  ggplot() +
  geom_line(aes(pca1,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca1,metric_val,color=ba_partner, group=interaction(elast,species_combination)),size=1)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15))+
  theme_bw()+
  facet_wrap(metric~species,ncol=4)+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Total basal area of competitors")



# see impact of different community composition on


# explore effect variation of competition with cliamte
performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="resilience") %>% #select one metric only with equil ba
  filter(vr=="mean") %>% 
  group_by(species,species_combination) %>% 
  filter(n()>4) %>%
  ggplot(aes(pca1,ba_partner,color=species_combination))+
  geom_point()+
  geom_smooth()+
  theme(legend.position = "none")+
  facet_wrap(~species)

performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="inv_50") %>% #select one metric only with equil ba
  filter(vr=="mean") %>% 
  filter(ba_partner!=0) %>% 
  group_by(species) %>%
  # ungroup() %>% 
  mutate(ba_partner=scale(ba_partner)) %>%
  ggplot(aes(pca1,ba_partner,color=species))+
  geom_point()+
  geom_smooth(method="lm",se = FALSE)

# show elasticity
performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species_combination%in%c("Abies_alba","Abies_alba.Picea_abies","Abies_alba.Fagus_sylvatica",
                                  "Abies_alba.Picea_abies.Pinus_sylvestris",
                                  "Fagus_sylvatica","Fagus_sylvatica.Quercus_petraea",
                                  "Fagus_sylvatica.Pinus_sylvestris","Fagus_sylvatica.Quercus_robur")) %>%
  filter(!metric%in%c("inv_max","inv_mean")) %>% 
  filter(species=="Abies alba") %>% 
  mutate(elast_combi=as.factor(paste0(elast,species_combination))) %>% 
  ggplot() +
  geom_line(aes(pca1,metric_val,color=species_combination, group=interaction(elast,species_combination)),linewidth=0.1)+
  geom_hline(yintercept = 0)+
  # scale_color_gradientn(colours = viridis(15),trans="log")+
  theme_bw()+
  facet_wrap(species~metric,scales="free_y",ncol=4)+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Species combination")

# just look how basal area changes on a single simul
# performance %>%
#   left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>%
#   # filter(species=="Abies alba") %>%
#   filter(!metric%in%c("inv_max","inv_mean")) %>%
#   filter(vr=="mean") %>%
#   filter(species_combination=="Abies_alba.Picea_abies") %>%
#   ggplot() +
#   geom_point(aes(pca1,metric_val,color=ba_partner, group=interaction(elast,species_combination)))+
#   geom_line(aes(pca1,metric_val,color=ba_partner, group=interaction(elast,species_combination)),linewidth=1)+
#   geom_hline(yintercept = 0)+
#   scale_color_gradientn(colours = viridis(15),trans="log")+
#   facet_wrap(species~metric,scales="free_y",ncol=4)



# covariation between performances
performance %>% 
  filter(species_combination%in%c("Abies_alba","Abies_alba.Picea_abies","Abies_alba.Fagus_sylvatica",
                                  "Abies_alba.Picea_abies.Pinus_sylvestris",
                                  "Fagus_sylvatica","Fagus_sylvatica.Quercus_petraea",
                                  "Fagus_sylvatica.Pinus_sylvestris","Fagus_sylvatica.Quercus_robur")) %>%
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species=="Fagus sylvatica") %>%
  filter(!metric%in%c("inv_max","inv_mean")) %>% 
  filter(vr=="mean") %>% 
  filter(!grepl("Acer_pseudoplatanus",species_combination)) %>% 
  mutate(elast_combi=as.factor(paste0(elast,species_combination)),
         species_combination=factor(species_combination,levels=c("Fagus_sylvatica","Abies_alba.Fagus_sylvatica",
                                                                 "Fagus_sylvatica.Quercus_petraea","Fagus_sylvatica.Pinus_sylvestris",
                                                                 "Fagus_sylvatica.Quercus_robur"))) %>% 
  group_by(metric) %>% 
  mutate(metric_val=scale(metric_val,scale = TRUE)) %>% 
  ggplot() +
  geom_point(aes(pca1,metric_val,color=metric))+
  geom_smooth(aes(pca1,metric_val,color=metric),se = FALSE)+
  geom_hline(yintercept = 0)+
  facet_wrap(~species_combination,ncol=3)+
  theme_bw()+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Type of performance")




# Each metric as fonction of competition, clim is color 
performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(!metric%in%c("inv_max","inv_mean")) %>% 
  filter(vr=="mean") %>% 
  mutate(elast_combi=as.factor(paste0(elast,species_combination))) %>% 
  ggplot() +
  geom_point(aes(log(ba_partner+1),metric_val,color=clim_id,group=clim_id))+
  geom_line(aes(log(ba_partner+1),metric_val,color=clim_id,group=clim_id))+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15),trans="log")+
  facet_wrap(species~metric,scales="free_y",ncol=4)




### Analysis of resilience ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# figure that helps
performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="resilience") %>%
  filter(species%in% c("Abies alba","Carpinus betulus","Fagus sylvatica","Fraxinus excelsior",
                       "Picea abies","Quercus ilex","Quercus petraea","Quercus robur")) %>%
  filter(vr=="mean") %>%
  filter(is.na(smallcombi)|smallcombi=="absent") %>%
  ggplot() +
  geom_line(aes(pca1,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca1,metric_val,color=ba_partner, group=interaction(elast,species_combination)),size=1)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15))+
  theme_bw()+
  facet_wrap(metric~species,ncol=4)+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Total basal area of competitors")

# analysis of species exclusion
data_exclu<-performance %>% 
  filter(metric=="resilience") %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  mutate(species=forcats::fct_reorder(species, shade),
         simul_state=case_when(!is.na(smallcombi)~"CompetitorExclusion",
                               excluded=="excluded"~"SpeciesExclusion",
                               TRUE~"SpeciesCoex")) %>% 
  group_by(species,clim_id) %>% 
  mutate(n_comb=n()) %>% 
  group_by(species,clim_id,simul_state) %>% 
  mutate(prop=n()/n_comb,
         nih_mean=mean(nih,na.rm=TRUE)) %>% 
  dplyr::select(species,clim_id,simul_state,n_comb,prop,nih_mean,shade) %>% unique() %>% 
  ungroup() 
data_exclu%>% 
  ggplot(aes(clim_id,prop,fill=simul_state))+
  geom_col()+
  facet_wrap(~species)


performance %>% 
  filter(metric=="resilience") %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  group_by(species) %>% 
  mutate(nih_mean=mean(nih,na.rm=TRUE)) %>% 
  dplyr::select(species,nih_mean,shade) %>% unique() %>% 
  ungroup() %>% ggplot(aes(shade,nih_mean))+geom_point()+geom_smooth()
model=lm(prop~clim_id*nih_mean+clim_id:species,data=subset(data_exclu,simul_state=="CompetitorExclusion"))
summary(model)
model$model
  
#### model test ####
#%%%%%%%%%%%%%%%%%%%

# test effect of competition
data.compet<- performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="resilience") %>% 
  group_by(species) %>% 
  mutate(ba_tot=ba_target+ba_partner,
         ba_tot=scale(ba_tot,scale=TRUE)) 
summary(lmer(ba_tot~pca1+(pca1|species),data=data.compet))

performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(vr=="mean") %>% 
  # filter(species=="Abies alba") %>% 
  filter(metric=="inv_50")->data.test
data.test %>% ggplot(aes(ba_partner,metric_val,color=pca1))+geom_point()+facet_wrap(~species)+scale_color_gradientn(colours = viridis(n=11))
summary(lmer(metric_val~(1|species)+(1|species:species_combination)+pca1*ba_partner,data=data.test))
