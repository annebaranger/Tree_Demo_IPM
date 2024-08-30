library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
# get basic data
mean_pca<-tar_read(climate.cat)$FUNDIV_plotcat %>% 
  group_by(species,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2))
species.combination.select<- tar_read(species.combination.select) %>% 
  # mutate(pca1=(clim_up-clim_low)/2) %>% 
  select(species,species_combination,n_species,clim_id,wai,sgdd,clim_up,clim_low) %>% 
  left_join(mean_pca,by=c("species","clim_id"))


species.combination.select %>% 
  filter(species%in%c("Abies alba","Fagus sylvatica")) %>% 
  ggplot(aes(wai,sgdd,color=clim_id))+
  geom_point()+
  facet_wrap(~species)
# get disturbance metrics from mean 
disturbance_metric<-tar_read(disturbance_metric) %>% 
  mutate(elast="Abies_alba-amean",vr="mean")

# get disturbance metrics from elasticity analysis
tar_load(sim_disturbance_elast)
disturbance_metric_elast<-bind_rows(sim_disturbance_elast[grep( pattern = "_out", x = names(sim_disturbance_elast))]) %>% 
  select(-simul_eq_elast) %>% 
  rowwise() %>% 
  mutate(species_combination=gsub(pattern=elast,
                                  replacement=gsub(" ","_",species),
                                  x=species_combination),
         vr=str_split(elast,pattern="-",simplify = T)[2]) %>% 
  ungroup()


# compute mean elasticity per vital rates, for each species combi/clim
disturbance <- bind_rows(disturbance_metric,
                         disturbance_metric_elast) %>% 
  arrange(species,clim_id,species_combination,elast) %>% 
  pivot_longer(cols=c("resilience","recovery","resistance"),
               names_to="metric",
               values_to="metric_val") %>% 
  group_by(species,clim_id,species_combination,metric) %>% 
  mutate(dmetric=(metric_val[1]-metric_val)/metric_val[1]) %>%
  ungroup() %>%
  filter(vr!="mean") %>% 
  group_by(species,clim_id,species_combination,metric,vr) %>% 
  summarise(elast_mean=mean(dmetric,na.rm=TRUE),
            elast_sd=sd(dmetric,na.rm=TRUE)) %>% 
  ungroup()


# plots
disturbance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species=="Abies alba") %>% 
  # filter(clim_id!=10) %>% 
  ggplot(aes(pca1,elast_mean,color=metric))+
  geom_point()+
  geom_smooth()+
  # ylim(c(-0.4,0.4))+
  facet_grid(vr~as.factor(n_species),scales="free_y")
