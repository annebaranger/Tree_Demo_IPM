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

# get invasion metrics
invasion_metric_elast<-tar_read(invasion_metric_elast) %>% 
  select(species,elast,species_combination,clim_id,inv_mean,inv_max,inv_50) %>% 
  rowwise() %>% 
  mutate(species_combination=gsub(pattern=elast,
                                  replacement=gsub(" ","_",species),
                                  x=species_combination),
         vr=str_split(elast,pattern="-",simplify = T)[2]) %>% 
  ungroup()


# get disturbance metrics from mean 
disturbance_metric<-tar_read(disturbance_metric) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>%
  select(species,elast,species_combination,clim_id,vr,resilience,resistance,recovery)
invasion_metric<-tar_read(invasion_metric) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>%
  select(species,elast,species_combination,clim_id,vr,inv_mean,inv_max,inv_50) 
# get disturbance metrics from elasticity analysis
tar_load(disturbance_metric_elast)
disturbance_metric_elast<-disturbance_metric_elast %>%
  select(species,elast,species_combination,clim_id,resilience,resistance,recovery) %>% 
  rowwise() %>% 
  mutate(species_combination=gsub(pattern=elast,
                                  replacement=gsub(" ","_",species),
                                  x=species_combination),
         vr=str_split(elast,pattern="-",simplify = T)[2]) %>% 
  ungroup()


# compute mean elasticity per vital rates, for each species combi/clim
disturbance <- bind_rows(disturbance_metric,
                         disturbance_metric_elast) %>% 
  arrange(species,clim_id,species_combination,elast) 
invasion <- bind_rows(invasion_metric,
                      invasion_metric_elast) %>% ungroup() %>% 
  arrange(species,species_combination,clim_id) 
elasticity<- invasion %>% 
  left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>% 
  arrange(species,clim_id,species_combination,elast,vr) %>% 
  pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max","inv_50"),
               names_to="metric",
               values_to="metric_val") %>% 
  group_by(species,clim_id,species_combination,metric) %>% 
  mutate(dmetric=(metric_val[1]-metric_val)/metric_val[1]) %>%
  ungroup() %>%
  filter(vr!="mean") %>% 
  group_by(species,clim_id,species_combination,metric,vr) %>% 
  summarise(elast_mean=mean(dmetric,na.rm=TRUE),
            elast_sd=sd(dmetric,na.rm=TRUE)) %>% 
  ungroup() %>% 
  filter(metric%in%c("resilience","recovery","resistance","inv_50"))

performance <-invasion %>% 
  left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>% 
  arrange(species,clim_id,species_combination,elast,vr) %>% 
  pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max","inv_50"),
               names_to="metric",
               values_to="metric_val") %>% 
  arrange(species,clim_id,elast,species_combination)
# plots
elasticity %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species=="Abies alba") %>% 
  # filter(clim_id!=10) %>% 
  ggplot(aes(pca1,elast_mean,color=metric))+
  geom_point()+
  geom_smooth()+
  # ylim(c(-0.4,0.4))+
  facet_grid(vr~as.factor(n_species),scales="free_y")


performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species=="Fagus sylvatica") %>%
  group_by(clim_id,pca1,elast,n_species,metric) %>% 
  summarise(metric_val=mean(metric_val)) %>% 
  ggplot(aes(pca1,metric_val,color=elast)) +
  geom_line()+
  theme(legend.position = "none")+
  facet_grid(metric~n_species,scales="free_y")
  
## add traits of competitors
