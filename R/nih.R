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

# compute mean climate var by categories
mean_pca<-tar_read(climate.cat)$FUNDIV_plotcat %>% 
  group_by(species,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2))

# clean list of combinations
species.combination.select<- tar_read(species.combination.select) %>% 
  # mutate(pca1=(clim_up-clim_low)/2) %>% 
  select(species,species_combination,n_species,clim_id,wai,sgdd,clim_up,clim_low) %>% 
  left_join(mean_pca,by=c("species","clim_id"))

shade=read.csv2("data/data_Niinemets&Valladares_2006.csv") |> 
  rename(species=Species) |> 
  filter(species %in% gsub("_"," ",species.list.disturbance)) |> 
  mutate(s_p=gsub(" ","_",species),
         shade_tolerance.mean=as.numeric(shade_tolerance.mean)) |> 
  select(s_p,shade_tolerance.mean)
#### get invasion metrics ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get disturbance metrics from mean 
disturbance_metric<-tar_read(disturbance_ba) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>%
  select(species,elast,species_combination,clim_id,vr,resilience,resistance,recovery,matches(species.list.ipm))
invasion_metric<-tar_read(invasion_ba) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>% ## attention erruer
  select(species,elast,species_combination,clim_id,vr,inv_mean,inv_max,inv_50,matches(species.list.ipm)) 


# compute mean elasticity per vital rates, for each species combi/clim
invasion<-invasion_metric %>%
  arrange(species, clim_id, species_combination, elast) %>%
  pivot_longer(cols = matches(species.list.ipm)) %>%
  left_join(shade, by = c("name" = "s_p")) %>% 
  filter(!(is.na(value)&(name!=gsub(" ","_",species)))) %>%
  mutate(s_p = gsub(" ", "_", species),
         value=if_else(is.na(value),0,value)) %>% 
  left_join(shade, by = c("s_p")) |> 
  rename(shade_target = shade_tolerance.mean.y,
         shade_partner = shade_tolerance.mean.x) %>% 
  group_by(species, elast, species_combination, clim_id, vr) %>%
  mutate(ba_multiple = if_else(name == gsub(" ", "_", species), 1, 0),
         ba_partner = sum(value * abs(1 - ba_multiple), na.rm = TRUE),
         ba_target = sum(value * ba_multiple, na.rm = TRUE),
         rel_ba_partner = ba_partner / (ba_target + ba_partner),
         nih = if_else(name == species_combination, 0, 
                       sum((shade_target - shade_partner) * (value)) / ba_partner),
         nid = if_else(name == species_combination, 0, 
                       sum(abs(shade_target - shade_partner) * (value)) / ba_partner)) %>%
  ungroup() %>% 
  select(-c("ba_multiple","shade_partner")) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  select(!matches(species.list.ipm)) 


disturbance<-disturbance_metric %>%
  arrange(species, clim_id, species_combination, elast) %>%
  mutate(id_metric=row_number()) |> 
  pivot_longer(cols = matches(species.list.ipm)) %>%
  left_join(shade, by = c("name" = "s_p")) %>% 
  filter(!(is.na(value)&(name!=gsub(" ","_",species)))) %>%
  mutate(s_p = gsub(" ", "_", species),
         value=if_else(is.na(value),0,value)) %>% 
  left_join(shade, by = c("s_p")) |> 
  rename(shade_target = shade_tolerance.mean.y,
         shade_partner = shade_tolerance.mean.x) %>% 
  group_by(species, elast, species_combination, clim_id, vr) %>%
  mutate(ba_multiple = if_else(name == gsub(" ", "_", species), 1, 0),
         ba_partner = sum(value * abs(1 - ba_multiple), na.rm = TRUE),
         ba_target = sum(value * ba_multiple, na.rm = TRUE),
         rel_ba_partner = ba_partner / (ba_target + ba_partner),
         nih = if_else(name == species_combination, 0, 
                       sum((shade_target - shade_partner) * (value)) / ba_partner),
         nid = if_else(name == species_combination, 0, 
                       sum(abs(shade_target - shade_partner) * (value)) / ba_partner)) %>%
  ungroup() %>% 
  select(-c("ba_multiple","shade_partner")) %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  select(!matches(species.list.ipm)) 


performance <-invasion %>% 
  left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>% 
  arrange(species,clim_id,species_combination,elast,vr) %>% 
  pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max","inv_50"),
               names_to="metric",
               values_to="metric_val") %>% 
  arrange(species,clim_id,elast,species_combination)%>% 
  mutate(ba_partner=case_when(grepl("inv",metric)~ba_partner.x,
                              TRUE~ba_partner.y),
         ba_target=case_when(grepl("inv",metric)~ba_target.x,
                             TRUE~ba_target.y),
         rel_ba_partner=case_when(grepl("inv",metric)~rel_ba_partner.x,
                                  TRUE~rel_ba_partner.y),
         nih=case_when(grepl("inv",metric)~nih.x,
                                  TRUE~nih.y),
         nid=case_when(grepl("inv",metric)~nid.x,
                       TRUE~nid.y)) %>% 
  select(!matches(".x")) %>% select(!matches(".y")) |> 
  unique()

#### plot ####
#%%%%%%%%%%%%%

performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="inv_50") %>% 
  filter(vr=="mean") %>%
  filter(species %in% gsub("_"," ",species.list.disturbance)) %>%   
  # group_by(species,species_combination,metric) %>% 
  # filter(n()>8) %>%
  ggplot() +
  geom_line(aes(pca1,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca1,metric_val,color=nid, group=interaction(elast,species_combination)),size=1)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(50),trans="log")+
  theme_bw()+
  facet_wrap(species~metric,scale="free",ncol=4)+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Total basal area of competitors")


#### model ####
#%%%%%%%%%%%%%
library(semEff)
