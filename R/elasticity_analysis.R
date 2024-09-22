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

# compute mean climate var by categories
mean_pca<-tar_read(climate.cat)$FUNDIV_plotcat %>% 
  group_by(species,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2))

# clean list of combinations
species.combination.select<- tar_read(species.combination.select) %>% 
  # mutate(pca1=(clim_up-clim_low)/2) %>% 
  select(species,species_combination,n_species,clim_id,wai,sgdd,clim_up,clim_low) %>% 
  left_join(mean_pca,by=c("species","clim_id"))

#### get invasion metrics ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get disturbance metrics from mean 
disturbance_metric<-tar_read(disturbance_ba) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>%
  select(species,elast,species_combination,clim_id,vr,resilience,resistance,recovery,matches(species.list.ipm))
invasion_metric<-tar_read(invasion_ba) %>% 
  mutate(elast="Abies_alba-amean",vr="mean") %>% ## attention erruer
  select(species,elast,species_combination,clim_id,vr,inv_mean,inv_max,inv_50,matches(species.list.ipm)) 


# get disturbance metrics from elasticity analysis
tar_load(disturbance_ba_elast)
disturbance_ba_elast<-disturbance_ba_elast %>%
  select(species,elast,species_combination,clim_id,resilience,resistance,recovery,matches(species.list.ipm)) %>% 
  rowwise() %>% 
  mutate(species_combination=gsub(pattern=elast,
                                  replacement=gsub(" ","_",species),
                                  x=species_combination),
         vr=str_split(elast,pattern="-",simplify = T)[2]) %>% 
  ungroup()
invasion_metric_elast<-tar_read(invasion_ba_elast) %>% 
  select(species,elast,species_combination,clim_id,inv_mean,inv_max,inv_50,matches(species.list.ipm)) %>% 
  rowwise() %>% 
  mutate(species_combination=gsub(pattern=elast,
                                  replacement=gsub(" ","_",species),
                                  x=species_combination),
         vr=str_split(elast,pattern="-",simplify = T)[2]) %>% 
  ungroup()

# compute mean elasticity per vital rates, for each species combi/clim
disturbance <- bind_rows(disturbance_metric,
                         disturbance_ba_elast) %>% 
  arrange(species,clim_id,species_combination,elast) %>% 
  pivot_longer(cols=matches(species.list.ipm)) %>% 
  group_by(species,elast,species_combination,clim_id,vr) %>% 
  mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                               TRUE~0)) %>% 
  mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
         ba_target=sum(value*ba_multiple,na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(-ba_multiple) %>% 
  pivot_wider(names_from = name,
              values_from = value) %>% 
  select(!matches(species.list.ipm))
invasion <- bind_rows(invasion_metric,
                      invasion_metric_elast) %>% ungroup() %>% 
  arrange(species,species_combination,clim_id) %>% 
  pivot_longer(cols=matches(species.list.ipm)) %>% 
  group_by(species,elast,species_combination,clim_id,vr) %>% 
  mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                                TRUE~0)) %>% 
  mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
         ba_target=sum(value*ba_multiple,na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(-ba_multiple) %>% 
  pivot_wider(names_from = name,
              values_from = value) %>% 
  select(!matches(species.list.ipm))

elasticity<- invasion %>% 
  left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>% 
  arrange(species,clim_id,species_combination,elast,vr) %>% 
  pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max","inv_50"),
               names_to="metric",
               values_to="metric_val") %>% 
  filter(!grepl("Fagus_sylvatica.Picea.abies",species_combination)) %>% 
  mutate(ba_partner=case_when(grepl("inv",metric)~ba_partner.x, 
                              TRUE~ba_partner.y),
         ba_target=case_when(grepl("inv",metric)~ba_target.x,
                              TRUE~ba_target.y)) %>% 
  group_by(species,clim_id,species_combination,metric) %>% 
  mutate(dmetric=100*(metric_val[1]-metric_val)/metric_val[1]) %>%
  ungroup() %>%
  filter(vr!="mean") %>% 
  group_by(species,clim_id,species_combination,metric,vr) %>% 
  summarise(elast_mean=mean(dmetric),
            elast_sd=sd(dmetric),
            ba_partner_mean=mean(ba_partner),
            ba_target_mean=mean(ba_target)) %>% 
  ungroup() %>% 
  filter(metric%in%c("resilience","recovery","resistance","inv_50"))

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
                             TRUE~ba_target.y)) %>% 
  select(!matches(".x")) %>% select(!matches(".y"))


#### plots ####
#%%%%%%%%%%%%%%

elasticity %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species=="Abies alba") %>% 
  filter(ba_partner_mean<100) %>% 
  # filter(clim_id!=10) %>% 
  ggplot(aes(pca1,elast_mean,color=ba_partner_mean,group=species_combination))+
  geom_point()+
  geom_smooth(aes(pca1,elast_mean,color=ba_partner_mean,group=species_combination),
              method="gam")+
  # geom_smooth(aes(group=species_combination))+
  scale_color_gradientn(colours = viridis(15),trans="log")+
  # ylim(c(-0.4,0.4))+
  geom_hline(yintercept = 0)+
  facet_wrap(vr~metric,scales="free_y")

# example of vital rates sensitivity, for abies alba, inv_50
elasticity %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species=="Abies alba") %>% 
  filter(metric=="inv_50") %>% 
  filter(ba_partner_mean<100) %>% 
  # filter(clim_id!=10) %>% 
  ggplot(aes(pca1,elast_mean,color=ba_partner_mean,group=species_combination))+
  geom_line(color="grey")+
  geom_point()+
  # geom_smooth(aes(group=species_combination))+
  scale_color_gradientn(colours = viridis(15),trans="log")+
  # ylim(c(-0.4,0.4))+
  geom_hline(yintercept = 0)+
  facet_wrap(vr~metric,scales="free_y")+
  theme_bw()+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Mean elasticity of performance to vital rate",
       color="Total basal area of competitors",
       title="Abies alba")

# variation of performances with climate and compet

performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(species_combination%in%c("Abies_alba","Abies_alba.Picea_abies","Abies_alba.Fagus_sylvatica",
                                  "Abies_alba.Picea_abies.Pinus_sylvestris",
                                  "Fagus_sylvatica","Fagus_sylvatica.Quercus_petraea",
                                  "Fagus_sylvatica.Pinus_sylvestris","Fagus_sylvatica.Quercus_robur")) %>%
  filter(!metric%in%c("inv_max","inv_mean")) %>% 
  filter(vr=="mean") %>%
  mutate(elast_combi=as.factor(paste0(elast,species_combination))) %>% 
  ggplot() +
  geom_line(aes(pca1,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca1,metric_val,color=ba_partner, group=interaction(elast,species_combination)),size=0.5)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15),trans="log")+
  theme_bw()+
  facet_wrap(species~metric,scales="free_y",ncol=4)+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Total basal area of competitors")

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

# model test
performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(vr=="mean") %>% 
  # filter(species=="Abies alba") %>% 
  filter(metric=="inv_50")->data.test
data.test %>% ggplot(aes(ba_partner,metric_val,color=pca1))+geom_point()+facet_wrap(~species)+scale_color_gradientn(colours = viridis(n=11))
summary(lmer(metric_val~(1|species)+pca1*ba_partner,data=data.test))


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


## add traits of competitors



#open simul of abal.piab.pisy & abal.piab for the same climate
# app<-readRDS("rds/Abies_alba/clim_1/sim_disturbance/Abies_alba.Picea_abies.Pinus_sylvestris.rds")
# ap<-readRDS("rds/Abies_alba/clim_1/sim_disturbance/Abies_alba.Picea_abies.rds")
# readRDS("rds/Fagus_sylvatica/clim_5/sim_disturbance/Fagus_sylvatica.rds")->sim
# 
# sim %>% filter(var=="BAsp",species!="Picea_abies") %>% 
#   ggplot(aes(time,value,color=species))+
#   geom_line()+
#   xlim(c(490,503))
# sim %>% filter(var=="H",time==500)


tar_load(sim_dist.id)
tar_load(sim_disturbance)
sim_disturbance<-sim_disturbance[grepl("Fagus_sylvatica/",sim_disturbance)]
readRDS(sim_disturbance[1])

out<-setNames(data.frame(matrix( nrow = 0,ncol=9)),
              nm=c("species","var","time","mesh","size","equil","value","sim","ba_partner"))
for(i in 1:length(sim_disturbance)){
  print(i)
  species_partner=str_remove(strsplit(sim_disturbance[i],"/")[[1]][[5]],".rds")
  clim=as.numeric(strsplit(strsplit(sim_disturbance[i],"/")[[1]][[3]],"_")[[1]][[2]])
  sp=strsplit(sim_disturbance[i],"/")[[1]][[2]]
  ba_partner=disturbance %>% 
    filter(species==gsub("_"," ",sp)&
           species_combination==species_partner&
           clim_id==clim&
           grepl("amean",elast)) %>% 
    pull(ba_partner)
  sim.i<-readRDS(sim_disturbance[i]) %>% 
    filter(species=="Fagus_sylvatica") %>% 
    filter(var=="BAsp") %>% 
    mutate(sim=i,
           ba=ba_partner)
  out<-rbind(out,sim.i)
}
out %>% 
  filter(var=="BAsp") %>% 
  ggplot(aes(time,value,color=ba,group=as.factor(sim)))+
  geom_line()
