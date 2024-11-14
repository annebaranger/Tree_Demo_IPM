library(targets)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(glmmTMB)
library(semEff)


tar_load(performance)
species.combination<-tar_read(species.combination) |> 
  dplyr::select(species,clim_id,species_combination,n_species,n,is.sim,is.sp,prop,prop_all)



mean_pca<-tar_read(climate.cat)$FUNDIV_plotcat %>% 
  group_by(species) |> 
  mutate(pca1_opt=mean(pca1)) |> 
  group_by(species,pca1_opt,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(pca_sc=as.numeric(scale(pca1)),
         dist_pca1=abs(pca1-pca1_opt)/(max(pca1)-min(pca1))) |> 
  ungroup()

mean_pca |> ggplot(aes(dist_pca1,species))+geom_boxplot()

performance_mod<-performance |> 
  # filter(gsub(" ","_",species)!=species_combination) |> 
  left_join(mean_pca, by=c("species","clim_id")) |> 
  left_join(species.combination,by=c("species","clim_id","species_combination")) |>
  dplyr::select(species,clim_id,species_combination,excluded,competexcluded,metric, metric_val,ba_target,ba_partner,pca1,pca1_opt,dist_pca1,prop_all) |> 
  filter(metric%in%c("resilience","ba_dif","BA_100","inv_50")) |> 
  group_by(species,clim_id,species_combination) |> 
  arrange(species,clim_id,species_combination,metric) |> 
  mutate(ba_partner=case_when(metric%in%c("BA_100","inv_50")~ba_partner[2],
                                 TRUE~ba_partner),
         ba_target=case_when(metric%in%c("BA_100","inv_50")~ba_target[2],
                              TRUE~ba_target),
         metric_val=case_when(metric=="resilience"~log(metric_val*ba_target),
                              TRUE~metric_val)) |> 
  ungroup() |> 
  group_by(metric) |> 
  mutate(metric_val=as.numeric(scale(metric_val))) |>
  ungroup() |> 
  pivot_wider(names_from="metric",
              values_from = "metric_val") |> 
  mutate(pca1=as.numeric(scale(pca1)),
         dist_pca1=as.numeric(scale(dist_pca1)),
         ba_partner=as.numeric(scale(ba_partner)),
         ba_target=as.numeric(scale(ba_target)),
         pca1_2=pca1^2) 
  
  # filter(species%in%c("Abies alba","Carpinus betulus","Fagus sylvatica","Quercus ilex","Picea abies","Quercus robur"))

## plots ####
#%%%%%%%%%%%%
performance_mod |>
  # filter(gsub(" ","_",species)!=species_combination) |> 
  # filter(excluded=="notExcluded") |> 
  ggplot(aes(dist_pca1,ba_partner,color=log(prop_all)))+
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_wrap(~species)

performance_mod |>
  filter(metric=="ba_dif") |>
  ggplot(aes(pca1,ba_partner,color=log(prop_all)))+
  geom_point()+
  geom_smooth(se=FALSE)


performance_mod |> 
  filter(metric=="resilience") |> 
  ggplot(aes(clim_id,prop_all))+
  geom_point()+
  facet_wrap(~species)


## model ####
#%%%%%%%%%%%%
# mod_compet<-lmer(ba_partner ~ (1|species)+(1|species:species_combination)+ pca1 + I(pca1^2), 
#                  weights = prop_all,
#                  data=subset(performance_mod,metric=="resilience"))
# summary(mod_compet)
# confint(mod_compet)
# 
# sim_res <- simulateResiduals(mod_compet)
# plot(sim_res)
# 

mod_sem = list(
  lmer(ba_partner ~ (1|species) + pca1 + dist_pca1, 
       weights = prop_all,
       data=performance_mod),
  # lmer(inv_50 ~ (1|species) + (1|species:species_combination)+ pca1 + dist_pca1 + ba_partner + ba_target + BA_100,
  #      data=performance_mod),
  lmer(inv_50 ~ (1|species)+ pca1 + dist_pca1 + ba_partner + ba_target ,
       data=performance_mod),
  lmer(resilience ~ (1|species) + pca1 + dist_pca1 + ba_partner + ba_target ,
       data=performance_mod),
  lmer(ba_target ~ (1|species) + pca1 + dist_pca1 + ba_partner,
       data=performance_mod))

# Extract the effects
boot_sem = bootEff(mod_sem, R = 100, seed = 13, parallel = "no", ran.eff = "species")

psem_model <- as.psem(mod_sem)
summary(psem_model)
piecewiseSEM:::plot.psem(
  piecewiseSEM::as.psem(mod_sem),
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "grey"),
  layout = "tree"
)


out = list()
out$sem = mod_sem
out$boot = boot_sem

GGally::ggpairs(performance_mod[,c("inv_50","ba_dif","resilience")])
performance |> ggplot(aes(metric_val))+geom_density()+facet_wrap(~metric,scales="free")
