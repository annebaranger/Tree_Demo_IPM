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
  group_by(species,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(pca_sc=scale(pca1))

performance_mod<-performance |> 
  left_join(mean_pca, by=c("species","clim_id")) |> 
  left_join(species.combination,by=c("species","clim_id","species_combination"))

## plots 
performance_mod |>
  filter(metric=="ba_dif") |>
  # filter(gsub(" ","_",species)!=species_combination) |> 
  # filter(excluded=="notExcluded") |> 
  ggplot(aes(pca1,ba_partner,color=log(prop_all)))+
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_wrap(~species)


performance_mod |> 
  filter(metric=="resilience") |> 
  ggplot(aes(clim_id,prop_all))+
  geom_point()+
  facet_wrap(~species)
## model
mod_compet<-lmer(ba_partner ~ (1|species)+(1|species:species_combination)+ pca1 + I(pca1^2), 
                 weights = prop_all,
                 data=subset(performance_mod,metric=="resilience"))
summary(mod_compet)
confint(mod_compet)

sim_res <- simulateResiduals(mod_compet)
plot(sim_res)


mod_data<-performance |> 
  filter(metric=="resilience") |> 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) |> 
  mutate(pca12=pca1^2)
mod_data <- na.omit(mod_data)
# submodel for climate and competition
mod_data |> ggplot(aes(ba_partner,pca1,color=species_combination))+
  geom_point()+
  geom_smooth(method="lm")+
  theme(legend.position = "none")+facet_wrap(~species,scales="free")
mod_compet<-lmer(nid~(1|species)+(1|species:species_combination)+
                   pca1,data=mod_data)
confint(mod_compet)


mod_resi<-lmer(metric_val ~ (1|species)+(1|species:species_combination)+ pca1 + pca12 +
                 pca1*ba_partner + pca1*nih , data=mod_data)
confint(mod_resi)
summary(mod_resi)
# -- Make model
mod_sem = list(
  lmer(ba_partner ~ (1|species)+(1|species:species_combination)+ pca1 + I(pca1^2), 
       weights = prop_all,
       data=subset(performance_mod,metric=="ba_dif")),
  # lmer(nih~(1|species)+(1|species:species_combination)+
  #       pca1,data=mod_data), 
  lmer(metric_val ~  (1|species)+(1|species:species_combination)+ pca1 + I(pca1^2) + pca1*ba_partner, 
     data=subset(performance_mod,metric=="ba_dif")))



# Extract the effects
boot_sem = bootEff(mod_sem, R = 5, seed = 13, parallel = "no", ran.eff = "species")

piecewiseSEM:::plot.psem(
  piecewiseSEM::as.psem(mod_sem),
  node_attrs = data.frame(shape = "rectangle", color = "black", fillcolor = "grey"),
  layout = "tree"
)
# Output list
out = list()
out$sem = mod_sem
out$boot = boot_sem

# Return output