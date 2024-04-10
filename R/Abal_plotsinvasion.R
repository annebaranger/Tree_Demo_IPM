library(matreex)
library(dplyr)
library(tidyr)
library(ggplot2)
library(targets)

# get ba eq 
sim_forest_list=tar_read(sim_forest_list)$list.forest
tar_load(sim_equil)
invasion_metric_ba=invasion_metric |> 
  mutate(ba_init=NA)
for (i in 1:dim(invasion_metric)[1]){
  print(i)
  sp=invasion_metric[i,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=invasion_metric[i,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=invasion_metric[i,"ID.spclim"][[1]]
  
  if(length(species.in)>1){
    simul_eq.partner=sim_forest_list |> 
      filter(species==sp &
               species_combination==sub(paste0(s_p,"\\."),"",species.comb) &
               ID.spclim == clim) |> 
      pull(simul_eq)
    
    # Read the simulation at equilibrium
    sim_equilibrium.partner.in = readRDS(sim_equil[simul_eq.partner])
    
    BAeq=sum(sim_equilibrium.partner.in |> filter(var=="BAsp",equil) |> pull(value))
    invasion_metric_ba$ba_init[i]=BAeq
  }else{
    # simul_eq.species=sim_forest_list |> 
    #   filter(species==sp &
    #            species_combination==s_p &
    #            ID.spclim == clim) |> 
    #   pull(simul_eq)
    # sim_equilibrium.species.in = readRDS(sim_equil[simul_eq.species])
    # equil.i = sim_equilibrium.species.in %>%
    #   filter(var == "n", equil) %>%  
    #   mutate(value=case_when(size>0~0,
    #                          TRUE~value)) |>
    #   mutate(BAtot=sum((100/2000)^2*pi*value)) |> View() 
    #   pull(value)
    invasion_metric_ba$ba_init[i]=0
  }
  
}

invasion_metric_ba |>
  rowwise() |> 
  mutate(ncomb=length(unlist(strsplit(species_combination,"[.]")))) |> 
  ungroup() |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  ggplot(aes(ba_init,value,color=as.factor(ncomb)))+geom_point()+
  facet_wrap(~inv,scales="free")


invasion_metric_ba |>
  rowwise() |> 
  mutate(ncomb=length(unlist(strsplit(species_combination,"[.]")))) |> 
  ungroup() |> 
  ggplot(aes(as.factor(ncomb),ba_init))+geom_boxplot()

# explo invasion_rate
inv_test=invasion_metric |>
  rowwise() |> 
  mutate(ncomb=length(unlist(strsplit(species_combination,"[.]")))) |> 
  ungroup()
inv_test|> 
  filter(ncomb!=1) |>
  ggplot(aes(wai_id,sgdd_id,color=inv_50))+
  geom_point(size=6,alpha=2)+
  facet_wrap(~ncomb)+
  scale_color_gradientn(colours = viridis(15))


inv_test |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  group_by(ncomb,wai_id,inv) |> 
  summarize(mean.inv=mean(value)) |> 
  ggplot(aes(ncomb,mean.inv,color=inv))+
  geom_point()+
  facet_wrap(~wai_id)

inv_test |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  group_by(ncomb,wai_id,sgdd_id,inv) |> 
  summarize(mean.inv=mean(value)) |>
  filter(inv=="inv_50") |> 
  ggplot(aes(wai_id,mean.inv,color=sgdd_id))+
  geom_point()+
  facet_wrap(~ncomb)+
  scale_color_gradientn(colours = viridis(15))


inv_test |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  group_by(ncomb,wai_id,sgdd_id,inv) |> 
  summarize(mean.inv=mean(value)) |> 
  ggplot(aes(ncomb,mean.inv,color=inv))+
  geom_point()+
  facet_grid(sgdd_id~wai_id)


summary(lm(inv_50~ncomb+wai+sgdd,data=inv_test))
car::Anova(lm(inv_50~ncomb+wai+sgdd,data=inv_test))


tar_load(disturbance_metric)
dist_test=disturbance_metric |>
  rowwise() |> 
  mutate(ncomb=length(unlist(strsplit(species_combination,"[.]")))) |> 
  ungroup()
dist_test|> 
  filter(ncomb!=1) |> 
  ggplot(aes(wai_id,sgdd_id,color=recovery))+
  geom_point(size=6,alpha=2)+
  facet_wrap(~ncomb)+
  scale_color_gradientn(colours = viridis(15))


dist_test |> 
  pivot_longer(cols=c("recovery","resistance","resilience"),names_to = "resil") |> 
  group_by(ncomb,wai_id,resil) |> 
  summarize(mean.inv=mean(value)) |> 
  ggplot(aes(ncomb,mean.inv,color=resil))+
  geom_point()+
  facet_wrap(~wai_id)

dist_test |> 
  pivot_longer(cols=c("recovery","resistance","resilience"),names_to = "resil") |> 
  group_by(ncomb,wai_id,sgdd_id,resil) |> 
  summarize(mean.inv=mean(value)) |>
  filter(resil=="resilience") |> 
  ggplot(aes(sgdd_id,mean.inv,color=wai_id))+
  geom_point()+
  facet_wrap(~ncomb)+
  scale_color_gradientn(colours = viridis(15))


dist_test |> 
  pivot_longer(cols=c("recovery","resistance","resilience"),names_to = "resil")|> 
  group_by(ncomb,wai_id,sgdd_id,resil) |> 
  summarize(mean.inv=mean(value)) |> 
  ggplot(aes(ncomb,mean.inv,color=resil))+
  geom_point()+
  facet_grid(sgdd_id~wai_id)

dist_test |> 
  pivot_longer(cols=c("recovery","resistance","resilience"),names_to = "resil")|> 
  group_by(resil) |>
  mutate(value=scale(value)) |>
  ungroup() |>
  ggplot(aes(wai_id,sgdd_id,color=value))+
  geom_point(size=6)+
  facet_grid(resil~ncomb)+
  scale_color_gradientn(colours = viridis(15))

dist_test |> 
  pivot_longer(cols=c("resistance","resilience"),names_to = "resil")|> 
  group_by(resil) |>
  mutate(value=scale(value)) |>
  ungroup() |>
  ggplot(aes(recovery,value,color=resil))+
  geom_point()+
  facet_wrap(~ncomb)


dist_test.lm=dist_test |> group_by(species_combination) |> mutate(n=n()) |> filter(n>2)

summary(lm(resilience~ncomb+ncomb:species_combination+wai+sgdd,data=dist_test.lm))
anova(lm(resilience~ncomb+ncomb:species_combination+wai+sgdd,data=dist_test.lm))
car::Anova(lm(resilience~ncomb+ncomb:species_combination+wai+sgdd,data=dist_test.lm))

