{library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(lme4)
library(glmmTMB)
library(effects)
library(DHARMa)

#### Load basic data ####
#%%%%%%%%%%%%%%%%%%%%%%%
tar_load(species.list.ipm)
tar_load(species.list.disturbance)

## get demo traits
tar_load(mean_demo)

pca=prcomp(mean_demo[,c("inv_50","ba_equil")],scale=TRUE,center=TRUE)
factoextra::fviz_pca_var(pca)

mean_demo$pca1<-pca$x[,1]


load("inv_30.RData")
# traits<-read.csv("data/traits_complete.csv") %>%
#   select(species,shade) %>%
#   mutate(species=gsub(" ","_",species)) %>%
#   rename(shade=shade) %>% 
#   left_join(mean_demo[,c("species","ba_equil","inv_50","pca1")] %>% 
#               mutate(species=gsub(" ","_",species))) %>% 
#   rename(inv_sp=inv_50)

traits<-read.csv("data/traits_complete.csv") %>%
  select(species,HM,WD) %>%
  mutate(species=gsub(" ","_",species)) %>%
  left_join(inv_30) |> 
  rename(inv_sp=inv_30)
GGally::ggpairs(traits[,c(2:4)])
traits |> ggplot(aes(shade,pca1))+geom_point()+geom_smooth()
rm(pca,mean_demo)

# compute mean climate var by categories
mean_pca<-tar_read(climate.cat)$FUNDIV_plotcat %>% 
  group_by(species,clim_id) %>% 
  summarise(pca1=mean(pca1),pca2=mean(pca2)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate(pca_sc=scale(pca1))

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
  
load("performance_2.RData")
}

performance |> left_join(traits |> mutate(species=gsub("_"," ",species))) |> 
  group_by(species,shade,inv_sp,ba_equil) |> filter(metric=="inv_50") |> summarise(inv_max=max(metric_val,na.rm=TRUE)) ->tt 
pca=prcomp(tt[,c("inv_max","inv_sp")],scale=TRUE,center=TRUE)
factoextra::fviz_pca_var(pca)
# ### Species exclusion ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# species.combination.excl$smallcombi<-NA
# for(i in 1:dim(species.combination.excl)[1]){
#   if(species.combination.excl$competexcluded[i]!=""){
#     s_p=gsub(" ","_",species.combination.excl$species[i])
#     sp_combi=unlist(strsplit(species.combination.excl$species_combination[i], split = "\\."))
#     # sp_combi=sp_combi[sp_combi!=s_p]
#     sp_ex=unlist(strsplit(species.combination.excl$competexcluded[i],"\\."))[-1]
# 
#     new_combi=paste(sort(sp_combi[!sp_combi%in%sp_ex]),collapse=".")
#     ex<-species.combination.excl %>%
#       filter(species==species.combination.excl$species[i]&
#                clim_id==species.combination.excl$clim_id[i]&
#                species_combination==new_combi)
#     if(dim(ex)[1]==1){
#       species.combination.excl$smallcombi[i]="present"
#     }else{
#       species.combination.excl$smallcombi[i]="absent"
#     }
#   }
# 
# }
# # for checking the proportion of combination represented after species exclusion
# # species.combination.excl %>%
# #   filter(competexcluded!="") %>%
# #   filter(smallcombi=="absent")
# 
# ### Performance metrics ####
# #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# # get disturbance metrics from mean
# disturbance_metric<-tar_read(disturbance_ba) %>%
#   mutate(elast="Abies_alba-amean",vr="mean") %>%
#   dplyr::select(species,elast,species_combination,clim_id,vr,resilience,resistance,recovery,matches(species.list.ipm))
# invasion_metric<-tar_read(invasion_ba) %>%
#   mutate(elast="Abies_alba-amean",vr="mean") %>% ## attention erruer
#   dplyr::select(species,elast,species_combination,clim_id,vr,inv_mean,inv_max,inv_50,BA_100,BA_500,BA_1000,matches(species.list.ipm))
# 
# 
# # compute mean elasticity per vital rates, for each species combi/clim
# disturbance <- disturbance_metric %>%
#   arrange(species,clim_id,species_combination,elast) %>%
#   pivot_longer(cols=matches(species.list.ipm)) %>%
#   left_join(traits, by = c("name" = "species")) %>%
#   group_by(species,elast,species_combination,clim_id,vr) %>%
#   mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
#                                TRUE~0)) %>%
#   pivot_longer(cols=colnames(traits)[-1],
#                values_to = "trait_val",
#                names_to = "trait_name") %>%
#   group_by(species,elast,species_combination,clim_id,vr,trait_name) %>%
#   mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
#          ba_target=sum(value*ba_multiple,na.rm = TRUE),
#          rel_ba_partner=ba_partner/(ba_target+ba_partner),
#          nih=sum((sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner,
#          nid=sum(abs(sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner) %>%
#   ungroup() %>%
#   dplyr::select(-c("ba_multiple","trait_val")) %>%
#   pivot_wider(values_from = c("nih","nid"),
#               names_from = "trait_name") %>%
#   pivot_wider(names_from = name,
#               values_from = value) %>%
#   dplyr::select(!matches(species.list.ipm))
# invasion <- invasion_metric %>% ungroup() %>%
#   arrange(species,species_combination,clim_id) %>%
#   pivot_longer(cols=matches(species.list.ipm)) %>%
#   left_join(traits, by = c("name" = "species")) %>%
#   group_by(species,elast,species_combination,clim_id,vr) %>%
#   mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
#                                TRUE~0),
#          value=case_when(!is.na(value)~1,
#                          gsub("_"," ",name)==species~1,
#                          TRUE~value)) %>%  # *1 because invasion are made at constant basal area of compet
#   pivot_longer(cols=colnames(traits)[-1],
#                values_to = "trait_val",
#                names_to = "trait_name") %>%
#   group_by(species,elast,species_combination,clim_id,vr,trait_name) %>%
#   mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
#          ba_target=sum(value*ba_multiple,na.rm = TRUE),
#          rel_ba_partner=ba_partner/(ba_target+ba_partner),
#          nih=sum((sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner,
#          nid=sum(abs(sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner) %>%
#   ungroup() %>%
#   dplyr::select(-c("ba_multiple","trait_val")) %>%
#   pivot_wider(values_from = c("nih","nid"),
#               names_from = "trait_name") %>%
#   pivot_wider(names_from = name,
#               values_from = value) %>%
#   dplyr::select(!matches(species.list.ipm))%>%
#   left_join(species.combination.excl,by=c("species","clim_id","species_combination"))
# ba_dif<-tar_read(ba_dif) %>%
#   select(!matches(c('ba_target',"ba_partner","n_species"))) %>%
#   mutate(elast="Abies_alba-amean",vr="mean")
# 
# 
# # gather all data
# performance <-invasion %>%
#   left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>%
#   left_join(ba_dif,by=c("species","elast","species_combination","clim_id","vr"))%>%
#   arrange(species,clim_id,species_combination,elast,vr) %>%
#   pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max",
#                       "inv_50","BA_100","BA_500","BA_1000","ba_dif"),
#                names_to="metric",
#                values_to="metric_val") %>%
#   arrange(species,clim_id,elast,species_combination)%>%
#   pivot_longer(cols=matches("\\.x|\\.y")) %>%
#   mutate(name_upd=gsub("\\.x$|\\.y$", "", name),
#          name=str_sub(name,-1)) %>%
#   pivot_wider(names_from = name,
#               values_from = value) %>%
#   mutate(compet_val=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~x,
#                               TRUE~y)) %>%
#   select(-c("x","y")) %>%
#   pivot_wider(names_from = name_upd,
#               values_from = compet_val)
#   # filter(excluded!="excluded") %>%
#   # filter(is.na(smallcombi)|smallcombi!="present")
#   save(performance,file = "performance_2.RData")


#### Plots ####
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
  # filter(species%in% c("Abies alba","Carpinus betulus","Fagus sylvatica","Fraxinus excelsior",
  #                      "Picea abies","Quercus ilex","Quercus petraea","Quercus robur")) %>%
  filter(vr=="mean") %>%
  filter(is.na(smallcombi)|smallcombi=="absent") %>% 
  ggplot() +
  geom_line(aes(pca1,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca1,metric_val,color=ba_partner, group=interaction(elast,species_combination)),size=1)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15))+
  theme_bw()+
  facet_wrap(metric~species,scales="free",ncol=4)+
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




### Coexistence ####
#%%%%%%%%%%%%%%%%%%%

# figure that helps

# analysis of percentage of species exclusion
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
         nih_mean=mean(nih_shade,na.rm=TRUE)) %>% 
  dplyr::select(species,clim_id,simul_state,n_comb,prop,nih_mean,shade) %>% unique() %>% 
  ungroup() 
data_exclu%>% 
  # filter(simul_state=="CompetitorExclusion") %>% 
  ggplot(aes(clim_id,prop,fill=simul_state))+
  geom_col()+
  facet_wrap(~species)

model1<-glmmTMB(prop~(1|species)+clim_id*nih_mean,
               data=subset(data_exclu,simul_state=="CompetitorExclusion"),
               family = beta_family(link = "logit"))

summary(model1)

# Simulate residuals
sim_res <- simulateResiduals(model1)
plot(sim_res)
ae1<-allEffects(model1)
as.data.frame(ae1$`clim_id:nih_mean`) |> 
  filter(nih_mean%in%c(-3,-0.5,2)) %>% 
  ggplot()+
  geom_ribbon(aes(x=clim_id,ymin=lower,ymax=upper,fill=as.factor(nih_mean)),alpha=0.2)+
  geom_line(aes(clim_id,fit,color=as.factor(nih_mean)))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  labs(x="Relative niche of species (cold/humid -> hot/dry)",
       y="BA at equilibrium / BA of pure stand",
       color="Competitive hierarchy",
       fill="Competitive hierarchy")+
  theme_classic()


rm(model1)
rm(data_exclu)

# #### Maintainance ####
#%%%%%%%%%%%%%%%%%%%%

# species alone ###
## data
data_maint_sp<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  # rename(pca_sc=pca1) |> 
  filter(species==gsub("_"," ",species_combination)) |> 
  filter(metric=="ba_dif") |> 
  left_join(traits %>% mutate(species=gsub("_"," ",species)))  

## plot
data_maint_sp |> 
  ggplot(aes(pca_sc,ba_target))+
  geom_line(aes(group=species))+
  geom_point(aes(color=species))

## model
trait=colnames(traits)[-1]

model2.0<-lmer(ba_target~(1|species)+pca_sc,
               data=data_maint_sp)
model2.1<-lmer(ba_target~(pca_sc|species)+pca_sc,
               data=data_maint_sp)
model2.2<-lmer(ba_target~(1|species)+pca_sc+I(pca_sc^2),
               data=data_maint_sp)
model2.3<-lmer(ba_target~(pca_sc|species)+pca_sc+I(pca_sc^2),
               data=data_maint_sp)
model2.4<-lmer(ba_target~(I(pca_sc^2)|species)+pca_sc+I(pca_sc^2),
               data=data_maint_sp)
## Calculate AIC and BIC
AIC_values <- AIC(model2.0,model2.1, model2.2,model2.3,model2.4)
BIC_values <- BIC(model2.0,model2.1, model2.2,model2.3,model2.4)

## Combine AIC and BIC into a data frame
model_comparison <- data.frame(
  Model = c("Model 2.0","Model 2.1", "Model 2.2", "Model 2.3", 
            "Model 2.4"),
  AIC = AIC_values$AIC,
  BIC = BIC_values$BIC
)

for(t in trait){
  print(t)
  data_maint_sp$trait<-data_maint_sp[[t]]
  
  model2.5<-lmer(ba_target~(pca_sc|species)+trait+pca_sc+I(pca_sc^2),
                 data=data_maint_sp)
  model2.6<-lmer(ba_target~(pca_sc|species)+pca_sc*trait+I(pca_sc^2),
                 data=data_maint_sp)
  model2.7<-lmer(ba_target~(pca_sc|species)+pca_sc+I(pca_sc^2)*trait,
                 data=data_maint_sp)
  model2.8<-lmer(ba_target~(pca_sc|species)+pca_sc*trait+I(pca_sc^2)*trait,
                 data=data_maint_sp)
  
  AIC_values <- AIC(model2.5,model2.6,model2.7,model2.8)
  BIC_values <- BIC(model2.5,model2.6,model2.7,model2.8)
  model_comp_t<- data.frame(
    Model = paste0(c( "Model 2.5","Model 2.6","Model 2.7","Model 2.8"),"_",t),
    AIC = AIC_values$AIC,
    BIC = BIC_values$BIC
  )
  model_comparison<-rbind(model_comparison,
                          model_comp_t)
  rm(model_comp_t)
}


model_comparison |> 
  arrange(AIC,BIC)


# Examine random effects variance
print(VarCorr(model2.1))
print(VarCorr(model2.2))
print(VarCorr(model2.7))

model2=lmer(ba_target~(pca_sc|species)+pca_sc+I(pca_sc^2)*pca1,
            data=data_maint_sp)
rm(model2.0,model2.1, model2.2,model2.3,model2.4,model2.5,model2.6,model2.7,model2.8,
   AIC_values,BIC_values,model_comparison)

## residuals
sim_res <- simulateResiduals(model2)
plot(sim_res)

summary(model2)
ae2<-allEffects(model2)
plot(ae2)
plot_2<-as.data.frame(ae2$`I(pca_sc^2):pca1`) |> 
  filter(pca1%in%c(-2,0.7,3)) |>
  mutate(`Competitive ability`=factor(case_when(pca1==-2~"Low",
                                              pca1==0.7~"Medium",
                                              pca1==3~"High"),
                                    levels=c("Low","Medium","High"))) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper,fill=`Competitive ability`),alpha=0.2)+
  geom_line(aes(pca_sc,fit,color=`Competitive ability`))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Basal area at equilibrium \n (m2/Ha)")
plot_2
ggsave(plot=plot_2,
       filename= "figure/sfe/Ba_quil_sp.png",
       dpi=600,
       width=10,
       height = 6,
       units = "cm")
rm(plot_2,model2,sim_res,ae2,data_maint_sp)       
  

# without traits
model2=lmer(ba_target~(pca_sc|species)+pca_sc+I(pca_sc^2),
            data=data_maint_sp)
## residuals
sim_res <- simulateResiduals(model2)
plot(sim_res)

summary(model2)
ae2<-allEffects(model2)
plot(ae2)

# Create a sequence of 'pca_sc' values
pca_seq <- seq(min(data_maint_sp$pca_sc), max(data_maint_sp$pca_sc), length.out = 100)
newdata <- data.frame(pca_sc = pca_seq) 
predictions <- predict(model2, newdata = newdata,re.form=~0, se.fit = TRUE)
pred_df <- cbind(newdata, fit = predictions$fit, se.fit = predictions$se.fit)

# Calculate confidence intervals
pred_df$lower <- pred_df$fit - 1.96 * pred_df$se.fit
pred_df$upper <- pred_df$fit + 1.96 * pred_df$se.fit

plot_222<-pred_df |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper),alpha=0.2)+
  geom_line(aes(pca_sc,fit))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Abundance in late successionnal stage \n (m2/Ha)")


plot_22<-as.data.frame(ae2$`I(pca_sc^2)`) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper),alpha=0.2)+
  geom_line(aes(pca_sc,fit))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Abundance in late successionnal stage \n (m2/Ha)")
plot_22
rm(plot_2,model2,sim_res,ae2,data_maint_sp)     

# Plot effect for shade level
# shade_levels <- c(1, 3, 5)
# 
# effects_list <- list()
# 
# for (i in seq_along(shade_levels)) {
#   shade_value <- shade_levels[i]
#   
#   effect_i <- Effect(
#     focal.predictors = "pca_sc",
#     mod = model2,
#     xlevels =100,
#     fixed.predictors = list(given.values = c(shade=shade_value))
#   )
#   
#   effects_list[[i]] <- as.data.frame(effect_i)
#   effects_list[[i]]$shade <- shade_value
# }
# 
# # Combine all effects into one data frame
# effects_df <- do.call(rbind, effects_list)
# 
# ggplot(effects_df, aes(x = pca_sc, y = fit, color = factor(shade))) +
#   geom_line(size = 1) +
#   geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(shade)), alpha = 0.2, color = NA) +
#   theme_minimal()

# with species combination
## data
data_maint<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>%
  filter(metric=="ba_dif") %>% 
  filter(!is.nan(nih_pca1)) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  mutate(species=forcats::fct_reorder(species, shade),
         simul_state=case_when(!is.na(smallcombi)~"CompetitorExclusion",
                               excluded=="excluded"~"SpeciesExclusion",
                               TRUE~"SpeciesCoex"),
         metric_val=case_when(metric_val>1~1,
                              TRUE~metric_val),
         metric_val=(metric_val * (dim(.)[1] - 1) + 0.5) / dim(.)[1]) 

## plots
data_maint %>% 
  filter(!is.nan(nih_pca1)) %>% 
  # filter(clim_id%in%c(1,5,10)) %>%
  # filter(species=="Fagus sylvatica") %>% View()
  ggplot(aes(pca_sc,metric_val,color=nih_pca1))+
  geom_point()+
  geom_smooth(method="lm",se=FALSE)+
  facet_wrap(~species)
data_maint %>% 
  filter(!is.nan(nih_pca1)) %>% 
  ggplot(aes(pca_sc,nih_pca1))+
  geom_point()+
  geom_smooth(method="gam",se=FALSE)

## model selections
model3.0<- glmmTMB(metric_val ~ (1|species) + pca_sc,
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.1<- glmmTMB(metric_val ~ (pca_sc|species)+pca_sc,
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.2<- glmmTMB(metric_val ~ (1|species)+pca_sc+I(pca_sc^2),
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.3<- glmmTMB(metric_val ~ (pca_sc|species)+pca_sc+I(pca_sc^2),
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.4<- glmmTMB(metric_val ~ (I(pca_sc^2)|species)+pca_sc+I(pca_sc^2),
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.5<- glmmTMB(metric_val ~ (pca_sc|species)+nih_pca1+pca_sc+I(pca_sc^2),
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.6<- glmmTMB(metric_val ~ (pca_sc|species)+pca_sc*nih_pca1+I(pca_sc^2),
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.7<- glmmTMB(metric_val ~ (pca_sc|species)+pca_sc+I(pca_sc^2)*nih_pca1,
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.8<- glmmTMB(metric_val ~ (pca_sc|species)+pca_sc*nih_pca1+I(pca_sc^2)*nih_pca1,
                   data = data_maint,
                   family = beta_family(link = "logit"))
model3.9<- glmmTMB(metric_val ~ pca_sc * nih_pca1 + I(pca_sc^2) * nih_pca1 + (nih_pca1 | species),
                   data = data_maint,
                   family = beta_family(link = "logit"))



# Calculate AIC and BIC
AIC_values <- AIC(model3.0,model3.1, model3.2,model3.3,model3.4,model3.5,model3.6,model3.7,model3.8,model3.9)
BIC_values <- BIC(model3.0,model3.1, model3.2,model3.3,model3.4,model3.5,model3.6,model3.7,model3.8,model3.9)

# Combine AIC and BIC into a data frame
model_comparison <- data.frame(
  Model = c("Model 3.0","Model 3.1", "Model 3.3", "Model 3.3","Model 3.4", 
            "Model 3.5","Model 3.6","Model 3.7","Model 3.8","Model 3.9"),
  AIC = AIC_values$AIC,
  BIC = BIC_values$BIC
)

model_comparison |> 
  arrange(AIC,BIC)

model3=model3.9

rm(model_comparison,model3.0,model3.1, model3.2,model3.3,model3.4,model3.5,model3.6,
   model3.7,model3.8,model3.9,AIC_values,BIC_values)
# Simulate residuals
sim_res <- simulateResiduals(model3)
plot(sim_res)


# draw effects
summary(model)
ae3<-allEffects(model3)
plot_3<-as.data.frame(ae3$`pca_sc:nih_pca1`) |> 
  filter(nih_pca1%in%c(-5,-0.1,4)) %>% 
  mutate(nih_pca1=factor(case_when(nih_pca1==-5~"High competition",
                                   nih_pca1==-0.1~"Medium competition",
                                   nih_pca1==4~"Low competition"),
                          levels=c("High competition","Medium competition","Low competition"))) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper,fill=as.factor(nih_pca1)),alpha=0.2)+
  geom_line(aes(pca_sc,fit,color=as.factor(nih_pca1)))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.25,0.5,0.75,1))+
  geom_hline(yintercept = 1,linetype="dashed")+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n (hot/dry -> cold/humid)",
       y="BA at equilibrium / BA of pure stand",
       color="Target species\n submitted to : ",
       fill="Target species\n submitted to : ")

ggsave(plot=plot_3,
       filename= "figure/sfe/Ba_quil_allsp-2.png",
       dpi=600,
       width=13,
       height = 8,
       units = "cm")
rm(plot_3,model3,sim_res,ae3,data_maint)       


# effect of species niche position


#### Resilience ####
#%%%%%%%%%%%%%%%%%%%

# species alone ###
## data
data_resilience_sp<- performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(!is.na(metric_val)) |> 
  filter(metric=="resilience") %>%
  filter(vr=="mean") %>%
  filter(species==gsub("_"," ",species_combination)) |> 
  mutate(metric_val=ba_target*metric_val,
         metric_val=log(metric_val)) 
## plots
data_resilience_sp |> 
  group_by(species) |> 
  # mutate(metric_val=scale(metric_val,scale = TRUE,center=FALSE)) |>
  ggplot(aes(pca_sc,metric_val))+
  geom_point(aes(group=species,color=pca1))+
  geom_smooth()+
  facet_wrap(~species)
data_resilience_sp |> 
  group_by(species) |> 
  ggplot(aes(pca_sc,metric_val))+
  geom_point(aes(group=species,color=pca1))+
  geom_line(aes(group=species,color=pca1))

data_resilience_sp |> 
  group_by(species) |> 
  mutate(metric_val=scale(metric_val,scale = TRUE,center=FALSE)) |>
  ggplot(aes(metric_val))+
  geom_density()
summary(data_resilience_sp)


## models
trait=colnames(traits)[-1]
model4.0<-glmmTMB(metric_val~(1|species)+pca_sc,
                  data=data_resilience_sp)
model4.1<-glmmTMB(metric_val~(pca_sc|species)+pca_sc,
                  data=data_resilience_sp)
model4.2<-glmmTMB(metric_val~(1|species)+pca_sc+I(pca_sc^2),
                  data=data_resilience_sp)
model4.3<-glmmTMB(metric_val~(pca_sc|species)+pca_sc+I(pca_sc^2),
                  data=data_resilience_sp)
model4.4<-glmmTMB(metric_val~(I(pca_sc^2)|species)+pca_sc+I(pca_sc^2),
                  data=data_resilience_sp)



## Calculate AIC and BIC
AIC_values <- AIC(model4.0,model4.1, model4.2,model4.3,model4.4)
BIC_values <- BIC(model4.0,model4.1, model4.2,model4.3,model4.4)

## Combine AIC and BIC into a data frame
model_comparison <- data.frame(
  Model = c("Model 4.0","Model 4.1", "Model 4.2", "Model 4.3", 
            "Model 4.4"),
  AIC = AIC_values$AIC,
  BIC = BIC_values$BIC
)



for(t in trait){
  print(t)
  data_resilience_sp$trait<-data_resilience_sp[[t]]
  
  model4.5<-glmmTMB(metric_val~(pca_sc|species)+trait+pca_sc+I(pca_sc^2),
                    data=data_resilience_sp)
  model4.6<-glmmTMB(metric_val~(pca_sc|species)+pca_sc*trait+I(pca_sc^2),
                    data=data_resilience_sp)
  model4.7<-glmmTMB(metric_val~(pca_sc|species)+pca_sc+I(pca_sc^2)*trait,
                    data=data_resilience_sp)
  model4.8<-glmmTMB(metric_val~(pca_sc|species)+pca_sc*trait+I(pca_sc^2)*trait,
                    data=data_resilience_sp)
  AIC_values <- AIC(model4.5,model4.6,model4.7,model4.8)
  BIC_values <- BIC(model4.5,model4.6,model4.7,model4.8)
  model_comp_t<- data.frame(
    Model = paste0(c("Model 4.5","Model 4.6","Model 4.7","Model 4.8"),"_",t),
    AIC = AIC_values$AIC,
    BIC = BIC_values$BIC
  )
  model_comparison<-rbind(model_comparison,
                          model_comp_t)
  rm(model_comp_t)
}


model_comparison |> 
  arrange(AIC,BIC)


model4<-glmmTMB(metric_val~(pca_sc|species)+pca_sc+I(pca_sc^2)*inv_sp,
                data=data_resilience_sp)
rm(model4.0,model4.1, model4.2,model4.3,model4.4,model4.5,model4.6,model4.7,model4.8,
   model_comparison, AIC_values,BIC_values)


## Simulate residuals
sim_res <- simulateResiduals(model4)
plot(sim_res)


## draw effects
summary(model4)
ae4<-allEffects(model4)
plot_4<-as.data.frame(ae4$`I(pca_sc^2):inv_sp`) |>
  filter(inv_sp %in% c(0.05,0.3,0.5)) |> 
  mutate(fit=exp(fit),
         lower=exp(lower),
         upper=exp(upper)) |> 
  mutate(`Competive ability`=factor(case_when(inv_sp==0.05~"Low",
                                              inv_sp==0.3~"Medium",
                                              inv_sp==0.5~"High"),
                                     levels=c("Low","Medium","High"))) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper,fill=`Competive ability`),alpha=0.2)+
  geom_line(aes(pca_sc,fit,color=`Competive ability`))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Resilience (1/m2)")

ggsave(plot=plot_4,
       filename= "figure/sfe/resilience_sp_2.png",
       dpi=600,
       width=13,
       height = 8,
       units = "cm")
rm(plot_4,model4,sim_res,ae4,data_resilience_sp)   

## alone
model4=glmmTMB(metric_val~(pca_sc|species)+pca_sc,
               data=data_resilience_sp)
ae4<-allEffects(model4) 
plot_44<-as.data.frame(ae4$pca_sc) |> 
  mutate(fit=exp(fit),
         lower=exp(lower),
         upper=exp(upper)) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper),alpha=0.2)+
  geom_line(aes(pca_sc,fit))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Resilience (1/m2)")
plot_44

# with species combination ###

## data

data_resilience<- performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(!is.na(metric_val)) |> 
  filter(metric=="resilience") %>%
  filter(vr=="mean") %>%
  filter(is.na(smallcombi)|smallcombi=="absent") |> 
  rowwise() %>% 
  mutate(n_species=length(strsplit(species_combination,"\\.")[[1]])) %>% 
  ungroup() %>% 
  group_by(species,clim_id) %>% 
  arrange(species,clim_id,n_species) %>% 
  mutate(metric_val=ba_target*metric_val,
         res_dif=log(metric_val/metric_val[1])) |> 
  filter(!species%in%c("Pinus pinaster","Pinus sylvestris","Pinus uncinata")) |> 
  filter(species!=gsub("_"," ",species_combination)) 
  

## plots
data_resilience%>% 
  mutate(ba_partner=case_when(ba_partner==0~NA,
                              TRUE~ba_partner)) %>% 
  # filter(species%in% c("Abies alba","Carpinus betulus","Fagus sylvatica","Fraxinus excelsior",
  #                      "Picea abies","Quercus ilex","Quercus petraea","Quercus robur")) %>%
  ggplot() +
  geom_line(aes(pca_sc,res_dif, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca_sc,res_dif,color=nih_inv_sp, group=interaction(elast,species_combination)),size=3)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15))+
  theme_bw()+
  facet_wrap(metric~species,ncol=4,scales = "free")+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Total basal area of competitors")

## models
model5.0<-glmmTMB(res_dif~(1|species)+pca_sc,
                  data=data_resilience)
model5.1<-glmmTMB(res_dif~(pca_sc|species)+pca_sc,
                  data=data_resilience)
model5.2<-glmmTMB(res_dif~(1|species)+pca_sc+I(pca_sc^2),
                  data=data_resilience)
model5.3<-glmmTMB(res_dif~(pca_sc|species)+pca_sc+I(pca_sc^2),
                  data=data_resilience)
model5.4<-glmmTMB(res_dif~(I(pca_sc^2)|species)+pca_sc+I(pca_sc^2),
                  data=data_resilience)
model5.5<-glmmTMB(res_dif~(pca_sc|species)+nih_inv_sp+pca_sc+I(pca_sc^2),
                  data=data_resilience)
model5.6<-glmmTMB(res_dif~(pca_sc|species)+pca_sc*nih_inv_sp+I(pca_sc^2),
                  data=data_resilience)
model5.7<-glmmTMB(res_dif~(pca_sc|species)+pca_sc+I(pca_sc^2)*nih_inv_sp,
                  data=data_resilience)
model5.8<-glmmTMB(res_dif~(pca_sc|species)+pca_sc*nih_inv_sp+I(pca_sc^2)*nih_inv_sp,
                  data=data_resilience)
model5.9<-glmmTMB(res_dif~(nih_inv_sp|species)+pca_sc*nih_inv_sp+I(pca_sc^2),
                  data=data_resilience)
model5.10<-glmmTMB(res_dif~(nih_inv_sp|species)+pca_sc+I(pca_sc^2)*nih_inv_sp,
                  data=data_resilience)
model5.11<-glmmTMB(res_dif~(nih_inv_sp|species)+pca_sc*nih_inv_sp+I(pca_sc^2)*nih_inv_sp,
                  data=data_resilience)



## Calculate AIC and BIC
AIC_values <- AIC(model5.0,model5.1, model5.2,model5.3,model5.4,model5.5,model5.6,
                  model5.7,model5.8,model5.9,model5.10,model5.11)
BIC_values <- BIC(model5.0,model5.1, model5.2,model5.3,model5.4,model5.5,model5.6,
                  model5.7,model5.8,model5.9,model5.10,model5.11)

## Combine AIC and BIC into a data frame
model_comparison <- data.frame(
  Model = c("Model 5.0","Model 5.1", "Model 5.2", "Model 5.3","Model 5.4","Model 5.5",
            "Model 5.6","Model 5.7","Model 5.8","Model 5.9","Model 5.10","Model 5.11"),
  AIC = AIC_values$AIC,
  BIC = BIC_values$BIC
)

model_comparison |> 
  arrange(AIC,BIC)


model5<-model5.6
rm(model5.0,model5.1, model5.2,model5.3,model5.4,model5.5,model5.6,model5.7,model5.8,
   model5.9,model5.10,model5.11,model_comparison, AIC_values,BIC_values)

## Simulate residuals
sim_res <- simulateResiduals(model5)
plot(sim_res)


## draw effects
summary(model5)
ae5<-allEffects(model5)
plot(ae5)
plot_5<-as.data.frame(ae5$`pca_sc:nih_inv_sp`) |> 
  mutate(fit=exp(fit),
         lower=exp(lower),
         upper=exp(upper)) |> 
  filter(nih_inv_sp%in%c(-0.4,0.04,0.5)) |>
  mutate(nih_inv_sp=factor(case_when(nih_inv_sp==-0.4~"High competition",
                                     nih_inv_sp==0.04~"Medium competition",
                                     nih_inv_sp==0.5~"Low competition"),
                          levels=c("High competition","Medium competition","Low competition"))) |>
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper,fill=nih_inv_sp),alpha=0.2)+
  geom_line(aes(pca_sc,fit,color=nih_inv_sp))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  theme_classic()+
  theme(text = element_text(size=10))+
  scale_y_continuous(limits = c(0, 3.1), breaks = c(0,0.25,0.5,0.75,1,1.5,2,3))+
  geom_hline(yintercept = 1,linetype="dashed")+
  labs(x="Climatic range scaled by species \n (hot/dry -> cold/humid)",
       y="Resilience in competition / pure stand ",
       color="Target species \n submitted to:",
       fill="Target species \n submitted to:")

ggsave(plot=plot_5,
       filename= "figure/sfe/resilience_allsp_2.png",
       dpi=600,
       width=13,
       height = 8,
       units = "cm")
rm(plot_5,model5,sim_res,ae5,data_resilience)   

# model2 <- lmer(metric_val ~(1|species) + pca_sc * ba_partner,
#               data = data_resilience)


#### Invasion ####
#%%%%%%%%%%%%%%%%%%%

# species alone ###
## data
data_inv_sp<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(metric=="inv_50") %>% 
  filter(species==gsub("_"," ",species_combination))


## plot
ggplot(data_inv_sp) +
  geom_line(aes(pca_sc,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca_sc,metric_val,color=pca1, group=interaction(elast,species_combination)),size=1)+
  theme_bw()+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Invasion rate in 50 first years")

## model
trait=colnames(traits)[-1]

model6.0<-glmmTMB(metric_val~(1|species)+pca_sc,
                  data=data_inv_sp)
model6.1<-glmmTMB(metric_val~(pca_sc|species)+pca_sc,
                  data=data_inv_sp)
model6.2<-glmmTMB(metric_val~(1|species)+pca_sc+I(pca_sc^2),
                  data=data_inv_sp)
model6.3<-glmmTMB(metric_val~(pca_sc|species)+pca_sc+I(pca_sc^2),
                  data=data_inv_sp)
model6.4<-glmmTMB(metric_val~(I(pca_sc^2)|species)+pca_sc+I(pca_sc^2),
                  data=data_inv_sp)

## Calculate AIC and BIC
AIC_values <- AIC(model6.0,model6.1, model6.2,model6.3,model6.4)
BIC_values <- BIC(model6.0,model6.1, model6.2,model6.3,model6.4)

## Combine AIC and BIC into a data frame
model_comparison <- data.frame(
  Model = c("Model 6.0","Model 6.1", "Model 6.2", "Model 6.3", 
            "Model 6.4"),
  AIC = AIC_values$AIC,
  BIC = BIC_values$BIC
)

for(t in trait){
  print(t)
  data_inv_sp$trait<-data_inv_sp[[t]]

  model6.5<-glmmTMB(metric_val~(pca_sc|species)+trait+pca_sc+I(pca_sc^2),
                    data=data_inv_sp)
  model6.6<-glmmTMB(metric_val~(pca_sc|species)+pca_sc*trait+I(pca_sc^2),
                    data=data_inv_sp)
  model6.7<-glmmTMB(metric_val~(pca_sc|species)+pca_sc+I(pca_sc^2)*trait,
                    data=data_inv_sp)
  model6.8<-glmmTMB(metric_val~(pca_sc|species)+pca_sc*trait+I(pca_sc^2)*trait,
                    data=data_inv_sp)
  
  AIC_values <- AIC(model6.5,model6.6,model6.7,model6.8)
  BIC_values <- BIC(model6.5,model6.6,model6.7,model6.8)
  model_comp_t<- data.frame(
    Model = paste0(c("Model 6.5","Model 6.6","Model 6.7","Model 6.8"),"_",t),
    AIC = AIC_values$AIC,
    BIC = BIC_values$BIC
  )
  model_comparison<-rbind(model_comparison,
                          model_comp_t)
  rm(model_comp_t)
}



model_comparison |> 
  arrange(AIC,BIC)

model6<-glmmTMB(metric_val~(pca_sc|species)+pca_sc+I(pca_sc^2)*pca1,
                data=data_inv_sp)
rm(model6.0,model6.1, model6.2,model6.3,model6.4,model6.5,model6.6,model6.7,model6.8,
   model_comparison, AIC_values,BIC_values)

## Simulate residuals
sim_res <- simulateResiduals(model6)
plot(sim_res)


## draw effects
summary(model6)
ae6<-allEffects(model6)
plot_6<-as.data.frame(ae6$`I(pca_sc^2):pca1`) |> 
  filter(pca1 %in% c(-2,0.7,3)) |> 
  mutate(`Competive ability`=factor(case_when(pca1==-2~"Low",
                                              pca1==0.7~"Medium",
                                              pca1==3~"High"),
                                    levels=c("Low","Medium","High"))) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper,fill=`Competive ability`),alpha=0.2)+
  geom_line(aes(pca_sc,fit,color=`Competive ability`))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Invasion rate \n (m2/year)")


ggsave(plot=plot_6,
       filename= "figure/sfe/invasion_sp.png",
       dpi=600,
       width=13,
       height = 8,
       units = "cm")
rm(plot_6,model6,sim_res,ae6,data_inv_sp)   

# without traits
model6=model6.4
ae6<-allEffects(model6)
plot_66<-as.data.frame(ae6$`I(pca_sc^2)`) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper),alpha=0.2)+
  geom_line(aes(pca_sc,fit))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  theme_classic()+
  theme(text = element_text(size=10))+
  labs(x="Climatic range scaled by species \n (hot/dry -> cold/humid)",
       y="Growth in early successionnal stage \n (m2/year)")


# model3 <- lmer(inv_dif ~(pca_sc|species) + pca_sc * nih + I(pca_sc^2) ,
#                data = data_inv)


# with combinations
## data
data_inv<- performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  filter(metric=="inv_50") %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  rowwise() %>% 
  mutate(n_species=length(strsplit(species_combination,"\\.")[[1]])) %>% 
  ungroup() %>% 
  group_by(species,clim_id) %>% 
  arrange(species,clim_id,n_species) %>% 
  mutate(inv_dif=metric_val/metric_val[1]) %>% 
  ungroup() %>% 
  mutate(nih=case_when(is.nan(nih_inv_sp)~0,
                       TRUE~nih_inv_sp),
         inv_dif=case_when(inv_dif<0~0,
                           TRUE~inv_dif),
         inv_dif=(inv_dif * (dim(.)[1] - 1) + 0.5) / dim(.)[1]) %>% 
  filter(!species%in%c("Pinus pinaster","Pinus pinea")) %>% 
  filter(species!=gsub("_"," ",species_combination))

## plots
ggplot(data_inv) +
  geom_line(aes(pca_sc,metric_val, group=interaction(elast,species_combination)),color="grey")+
  geom_point(aes(pca_sc,metric_val,color=nih, group=interaction(elast,species_combination)),size=1)+
  geom_hline(yintercept = 0)+
  scale_color_gradientn(colours = viridis(15))+
  theme_bw()+
  facet_wrap(metric~species,ncol=4)+
  labs(x="PCA first axis (cold/wet -> hot/dry)",
       y="Performance",
       color="Total basal area of competitors")


## models
model7.0 <- glmmTMB(inv_dif ~ (1|species) + pca_sc,
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.1 <- glmmTMB(inv_dif ~ (pca_sc|species) + pca_sc,
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.2 <- glmmTMB(inv_dif ~ (1|species) + pca_sc+I(pca_sc^2),
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.3 <- glmmTMB(inv_dif ~ (pca_sc|species) + pca_sc + I(pca_sc^2),
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.4 <- glmmTMB(inv_dif ~ (I(pca_sc^2)|species) + pca_sc + I(pca_sc^2),
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.5 <- glmmTMB(inv_dif ~ (pca_sc|species) + nih_pca1 + pca_sc + I(pca_sc^2),
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.6 <- glmmTMB(inv_dif ~ (pca_sc|species)+pca_sc*nih_pca1+I(pca_sc^2),
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.7 <- glmmTMB(inv_dif ~ (pca_sc|species)+pca_sc+I(pca_sc^2)*nih_pca1,
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.8 <- glmmTMB(inv_dif ~(pca_sc|species)+pca_sc*nih_pca1+I(pca_sc^2)*nih_pca1,
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.9 <- glmmTMB(inv_dif ~ (nih_pca1|species)+pca_sc*nih_pca1+I(pca_sc^2),
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.10 <- glmmTMB(inv_dif ~ (nih_pca1|species)+pca_sc+I(pca_sc^2)*nih_pca1,
                    data = data_inv,
                    family = beta_family(link = "logit"))
model7.11 <- glmmTMB(inv_dif ~ (nih_pca1|species)+pca_sc*nih_pca1+I(pca_sc^2)*nih_pca1,
                     data = data_inv,
                     family = beta_family(link = "logit"))

# model3 <- glmmTMB(inv_dif ~ pca_sc * ba_partner + I(pca_sc^2) * ba_partner+ (ba_partner | species),
#                  data = data_inv,
#                  family = beta_family(link = "logit"))

## Calculate AIC and BIC
AIC_values <- AIC(model7.0,model7.1, model7.2,model7.3,model7.4,model7.5,model7.6,
                  model7.7,model7.8,model7.9,model7.10,model7.11)
BIC_values <- BIC(model7.0,model7.1, model7.2,model7.3,model7.4,model7.5,model7.6,
                  model7.7,model7.8,model7.9,model7.10,model7.11)

## Combine AIC and BIC into a data frame
model_comparison <- data.frame(
  Model = c("Model 7.0","Model 7.1", "Model 7.2", "Model 7.3","Model 7.4","Model 7.5",
            "Model 7.6","Model 7.7","Model 7.8","Model 7.9","Model 7.10","Model 7.11"),
  AIC = AIC_values$AIC,
  BIC = BIC_values$BIC
)

model_comparison |> 
  arrange(AIC,BIC)


model7<-model7.9
rm(model7.0,model7.1, model7.2,model7.3,model7.4,model7.5,model7.6,model7.7,model7.8,
   model7.9,model7.10,model7.11,model_comparison, AIC_values,BIC_values)



# Simulate residuals
sim_res <- simulateResiduals(model7)
plot(sim_res)


ae7<-allEffects(model7)


plot_7 <- as.data.frame(ae7$`pca_sc:nih_pca1`) |>
  filter(nih_pca1%in%c(-5,-0.1,4)) |>
  mutate(nih_pca1=factor(case_when(nih_pca1==-5~"High competition",
                                   nih_pca1==-0.1~"Medium competition",
                                   nih_pca1==4~"Low competition"),
                         levels=c("High competition","Medium competition","Low competition"))) |>
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymin=lower,ymax=upper,fill=nih_pca1),alpha=0.2)+
  geom_line(aes(pca_sc,fit,color=nih_pca1))+
  scale_color_manual(values=c("burlywood1","tan1","tan4"))+
  scale_fill_manual(values=c("burlywood1","tan1","tan4"))+
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0,0.25,0.5,0.75,1))+
  geom_hline(yintercept = 1,linetype="dashed")+
  labs(x="Climatic range scaled by species \n ( hot/dry -> cold/humid)",
       y="Invasion rate in competition / pure stand",
       color="Target species \n submitted to:",
       fill="Target species \n submitted to:")+
  theme_classic()

ggsave(plot=plot_7,
       filename= "figure/sfe/invasion_allsp.png",
       dpi=600,
       width=13,
       height = 8,
       units = "cm")

rm(plot_7,model7,sim_res,ae7,data_inv)   


#### Common plots ####
#%%%%%%%%%%%%%%%%%%%

cowplot::plot_grid(plot_2+theme(legend.position = "none"),
                   plot_6+theme(legend.position = "none"),
                   plot_4+theme(legend.position = "none"),
                   ggpubr::get_legend(plot_2),
                   rel_widths = c(1,1,1,0.4),
                   ncol=4)->p
ggsave(plot=p,
       filename= "figure/sfe/pure_stand.png",
       dpi=600,
       width=30,
       height = 9,
       units = "cm")



cowplot::plot_grid(plot_3+theme(legend.position = "none"),
                   plot_7+theme(legend.position = "none"),
                   ggpubr::get_legend(plot_3),
                   rel_widths = c(1,1,0.3),
                   ncol=3)->pp
pp
ggsave(plot=pp,
       filename= "figure/sfe/inv_maintenance.png",
       dpi=600,
       width=30,
       height = 9,
       units = "cm")

cowplot::plot_grid(plot_22,
                   plot_66,
                   plot_44,
                   ncol=3)->ppp
ppp
ggsave(plot=ppp,
       filename= "figure/sfe/pure_stand_nocompet.png",
       dpi=600,
       width=30,
       height = 9,
       units = "cm")

#### competition ####
#%%%%%%%%%%%%%%%%%%%

# test effect of competition
data.compet<- performance %>% 
  left_join(species.combination.select,by=c("species","clim_id","species_combination")) %>% 
  filter(metric=="resilience") %>% 
  group_by(species) %>% 
  mutate(ba_tot=ba_target+ba_partner,
         ba_tot=scale(ba_tot,scale=TRUE)) 
summary(lmer(ba_tot~pca1+(pca1|species),data=data.compet))


