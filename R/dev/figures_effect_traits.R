### dataframes

{data_maint_sp<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  # rename(pca_sc=pca1) |> 
  filter(species==gsub("_"," ",species_combination)) |> 
  filter(metric=="ba_dif") |> 
  left_join(traits %>% mutate(species=gsub("_"," ",species)))  

data_maint<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>%
  filter(metric=="ba_dif") %>% 
  filter(!is.nan(nih_HM)) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  mutate(species=forcats::fct_reorder(species, HM),
         simul_state=case_when(!is.na(smallcombi)~"CompetitorExclusion",
                               excluded=="excluded"~"SpeciesExclusion",
                               TRUE~"SpeciesCoex"),
         metric_val=case_when(metric_val>1~1,
                              TRUE~metric_val),
         metric_val=(metric_val * (dim(.)[1] - 1) + 0.5) / dim(.)[1]) 

data_resilience_sp<- performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(!is.na(metric_val)) |> 
  filter(metric=="resilience") %>%
  filter(vr=="mean") %>%
  filter(species==gsub("_"," ",species_combination)) |> 
  mutate(metric_val=ba_target*metric_val,
         metric_val=log(metric_val)) 

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

data_inv_sp<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(metric=="inv_50") %>% 
  filter(species==gsub("_"," ",species_combination))


data_inv<- performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  filter(metric=="inv_50") %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  rowwise() %>% 
  mutate(n_species=length(strsplit(species_combination,"\\.")[[1]])) %>% 
  ungroup() %>% 
  group_by(species,clim_id) %>% 
  arrange(species,clim_id,n_species) %>% 
  mutate(inv_dif=metric_val/metric_val[1],
         linv_dif=log(inv_dif)) %>% 
  ungroup() %>% 
  mutate(nih=case_when(is.nan(nih_inv_sp)~0,
                       TRUE~nih_inv_sp)#,
         # inv_dif=case_when(inv_dif<0~0,
         #                   TRUE~inv_dif),
         # inv_dif=(inv_dif * (dim(.)[1] - 1) + 0.5) / dim(.)[1]
  ) %>% 
  # filter(!species%in%c("Pinus pinaster","Pinus pinea")) %>% 
  filter(species!=gsub("_"," ",species_combination)) |> 
  filter(!is.infinite(linv_dif))

### prepare datafrun

run_mod<-data.frame(data_name=c("data_maint_sp","data_maint","data_resilience_sp",
                       "data_resilience","data_inv_sp","data_inv"),
           response_name=c("ba_target","metric_val","metric_val",
                           "res_dif","metric_val","inv_dif"),
           response_var=c("response","response","response","response","response","response"),
           fixed_predictor=c("pca_sc","pca_sc","pca_sc","pca_sc","pca_sc","pca_sc"),
           predictors_list=c("","nih_","","nih_","","nih_"),
           group_var=c("species","species","species","species","species","species"),
           mod_extension=c("",",family = beta_family(link = \"logit\")","",
                           "","",",family = Gamma(link = \"log\")"),
           mod.type=c("lmer","glmmTMB","glmmTMB","glmmTMB","glmmTMB","glmmTMB"))





predict_traits<-setNames(data.frame(matrix(nrow = 0,ncol = 10)),
                         nm=c('data_name','response','pca_sc','trait','trait_value','predicted_mean','predicted_se','lwr','upr','pred'))
class(predict_traits$data_name)<-"character"
traits_effect<-setNames(data.frame(matrix(nrow = 0,ncol = 6)),
                        nm=c('data_name',"interaction",'trait','effect',"lwr","upr"))
best_models<-setNames(data.frame(matrix(nrow = 0,ncol = 6)),
                      nm=c('data_name',"formulas",'trait',"model","AIC","ncof"))

}
### RUN



traits_effect |> 
  filter(!grepl("sp",data_name)) |> 
  ggplot()+
  geom_segment(aes(x=lwr,xend=upr,y=trait,yend=trait))+
  geom_point(aes(x=effect,y=trait))+
  theme_classic()+
  labs(y="",x="Standardized trait effect")+
  facet_grid(interaction~data_name,scales="free_x")


predict_traits |> 
  filter(!grepl("sp",data_name)) |> 
  group_by(data_name,pred_name) |> 
  mutate(data_name=case_when(data_name=="data_inv"~"Invasion",
                             data_name=="data_maint"~"Maintenance",
                             data_name=="data_resilience"~"Resilience"),
         pred_name=case_when(pred_name=="nih_HM"~"NIH(Max Heigth)",
                             pred_name=="nih_WD"~"NIH(Wood Density)",
                             pred_name=="nih_inv_sp"~"NIH(Recruitment)",
                             pred_name=="none"~"none"),
         pred_cat=case_when(pred==min(pred)~"High compet",
                            pred==max(pred)~"Low compet",
                            TRUE~"Medium compet"),
         pred_cat=factor(pred_cat,levels=c("Low compet","Medium compet","High compet"))) |> 
  filter(pred_cat%in%c("Low compet","High compet")) |>
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymax=upr,ymin=lwr,
                  fill=pred_cat),
              alpha=0.2)+
  geom_line(aes(pca_sc,predicted_mean,color=pred_cat),
            size=1)+
  geom_hline(yintercept = 1)+
  scale_color_manual(values=c("darkseagreen","darksalmon","darkred"))+
  scale_fill_manual(values=c("darkseagreen","darksalmon","darkred"))+
  facet_grid(data_name~pred_name)+
  theme_classic()+
  theme(panel.background = element_rect(fill="white",colour="grey50"))+
  labs(color="Strength of competition \nexperienced by target species :",
       fill="Strength of competition \nexperienced by target species :",
       y="",
       x="Climate scaled by species")
