library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(lme4)

tar_load(models_traits)
models_traits$traits_effect |> 
  filter(!grepl("sp",data_name)) |> 
  ggplot()+
  geom_segment(aes(x=lwr,xend=upr,y=trait,yend=trait))+
  geom_point(aes(x=effect,y=trait))+
  theme_classic()+
  labs(y="",x="Standardized trait effect")+
  facet_grid(interaction~data_name,scales="free_x")




models_traits$predict_traits |>
  filter(grepl("sp",data_name)) |> 
  group_by(data_name,pred_name) |> 
  filter(is.na(pred)) |>
  mutate(data_name=case_when(data_name=="data_inv_sp"~"Invasion",
                             data_name=="data_maint_sp"~"Maintenance",
                             data_name=="data_resilience_sp"~"Resilience")) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymax=upr,ymin=lwr),
              alpha=0.2)+
  geom_line(aes(pca_sc,predicted_mean),
            size=1)+
  scale_color_manual(values=c("darkseagreen","darksalmon","darkred"))+
  scale_fill_manual(values=c("darkseagreen","darksalmon","darkred"))+
  facet_wrap(~data_name,scale="free")+
  theme_classic()+
  theme(panel.background = element_rect(fill="white",colour="grey50"))+
  labs(color="Strength of competition \nexperienced by target species :",
       fill="Strength of competition \nexperienced by target species :",
       y="",
       x="Climate scaled by species")


models_traits$predict_traits |> 
  filter(!grepl("sp",data_name)) |> 
  group_by(data_name,pred_name) |> 
  # filter(data_name=="data_inv",pred_name=="nih_HM") |> 
  mutate(data_name=case_when(data_name=="data_inv"~"Invasion",
                             data_name=="data_maint"~"Maintenance",
                             data_name=="data_resilience"~"Resilience"),
         pred_name=case_when(pred_name=="nih_HM"~"NIH(Max Heigth)",
                             pred_name=="nih_WD"~"NIH(Wood Density)",
                             pred_name=="nih_recruitment"~"NIH(Recruitment)",
                             pred_name=="none"~"none"),
         pred_cat=case_when(pred==min(pred)~"High compet",
                            pred==max(pred)~"Low compet",
                            TRUE~"Medium compet"),
         pred_cat=factor(pred_cat,levels=c("Low compet","Medium compet","High compet"))) |> 
  # filter(pred_cat%in%c("Low compet","High compet")) |> 
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
