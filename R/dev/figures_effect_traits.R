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


generate_formulas <- function(response_var, fixed_predictor, predictors_list,
                              group_var,mod_extension,mod.type) {
  formulas<-list()
  if(mod.type=="lmer"){
    mod.ran="lmer("
    mod.nran="lm("
  }else{
    mod.ran="glmmTMB("
    mod.nran="glmmTMB("
  }
  formulas[[1]] <- paste0(mod.ran,response_var,"~ (1|", group_var, ")",mod_extension,",data=data_mod)")
  for (pred in predictors_list) {
    formula_pred<-data.frame(model="model",
                             response=response_var,
                             clim=fixed_predictor) |> 
      crossing(quad_clim=c("",paste0("I(",fixed_predictor,"^2)"))) |> 
      crossing(predictor=c("",pred)) |>
      crossing(quad_pred=c("")) |> #,paste0("I(",pred,"^2)") 
      crossing(data.frame(random.eff=c("none","intercept","slope"),
                          formula.random=c("",paste0("(1|", group_var, ")"),paste0( "|", group_var, ")")))) |> 
      crossing(random.var=c("clim","quad_clim")) |> #,"predictor","quad_pred"
      crossing(interaction=c("none","clim","quad_clim","both")) |> 
      mutate(random.var=case_when(random.eff!="slope"~"",
                                  TRUE~random.var)) |> unique() |> 
      mutate(kp=case_when(random.var=="clim"&clim==""~FALSE,
                          random.var=="predictor"&predictor==""~FALSE,
                          random.var=="quad_clim"&quad_clim==""~FALSE,
                          random.var=="quad_pred"&quad_pred==""~FALSE,
                          random.var==""~TRUE,
                          TRUE~TRUE))|> 
      filter(kp) |> select(-kp) |> 
      mutate(interaction=case_when(predictor==""~NA,
                                   interaction=="both"&quad_clim==""~NA,
                                   interaction=="quad_clim"&quad_clim==""~NA,
                                   TRUE~interaction)) |> 
      unique() |> 
      mutate(interaction=case_when(interaction=="none"~"",
                                   interaction=="clim"~paste0(clim,"*",predictor),
                                   interaction=="quad_clim"~paste0(quad_clim,"*",predictor),
                                   interaction=="both"~paste0(clim,"*",predictor,"+",quad_clim,"*",predictor)),
             clim=case_when(random.eff=="slope"&random.var=="clim"~paste0(clim,"+(",clim,formula.random),
                            TRUE~clim),
             quad_clim=case_when(random.eff=="slope"&random.var=="quad_clim"~paste0(quad_clim,"+(",quad_clim,formula.random),
                                 TRUE~quad_clim),
             predictor=case_when(random.eff=="slope"&random.var=="predictor"~paste0(predictor,"+ (",predictor,formula.random),
                                 TRUE~predictor),
             quad_pred=case_when(random.eff=="slope"&random.var=="quad_pred"~paste0(quad_pred,"+(",quad_pred,formula.random),
                                 TRUE~quad_pred),
             rand=case_when(random.eff=="intercept"~formula.random,
                            TRUE~NA),
             across(c("clim","quad_clim","predictor","quad_pred","interaction"),
                    ~if_else(.=="",NA,.))) |> 
      rowwise() |> 
      mutate(list.pred=list(na.omit(c(rand,clim,quad_clim,predictor,quad_pred,interaction))),
             formula=case_when(random.eff=="none"~paste0(mod.nran,response,"~",
                                                         paste(list.pred,collapse = "+"),
                                                         mod_extension,
                                                         ",data=data_mod)"),
                               random.eff!="none"~paste0(mod.ran,response,"~",
                                                         paste(list.pred,collapse = "+"),
                                                         mod_extension,
                                                         ",data=data_mod)"))) |> 
      pull(formula)
    formulas<-c(formulas,formula_pred)
  }
  return(unique(formulas))
}


predict_traits<-setNames(data.frame(matrix(nrow = 0,ncol = 10)),
                         nm=c('data_name','response','pca_sc','trait','trait_value','predicted_mean','predicted_se','lwr','upr','pred'))
class(predict_traits$data_name)<-"character"
traits_effect<-setNames(data.frame(matrix(nrow = 0,ncol = 6)),
                        nm=c('data_name',"interaction",'trait','effect',"lwr","upr"))
}
### RUN

for (iter in 1:dim(run_mod)[1]){
  data_name=run_mod$data_name[iter]
  response_name=run_mod$response_name[iter]
  response_var=run_mod$response_var[iter]
  fixed_predictor=run_mod$fixed_predictor[iter]
  predictors_list=paste0(run_mod$predictors_list[iter],c("HM", "inv_sp","WD"))
  group_var=run_mod$group_var[iter] 
  mod_extension=run_mod$mod_extension[iter]
  mod.type=run_mod$mod.type[iter]
  
  eval(parse(text=paste0("data_mod<-",data_name)))
  data_mod[[response_var]]<-data_mod[[response_name]]
  
  formulas <- generate_formulas(response_var, fixed_predictor, predictors_list, group_var,mod_extension,mod.type)
  model_eval<-data.frame(formulas=unlist(formulas))|> 
    rowwise() |> 
    mutate(trait = purrr::map_chr(predictors_list, ~ ifelse(grepl(.x, formulas), .x, NA_character_)) %>% 
             purrr::discard(is.na) %>% 
             first()) %>%
    ungroup() |> 
    mutate(model=paste0("model_",row_number()),
           AIC=NA_real_,
           ncof=NA_real_)
  for(i in 1:dim(model_eval)[1]){
    print(paste0("model ",i,"/",dim(model_eval)[1]))
    eval(parse(text = paste0("model_i", "=",model_eval$formulas[i])))
    
    # Add AIC in the table
    model_eval$AIC[i]= AIC(model_i)
    if(mod.type=="lmer"){
      model_eval$ncof[i]=dim(summary(model_i)$coefficients)[1]
    }else{
      model_eval$ncof[i]=dim(summary(model_i)$coefficients$cond)[1]
      
    }
  }
  
  # best_mod
  min_aic=min(model_eval$AIC,na.rm=TRUE)
  max_aic=max(model_eval$AIC,na.rm=TRUE)
  best_model<-model_eval |> filter(AIC<min_aic+max(15,(max_aic-min_aic)/20)) |> 
    filter(ncof==min(ncof)) |> 
    filter(AIC==min(AIC)) |> 
    pull(formulas)
  
  
  best_model_trait<-model_eval |> 
    group_by(trait) |> 
    filter(AIC<min(AIC,na.rm=TRUE)+15) |> 
    filter(ncof==min(ncof)) |> 
    filter(AIC==min(AIC)) 
  
  for(i in 1:dim(best_model_trait)[1]){
    best_model=best_model_trait$formulas[i]
    
    mod_maint<-eval(parse(text = best_model))
    
    # Define new data for predictions
    if(grepl("pca_sc",best_model)){
      predictor=predictors_list[unlist(lapply(predictors_list,function(x)grepl(x,best_model)))]
      if(length(predictor)>0){
        new_data <- data.frame(
          response = 0,  # Dummy value, needed only for model.matrix
          pca_sc = seq(min(data_mod$pca_sc), max(data_mod$pca_sc), length.out = 15)  # Full range of pca_sc
        ) |> 
          crossing(pred = unname(quantile(data_mod[[predictor]],probs = c(0.05,0.5,0.95))))  # Use typical or specific values for pca1)
        colnames(new_data)[grepl("pred",colnames(new_data))]<-predictor
      }else{
        new_data <- data.frame(
          response = 0,  # Dummy value, needed only for model.matrix
          pca_sc = seq(min(data_mod$pca_sc), max(data_mod$pca_sc), length.out = 15)  # Full range of pca_sc
        ) 
      }
      
      
    }
    
    if(mod.type=="lmer"){
      fixed_effects <- fixef(mod_maint) # Extract fixed effects coefficients
      vcov_fixed <- vcov(mod_maint)
      X_fixed <- model.matrix(terms(mod_maint), new_data)
      new_data$predicted_mean <- X_fixed %*% fixed_effects
      fixed_var <- diag(X_fixed %*% vcov_fixed %*% t(X_fixed))
      new_data$predicted_se <- sqrt(fixed_var)
      new_data$lwr <- new_data$predicted_mean-1.96*new_data$predicted_se
      new_data$upr <- new_data$predicted_mean+1.96*new_data$predicted_se
    }else if(grepl("beta_family",mod_extension)){
      fixed_effects <- fixef(mod_maint)$cond  # Extract fixed effects coefficients
      vcov_fixed <- vcov(mod_maint)$cond  
      X_fixed <- model.matrix(terms(mod_maint), new_data)
      new_data$predicted_mean <- X_fixed %*% fixed_effects
      fixed_var <- diag(X_fixed %*% vcov_fixed %*% t(X_fixed))
      new_data$predicted_se <- sqrt(fixed_var)
      new_data$lwr <- plogis(new_data$predicted_mean-1.96*new_data$predicted_se)
      new_data$upr <- plogis(new_data$predicted_mean+1.96*new_data$predicted_se)
      new_data$predicted_mean <- plogis(X_fixed %*% fixed_effects)
    }else if(grepl("Gamma",mod_extension)){
      fixed_effects <- fixef(mod_maint)$cond  # Extract fixed effects coefficients
      vcov_fixed <- vcov(mod_maint)$cond  
      X_fixed <- model.matrix(terms(mod_maint), new_data)
      new_data$predicted_mean <- X_fixed %*% fixed_effects
      fixed_var <- diag(X_fixed %*% vcov_fixed %*% t(X_fixed))
      new_data$predicted_se <- sqrt(fixed_var)
      new_data$lwr <- exp(new_data$predicted_mean-1.96*new_data$predicted_se)
      new_data$upr <- exp(new_data$predicted_mean+1.96*new_data$predicted_se)
      new_data$predicted_mean <- exp(X_fixed %*% fixed_effects)
    }else{
      fixed_effects <- fixef(mod_maint)$cond  # Extract fixed effects coefficients  
      vcov_fixed <- vcov(mod_maint)$cond  
      X_fixed <- model.matrix(terms(mod_maint), new_data)
      new_data$predicted_mean <- X_fixed %*% fixed_effects
      fixed_var <- diag(X_fixed %*% vcov_fixed %*% t(X_fixed))
      new_data$predicted_se <- sqrt(fixed_var)
      new_data$lwr <- exp(new_data$predicted_mean-1.96*new_data$predicted_se)
      new_data$upr <- exp(new_data$predicted_mean+1.96*new_data$predicted_se)
      new_data$predicted_mean <- exp(X_fixed %*% fixed_effects)
    }
    if(length(predictor)>0){
      new_data$pred<-new_data[[predictor]]
      new_data$pred_name=predictor
      new_data<-new_data[,-match(predictor,colnames(new_data))]
    }else{
      new_data$pred<-NA
      new_data$pred_name<-"none"
    }
    
    predict_traits<-bind_rows(predict_traits,
                              new_data |> mutate(data_name=data_name) )
    if(length(predictor)>0){
      if(mod.type=="lmer"){
        effect_conf<-as.data.frame(confint(mod_maint)) |> 
          tibble::rownames_to_column(var="trait")
        mod_sum<-summary(mod_maint)$coefficients
        mean_effect<-mod_sum[as.logical(grepl(paste0("pca_sc:",predictor),rownames(mod_sum))+grepl(predictor,rownames(mod_sum))),
                             1]
        names_effect=rownames(mod_sum)[as.logical(grepl(paste0("pca_sc:",predictor),rownames(mod_sum))+grepl(predictor,rownames(mod_sum)))]
        vec_eff<-cbind(data_name,
                       "",
                       names_effect,
                       mean_effect,
                       unname(effect_conf[effect_conf$trait%in%names_effect,c(2,3)]))
        rownames(vec_eff)<-NULL
        colnames(vec_eff)<-colnames(traits_effect)[]
        vec_eff$interaction[grepl(":",vec_eff$trait)]<-TRUE
        vec_eff$interaction[!grepl(":",vec_eff$trait)]<-FALSE
        traits_effect<-rbind(traits_effect,
                             vec_eff)
        
      }else{
        effect_conf<-as.data.frame(confint(mod_maint)) |> 
          tibble::rownames_to_column(var="trait")
        vec_eff<-cbind(data_name,
                       "",
                       unname(effect_conf[as.logical(grepl(paste0("pca_sc:",predictor),effect_conf$trait)+grepl(predictor,effect_conf$trait)),
                                          c(1,4,2,3)]))
        rownames(vec_eff)<-NULL
        colnames(vec_eff)<-colnames(traits_effect)[]
        vec_eff$interaction[grepl(":",vec_eff$trait)]<-TRUE
        vec_eff$interaction[!grepl(":",vec_eff$trait)]<-FALSE
        vec_eff$trait=predictor
        traits_effect<-rbind(traits_effect,
                             vec_eff)
      }
    }
  }
}


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
