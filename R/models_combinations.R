## automatize test of different model forms and predictions of uncertainties


## maintenance ####
#%%%%%%%%%%%%%%%%%%

data_maint_sp<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  # rename(pca_sc=pca1) |> 
  filter(species==gsub("_"," ",species_combination)) |> 
  filter(metric=="ba_dif") |> 
  left_join(traits %>% mutate(species=gsub("_"," ",species)))  
data_name<-"data_maint_sp"
response_name <- "ba_target" 
response_var <- "response"
eval(parse(text=paste0("data_mod<-",data_name)))
data_mod[[response_var]]<-data_mod[[response_name]]

fixed_predictor <- "pca_sc"
predictors_list <- c("HM", "inv_sp","WD")
group_var <- "species"
mod_extension=""
mod.type="lmer"

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

data_name<-"data_maint"
response_name <- "metric_val" 
response_var <- "response"
eval(parse(text=paste0("data_mod<-",data_name)))
data_mod[[response_var]]<-data_mod[[response_name]]
fixed_predictor <- "pca_sc"
predictors_list <- c("nih_HM", "nih_inv_sp","nih_WD")
group_var <- "species"
mod_extension=",family = beta_family(link = \"logit\")"
mod.type="glmmTMB"

## resilience ####
#%%%%%%%%%%%%%%%%%
data_resilience_sp<- performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(!is.na(metric_val)) |> 
  filter(metric=="resilience") %>%
  filter(vr=="mean") %>%
  filter(species==gsub("_"," ",species_combination)) |> 
  mutate(metric_val=ba_target*metric_val,
         metric_val=log(metric_val)) 
data_name<-"data_resilience_sp"
response_name <- "metric_val" 
response_var <- "response"
eval(parse(text=paste0("data_mod<-",data_name)))
data_mod[[response_var]]<-data_mod[[response_name]]
fixed_predictor <- "pca_sc"
predictors_list <-  c("HM", "inv_sp","WD")
group_var <- "species"
mod_extension=""
mod.type="glmmTMB"


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
data_name<-"data_resilience"
response_name <- "res_dif" 
response_var <- "response"
eval(parse(text=paste0("data_mod<-",data_name)))
data_mod[[response_var]]<-data_mod[[response_name]]
fixed_predictor <- "pca_sc"
predictors_list <- c("nih_HM", "nih_inv_sp","nih_WD")
group_var <- "species"
mod_extension=""
mod.type="glmmTMB"

## invasion ####
#%%%%%%%%%%%%%%%
data_inv_sp<-performance %>% 
  left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
  left_join(traits %>% mutate(species=gsub("_"," ",species))) %>% 
  filter(metric=="inv_50") %>% 
  filter(species==gsub("_"," ",species_combination))
data_name<-"data_inv_sp"
response_name <- "metric_val" 
response_var <- "response"
eval(parse(text=paste0("data_mod<-",data_name)))
data_mod[[response_var]]<-data_mod[[response_name]]
fixed_predictor <- "pca_sc"
predictors_list <-  c("HM", "inv_sp","WD")
group_var <- "species"
mod_extension=""
mod.type="glmmTMB"

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
data_name<-"data_inv"
response_name <- "inv_dif" 
response_var <- "response"
eval(parse(text=paste0("data_mod<-",data_name)))
data_mod[[response_var]]<-data_mod[[response_name]]
fixed_predictor <- "pca_sc"
predictors_list <- c("nih_HM", "nih_inv_sp","nih_WD")
group_var <- "species"
mod_extension=",family = Gamma(link = \"log\")" #",family = beta_family(link = \"logit\")"
mod.type="glmmTMB"


## Functions ####
#%%%%%%%%%%%%%%%%
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

## RUN ####
#%%%%%%%%%#

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

predict_traits<-setNames(data.frame(matrix(nrow = 0,ncol = 9)),
                         nm=c('response','pca_sc','trait','trait_value','predicted_mean','predicted_se','lwr','upr','pred'))
traits_effect<-setNames(data.frame(matrix(nrow = 0,ncol = 4)),
                        nm=c('trait','effect',"lwr","upr"))
for(i in 1:dim(best_model_trait)[1]){
  best_model=best_model_trait$formulas[i]
  
  mod_maint<-eval(parse(text = best_model))
  
  # sim_res <- simulateResiduals(mod_maint)
  # plot(sim_res)
  
  
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
                            new_data )
  if(length(predictor)>0){
    effect_conf<-as.data.frame(confint(mod_maint)) |> 
      tibble::rownames_to_column(var="trait")
    vec_eff<-unname(effect_conf[effect_conf$trait==predictor,c(1,4,2,3)])
    rownames(vec_eff)<-NULL
    colnames(vec_eff)<-colnames(traits_effect)
    traits_effect<-rbind(traits_effect,
                         vec_eff)
    
  }
  # new_data |> 
  #   ggplot()+
  #   geom_ribbon(aes(x=pca_sc,ymax=upr,ymin=lwr,fill=as.factor(pred)),alpha=0.2)+
  #   geom_line(aes(pca_sc,predicted_mean,color=as.factor(pred)))+
  #   labs(color=paste0("Categories of ",predictor),
  #        fill=paste0("Categories of ",predictor),
  #        y=response_name,
  #        title=paste0(best_model))->plot
  # print(plot)
  
}

predict_traits |> 
  group_by(pred_name) |> 
  mutate(pred_cat=case_when(pred==min(pred)~"High compet",
                            pred==max(pred)~"Low compet",
                            TRUE~"Medium compet"),
         pred_cat=factor(pred_cat,levels=c("Low compet","Medium compet","High compet")),
         pred_name=factor(pred_name,levels=c("none",predictors_list))) |> 
  filter(pred_cat%in%c("Low compet","High compet")) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymax=upr,ymin=lwr,
                  fill=pred_cat),
              alpha=0.2)+
  geom_line(aes(pca_sc,predicted_mean,color=pred_cat),
            size=1)+
  labs(color=paste0("Categories of ",predictor),
       fill=paste0("Categories of ",predictor),
       y=response_name,
       title=paste0(best_model))+
  scale_color_manual(values=c("darkseagreen","darksalmon","darkred"))+
  scale_fill_manual(values=c("darkseagreen","darksalmon","darkred"))+
  facet_wrap(~pred_name)+
  theme_classic()

traits_effect |> 
  filter(grepl("sp",data_name)) |> 
  ggplot()+
  geom_segment(aes(x=lwr,xend=upr,y=trait,yend=trait))+
  geom_point(aes(x=effect,y=trait))+
  theme_classic()+
  labs(y="",x="Standardized trait effect")+
  facet_wrap()
 