## automatize test of different model forms and predictions of uncertainties

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

data_maint_sp |> 
  ggplot(aes(pca_sc,ba_target))+
  geom_point()+
  geom_smooth()


response_var <- "metric_val"
fixed_predictor <- "pca_sc"
predictors_list <- c("nih_shade", "nih_inv_sp","nih_pca1")
# predictors_list <- c("shade", "inv_sp","pca1")
group_var <- "species"
data_name<-"data_maint"
mod_extension=",family = beta_family(link = \"logit\")"
mod.type="glmmTMB"
generate_formulas <- function(response_var, fixed_predictor, predictors_list,
                              group_var,data_name,mod_extension,mod.type) {
  formulas<-list()
  if(mod.type=="lmer"){
    mod.ran="lmer("
    mod.nran="lm("
  }else{
    mod.ran="glmmTMB("
    mod.nran="glmmTMB("
  }
  formulas[[1]] <- paste0(mod.ran,response_var,"~ (1|", group_var, ")",mod_extension,",data=",data_name,")")
  for (pred in predictors_list) {
    formula_pred<-data.frame(model="model",
                         response=response_var,
                         clim=fixed_predictor) |> 
      crossing(quad_clim=c("",paste0("I(",fixed_predictor,"^2)"))) |> 
      crossing(predictor=c("",pred)) |>
      crossing(quad_pred=c("",paste0("I(",pred,"^2)") )) |> 
      crossing(data.frame(random.eff=c("none","intercept","slope"),
                          formula.random=c("",paste0("(1|", group_var, ")"),paste0( "|", group_var, ")")))) |> 
      crossing(random.var=c("clim","quad_clim","predictor","quad_pred")) |> 
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
                                   interaction=="clim"~paste0(clim,":",predictor),
                                   interaction=="quad_clim"~paste0(quad_clim,":",predictor),
                                   interaction=="both"~paste0(clim,":",predictor,"+",quad_clim,":",predictor)),
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
                                                         ",data=",data_name,")"),
                               random.eff!="none"~paste0(mod.ran,response,"~",
                                                         paste(list.pred,collapse = "+"),
                                                         mod_extension,
                                                         ",data=",data_name,")"))) |> 
      pull(formula)
    formulas<-c(formulas,formula_pred)
    }
  return(unique(formulas))
}

# Example usage

formulas <- generate_formulas(response_var, fixed_predictor, predictors_list, group_var,data_name,mod_extension,mod.type)
model_eval<-data.frame(formulas=unlist(formulas)) |> 
  mutate(model=paste0("model_",row_number()),
         AIC=NA_real_,
         ncof=NA_real_)
for(i in 1:dim(model_eval)[1]){
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
min_aic=min(model_eval$AIC)

best_model<-model_eval |> filter(AIC<min_aic+10) |> 
  filter(ncof==min(ncof)) |> 
  filter(AIC==min(AIC)) |> 
  pull(formulas)


mod_maint<-eval(parse(text = best_model))

sim_res <- simulateResiduals(mod_maint)
plot(sim_res)

# Define new data for predictions
if(grepl("pca_sc",best_model)){
  predictor=predictors_list[unlist(lapply(predictors_list,function(x)grepl(x,best_model)))]
  if(length(predictor)>0){
    new_data <- data.frame(
      metric_val = 0,  # Dummy value, needed only for model.matrix
      pca_sc = seq(min(data_maint$pca_sc), max(data_maint$pca_sc), length.out = 15)  # Full range of pca_sc
    ) |> 
      crossing(pred = unname(quantile(data_maint[[predictor]],probs = c(0.1,0.5,0.9))))  # Use typical or specific values for pca1)
    colnames(new_data)[grepl("pred",colnames(new_data))]<-predictor
  }else{
    new_data <- data.frame(
      metric_val = 0,  # Dummy value, needed only for model.matrix
      pca_sc = seq(min(data_maint$pca_sc), max(data_maint$pca_sc), length.out = 15)  # Full range of pca_sc
    ) 
  }
  
  
}

# Obtain fixed effect predictions
fixed_effects <- as.matrix(fixef(mod_maint))  # Extract fixed effects coefficients

# Design matrix for the new data based on fixed effects
X_fixed <- model.matrix(terms(mod_maint), new_data)

# Predict mean values using fixed effects only
new_data$predicted_mean <- X_fixed %*% fixed_effects

# Compute standard errors for the fixed effect predictions
# Get the variance-covariance matrix for the fixed effects
vcov_fixed <- vcov(mod_maint)
fixed_var <- diag(X_fixed %*% vcov_fixed %*% t(X_fixed))
new_data$predicted_se <- sqrt(fixed_var)

# Display the results
# (aes(x=pca_sc,ymin=lower,ymax=upper,fill=as.factor(nih_pca1)),alpha=0.2)+
new_data |> 
  mutate(lwr=predicted_mean-1.96*predicted_se,
         upr=predicted_mean+1.96*predicted_se) |> 
  ggplot()+
  geom_ribbon(aes(x=pca_sc,ymax=upr,ymin=lwr,fill=as.factor(pca1)),alpha=0.2)+
  geom_line(aes(pca_sc,predicted_mean,color=as.factor(pca1)))
