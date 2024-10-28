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


response_var <- "ba_target"
fixed_predictor <- "pca_sc"
predictors_list <- c("shade",  "inv_sp","pca1")
group_var <- "species"
data_name<-"data_maint_sp"
generate_formulas <- function(response_var, fixed_predictor, predictors_list, group_var,data_name) {
  formulas<-list()
  formulas[[1]] <- paste0("lmer(",response_var,"~ (1|", group_var, ")",",data=",data_name,")")
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
             clim=case_when(random.eff=="slope"&random.var=="clim"~paste0("(",clim,formula.random),
                            TRUE~clim),
             quad_clim=case_when(random.eff=="slope"&random.var=="quad_clim"~paste0("(",quad_clim,formula.random),
                            TRUE~quad_clim),
             predictor=case_when(random.eff=="slope"&random.var=="predictor"~paste0("(",predictor,formula.random),
                            TRUE~predictor),
             quad_pred=case_when(random.eff=="slope"&random.var=="quad_pred"~paste0("(",quad_pred,formula.random),
                            TRUE~quad_pred),
             rand=case_when(random.eff=="intercept"~formula.random,
                            TRUE~NA),
             across(c("clim","quad_clim","predictor","quad_pred","interaction"),
                    ~if_else(.=="",NA,.))) |> 
      rowwise() |> 
      mutate(list.pred=list(na.omit(c(rand,clim,quad_clim,predictor,quad_pred,interaction))),
             formula=case_when(random.eff=="none"~paste0("lm(",response,"~",
                                                         paste(list.pred,collapse = "+"),
                                                         ",data=",data_name,")"),
                               random.eff!="none"~paste0("lmer(",response,"~",
                                                         paste(list.pred,collapse = "+"),
                                                         ",data=",data_name,")"))) |> 
      pull(formula)
    formulas<-c(formulas,formula_pred)
    }
  return(unique(formulas))
}

# Example usage

formulas <- generate_formulas(response_var, fixed_predictor, predictors_list, group_var)
model_eval<-data.frame(formulas=unlist(formulas)) |> 
  mutate(model=paste0("model_",row_number()),
         AIC=NA_real_,
         ncof=NA_real_)
for(i in 1:dim(model_eval)[1]){
  eval(parse(text = paste0("model_i", "=",model_eval$formulas[i])))
  
  # Add AIC in the table
  model_eval$AIC[i]= AIC(model_i)
  model_eval$ncof[i]=dim(summary(model_i)$coefficients)[1]
}
