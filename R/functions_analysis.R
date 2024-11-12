#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analysis.R  
#' @description R script containing all functions relative to formatting
#               and analysis of models outputs
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%
#### GET TRAITS ####
#%%%%%%%%%%%%%%%%%%%


#' Function to load Shade tolerance from niinemets&Valladares
#' @description only outs evergreen and shadetol
#' @param file.shade file path of database

get_shade<-function(file.shade="data/data_Niinemets&Valladares_2006.csv"){
  shade=read.csv2(file.shade) |> 
    rename(species=Species,
           evergreen=Evergreen) |> 
    mutate(shade_tolerance.mean=as.numeric(shade_tolerance.mean)) |> 
    dplyr::select(species,evergreen,shade_tolerance.mean)
  return(shade)
}

#' Function to load leaf traits and wood density from a TRY query
#' @description Needs to perform a TRY query on TRY website to get the .csv file
#' @param file.try file of the try query csv

get_try<-function(file.try){
  
  ## load try query
  try_data<-read.csv2(file.try)
  
  ## 
  try_format<-try_data %>% 
    filter(!is.na(TraitID)) %>% 
    filter(TraitID!=18) %>% 
    dplyr::select(AccSpeciesName,ObservationID,ObsDataID,TraitID,TraitName,StdValue) %>% 
    unique() %>% 
    mutate(StdValue=as.numeric(StdValue))
  # table(try_format[try_format$TraitID==3116,]$UnitName) # checked if units were consistent across traits
  
  ## reformat
  leaf_trait<-try_format %>%
    mutate(species=case_when(AccSpeciesName=="Betula pendula"~"Betula",
                             AccSpeciesName=="Betula pubescens"~"Betula",
                             TRUE~AccSpeciesName)) %>% 
    group_by(species,TraitName) %>% 
    filter(!(TraitName=="Seed dry mass"&StdValue>quantile(StdValue,probs=0.95,na.rm=TRUE))) %>% 
    filter(!(TraitName=="Seed dry mass"&StdValue<quantile(StdValue,probs=0.05,na.rm=TRUE))) %>% 
    summarise(trait_mean=median(StdValue,na.rm=TRUE),
              trait_sd=sd(StdValue,na.rm=TRUE)) %>% 
    mutate(trait_mean=case_when(is.nan(trait_mean)~NA,
                                TRUE~trait_mean),
           trait_sd=case_when(is.nan(trait_sd)~NA,
                              TRUE~trait_sd),
           cv=trait_sd/trait_mean) %>% 
    ungroup() %>% 
    rename(trait=TraitName) |> 
    mutate(trait=case_when(grepl("SLA",trait)~"SLA",
                           grepl("Leaf nitrogen",trait)~"LN",
                           grepl("Leaf thickness",trait)~"LT",
                           grepl("Seed dry mass",trait)~"SDM",
                           grepl("SSD",trait)~"WD"))
  return(leaf_trait)
}


#' Function to load wood density from try and global wd database
#' @description if wd is missing in GWDD then takes the one of TRY
#' @param file.wd file path to GWDD
#' @param trait_try traits from try query

get_WD<-function(file.wd,
                 trait_try){
  wd_global<-read.csv2(file.wd) %>% 
    # filter(Binomial%in%sp_list) %>% 
    mutate(species=case_when(Binomial=="Betula pendula"~"Betula",
                             Binomial=="Betula pubescens"~"Betula",
                             TRUE~Binomial)) %>% 
    group_by(species) %>% 
    summarise(WD=mean(Wood.density))
  
  wd_try<-trait_try |> 
    filter(trait=="WD") |> 
    pivot_wider(names_from = trait,values_from = trait_mean) |> 
    dplyr::select(species,WD)
  
  sp_list<-sort(unique(c(wd_try$species,wd_global$species)))
  
  wd_tot<-data.frame(species=sp_list) %>% 
    left_join(wd_global,by='species') %>%
    rename(WD_glob=WD) |> 
    left_join(wd_try,by="species") %>% 
    rename(WD_try=WD) |> 
    mutate(WD=case_when(is.na(WD_glob)~WD_try,
                        TRUE~WD_glob)) |> 
    dplyr::select(species,WD)
  
  return(wd_tot)
}

#' Function to get maximum height from NFI
#' @param FUNDIV_data FUNDIV NFIs
get_maxH<-function(FUNDIV_data){
  max_height<- FUNDIV_data %>% unique() %>% 
    # filter(species%in%c(sp_list,"Betula")) %>% 
    filter(height1>0) %>% 
    group_by(species) %>% 
    filter(!is.na(height1)) %>% 
    filter(!is.na(weight1)) %>% 
    summarise(max_height=modi::weighted.quantile(height1,w=weight1,prob=0.99)) 
  return(max_height)
}

#' Function to get recruitment traits
#' @description applies functions of recruitment at mean ba 
#' @author Julien Barrere
#' @param FUNDIV_data FUNDIV NFIs
get_recruitment_traits = function(fit.list.allspecies, FUNDIV_data, comp.ref){
  
  # Remove in-growth
  FUNDIV_data = subset(FUNDIV_data, treestatus_th != 1)
  
  # Calculate the mean competition in the entire dataset
  BASP.mean = mean(FUNDIV_data$BAtotSP, na.rm = TRUE)
  BANONSP.mean = mean(FUNDIV_data$BAtotNONSP, na.rm = TRUE)
  
  # If we choose to attribute the same competition for all species: 
  if(comp.ref == "same"){
    # Provide mean value of competition across all species
    data_species = data.frame(species = unique(climate_species$sp), 
                              BATOTSP = BASP.mean, 
                              BATOTNonSP = BANONSP.mean, 
                              logBATOTSP = log(BASP.mean), 
                              intercept = 1) %>%
      # Add climatic data
      left_join((climate_species %>%
                   filter(N == 2) %>%
                   dplyr::select(species = sp, wai, wai2, 
                                 waib, sgdd, sgdd2, sgddb)), 
                by = "species")
  }
  
  # If we choose to attribute different competition per species
  if(comp.ref == "specific"){
    # Calculate mean competition per species
    data_species = FUNDIV_data %>%
      mutate(species = gsub("\\ ", "\\_", species)) %>%
      filter(species %in% unique(climate_species$sp)) %>%
      group_by(species) %>%
      summarize(BATOTSP = mean(BAtotSP, na.rm = TRUE), 
                BATOTNonSP = mean(BAtotNONSP, na.rm = TRUE)) %>%
      mutate(logBATOTSP = log(BATOTSP), 
             intercept = 1) %>%
      # Add climatic data
      left_join((climate_species %>%
                   filter(N == 2) %>%
                   dplyr::select(species = sp, wai, wai2, 
                                 waib, sgdd, sgdd2, sgddb)), 
                by = "species")
  }
  
  # Initialize the output dataset
  traits_rec = data.frame(species = names(fit.list.allspecies), 
                          recruitment = NA_real_, 
                          delay = NA_real_)
  
  # Loop on all species to gather traits
  for(i in 1:dim(traits_rec)[1]){
    
    # Species i
    sp.i = traits_rec$species[i]
    
    # Give the delay from fit list
    traits_rec$delay[i] = as.numeric(fit.list.allspecies[[i]]$info["delay"])
    
    # Recruitment parameters for species i
    vec.rec.i = fit.list.allspecies[[i]]$rec$params_m
    
    # Associated vector of variables
    vec.var.i = as.vector(subset(data_species, species == sp.i)[, names(vec.rec.i)])
    
    # Get the recruitment
    traits_rec$recruitment[i] = exp(sum(vec.rec.i*unlist(vec.var.i)))
    
  }
  
  # Return the trait dataset generated
  return(traits_rec)
}


#' Function to gather all traits
#' @description 
get_traits<-function(trait_shade,
                     trait_jul,
                     trait_try,
                     trait_wd,
                     trait_maxH,
                     trait_recruitment,
                     species.list.ipm){
  sp_list=gsub("_"," ",species.list.ipm)
  trait_jul$species=gsub("_"," ",trait_jul$species)
  trait_recruitment$species=gsub("_"," ",trait_recruitment$species)
  
  trait_try<-trait_try |> 
    dplyr::select(species,trait,trait_mean) |> 
    pivot_wider(names_from = trait,
                values_from = trait_mean)
  
  trait_tot<-data.frame(species=sp_list) |> 
    left_join(trait_jul[,c("species","height.dbh.ratio","max.growth")]) |> 
    left_join(trait_maxH) |>
    left_join(trait_recruitment) |> 
    left_join(trait_shade) |>  
    left_join(trait_try) |> 
    left_join(trait_wd) |> 
    rename(HM=max_height)
}


complete_trait<-function(trait_raw){
  imputed_data <- mice::mice(trait_raw, method = 'pmm', m = 5, maxit = 50, seed = 500)
  trait_complete <- mice::complete(imputed_data, 1)  # '1' refers to the first imputed dataset
  return(trait_complete)
}


#%%%%%%%%%%%%%%%%%%%%%%%%
#### GET PERFORMANCE ####
#%%%%%%%%%%%%%%%%%%%%%%%%

get_performance<-function(species.list.ipm,
                          disturbance_ba,
                          invasion_ba_2,
                          ba_dif,
                          sim_forest_excl,
                          climate.cat,
                          trait_complete,
                          trait_selection=c("HM","recruitment","WD")
                          ){
  
  traits<-trait_complete |> 
    dplyr::select(species,matches(trait_selection)) |> 
    mutate(species=gsub(" ","_",species))
  
  mean_pca<-climate.cat$FUNDIV_plotcat %>% 
    group_by(species,clim_id) %>% 
    summarise(pca1=mean(pca1),pca2=mean(pca2)) %>% 
    ungroup() %>% 
    group_by(species) %>% 
    mutate(pca_sc=scale(pca1))
  
  species.combination.excl<-sim_forest_excl$list.forests %>% 
    filter(forest.real) %>% 
    dplyr::select(species,clim_id,species_combination,excluded,competexcluded) %>% unique()
  
  species.combination.excl$smallcombi<-NA
  for(i in 1:dim(species.combination.excl)[1]){
    if(species.combination.excl$competexcluded[i]!=""){
      s_p=gsub(" ","_",species.combination.excl$species[i])
      sp_combi=unlist(strsplit(species.combination.excl$species_combination[i], split = "\\."))
      # sp_combi=sp_combi[sp_combi!=s_p]
      sp_ex=unlist(strsplit(species.combination.excl$competexcluded[i],"\\."))[-1]
      
      new_combi=paste(sort(sp_combi[!sp_combi%in%sp_ex]),collapse=".")
      ex<-species.combination.excl %>%
        filter(species==species.combination.excl$species[i]&
                 clim_id==species.combination.excl$clim_id[i]&
                 species_combination==new_combi)
      if(dim(ex)[1]==1){
        species.combination.excl$smallcombi[i]="present"
      }else{
        species.combination.excl$smallcombi[i]="absent"
      }
    }
    
  }
  
  ## ba_dif
  ba_dif<-ba_dif %>%
    dplyr::select(!matches(c('ba_target',"ba_partner","n_species"))) %>%
    mutate(elast="Abies_alba-amean",vr="mean")
  
  
  
  ## disturbance
  disturbance_metric<-disturbance_ba %>%
    mutate(elast="Abies_alba-amean",vr="mean") %>%
    dplyr::select(species,elast,species_combination,clim_id,vr,resilience,resistance,recovery,matches(species.list.ipm))
  
  disturbance <- disturbance_metric %>%
    arrange(species,clim_id,species_combination,elast) %>%
    pivot_longer(cols=matches(species.list.ipm)) %>%
    left_join(traits, by = c("name"="species")) %>%
    group_by(species,elast,species_combination,clim_id,vr) %>%
    mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                                 TRUE~0)) %>%
    pivot_longer(cols=colnames(traits)[-1],
                 values_to = "trait_val",
                 names_to = "trait_name") %>%
    group_by(species,elast,species_combination,clim_id,vr,trait_name) %>%
    mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
           ba_target=sum(value*ba_multiple,na.rm = TRUE),
           rel_ba_partner=ba_partner/(ba_target+ba_partner),
           nih=sum((sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner,
           nid=sum(abs(sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner) %>%
    ungroup() %>%
    dplyr::select(-c("ba_multiple","trait_val")) %>%
    pivot_wider(values_from = c("nih","nid"),
                names_from = "trait_name") %>%
    pivot_wider(names_from = name,
                values_from = value) %>%
    dplyr::select(!matches(species.list.ipm))
  
  ## invasion
  invasion_metric<-invasion_ba_2 %>%
    mutate(elast="Abies_alba-amean",vr="mean") %>% ## attention erruer
    dplyr::select(species,elast,species_combination,clim_id,vr,inv_mean,inv_max,inv_50,BA_100,BA_500,BA_1000,matches(species.list.ipm))
  
  invasion <- invasion_metric %>% ungroup() %>%
    arrange(species,species_combination,clim_id) %>%
    pivot_longer(cols=matches(species.list.ipm)) %>%
    left_join(traits, by = c("name" = "species")) %>%
    group_by(species,elast,species_combination,clim_id,vr) %>%
    mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                                 TRUE~0),
           value=case_when(!is.na(value)~1,
                           gsub("_"," ",name)==species~1,
                           TRUE~value)) %>%  # *1 because invasion are made at constant basal area of compet
    pivot_longer(cols=colnames(traits)[-1],
                 values_to = "trait_val",
                 names_to = "trait_name") %>%
    group_by(species,elast,species_combination,clim_id,vr,trait_name) %>%
    mutate(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
           ba_target=sum(value*ba_multiple,na.rm = TRUE),
           rel_ba_partner=ba_partner/(ba_target+ba_partner),
           nih=sum((sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner,
           nid=sum(abs(sum(trait_val*ba_multiple)-trait_val)*value,na.rm=TRUE)/ba_partner) %>%
    ungroup() %>%
    dplyr::select(-c("ba_multiple","trait_val")) %>%
    pivot_wider(values_from = c("nih","nid"),
                names_from = "trait_name") %>%
    pivot_wider(names_from = name,
                values_from = value) %>%
    dplyr::select(!matches(species.list.ipm))%>%
    left_join(species.combination.excl,by=c("species","clim_id","species_combination"))
  
  
  ## gather all
  performance <- invasion %>%
    left_join(disturbance,by=c("species","elast","species_combination","clim_id","vr")) %>%
    left_join(ba_dif,by=c("species","elast","species_combination","clim_id","vr"))%>%
    arrange(species,clim_id,species_combination,elast,vr) %>%
    pivot_longer(cols=c("resilience","recovery","resistance","inv_mean","inv_max",
                        "inv_50","BA_100","BA_500","BA_1000","ba_dif"),
                 names_to="metric",
                 values_to="metric_val") %>%
    arrange(species,clim_id,elast,species_combination) |> 
    # matches the correct values with metrics
    pivot_longer(cols=matches("\\.x|\\.y")) %>%
    mutate(name_upd=gsub("\\.x$|\\.y$", "", name),
           name=str_sub(name,-1)) %>% 
    pivot_wider(names_from = name,
                values_from = value) %>%
    mutate(compet_val=case_when(metric%in%c("inv_mean","inv_max","inv_50","BA_100","BA_500","BA_1000")~x,
                                TRUE~y)) %>%
    dplyr::select(-c("x","y")) %>%
    pivot_wider(names_from = name_upd,
                values_from = compet_val)
  return(performance)
}


#%%%%%%%%%%%%%%%%%%%%%%%%
#### GET PERFORMANCE ####
#%%%%%%%%%%%%%%%%%%%%%%%%
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
      filter(kp) |> dplyr::select(-kp) |> 
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

#' Function that fits model over all demo perf
#' @description description
run_model<-function(performance,
                    climate.cat,
                    trait_complete){
  mean_pca<-climate.cat$FUNDIV_plotcat %>% 
    group_by(species,clim_id) %>% 
    summarise(pca1=mean(pca1),pca2=mean(pca2)) %>% 
    ungroup() %>% 
    group_by(species) %>% 
    mutate(pca_sc=scale(pca1))
  
  data_maint_sp<-performance %>% 
    left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
    # rename(pca_sc=pca1) |> 
    filter(species==gsub("_"," ",species_combination)) |> 
    filter(metric=="ba_dif") |> 
    left_join(trait_complete %>% mutate(species=gsub("_"," ",species)))  
  
  data_maint<-performance %>% 
    left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>%
    filter(metric=="ba_dif") %>% 
    filter(!is.nan(nih_HM)) %>% 
    left_join(trait_complete %>% mutate(species=gsub("_"," ",species))) %>% 
    mutate(species=forcats::fct_reorder(species, HM),
           simul_state=case_when(!is.na(smallcombi)~"CompetitorExclusion",
                                 excluded=="excluded"~"SpeciesExclusion",
                                 TRUE~"SpeciesCoex"),
           metric_val=case_when(metric_val>1~1,
                                TRUE~metric_val),
           metric_val=(metric_val * (dim(.)[1] - 1) + 0.5) / dim(.)[1]) 
  
  data_resilience_sp<- performance %>% 
    left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
    left_join(trait_complete %>% mutate(species=gsub("_"," ",species))) %>% 
    filter(!is.na(metric_val)) |> 
    filter(metric=="resilience") %>%
    filter(vr=="mean") %>%
    filter(species==gsub("_"," ",species_combination)) |> 
    mutate(metric_val=ba_target*metric_val,
           metric_val=log(metric_val)) 
  
  data_resilience<- performance %>% 
    left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
    left_join(trait_complete %>% mutate(species=gsub("_"," ",species))) %>% 
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
    left_join(trait_complete %>% mutate(species=gsub("_"," ",species))) %>% 
    filter(metric=="inv_50") %>% 
    filter(species==gsub("_"," ",species_combination))
  
  
  data_inv<- performance %>% 
    left_join(mean_pca[,c("species","clim_id","pca_sc")]) %>% 
    filter(metric=="inv_50") %>% 
    left_join(trait_complete %>% mutate(species=gsub("_"," ",species))) %>% 
    rowwise() %>% 
    mutate(n_species=length(strsplit(species_combination,"\\.")[[1]])) %>% 
    ungroup() %>% 
    group_by(species,clim_id) %>% 
    arrange(species,clim_id,n_species) %>% 
    mutate(inv_dif=metric_val/metric_val[1],
           linv_dif=log(inv_dif)) %>% 
    ungroup() %>% 
    # mutate(nih=case_when(is.nan(nih_inv_sp)~0,
    #                      TRUE~nih_inv_sp)#,
    #        # inv_dif=case_when(inv_dif<0~0,
    #        #                   TRUE~inv_dif),
    #        # inv_dif=(inv_dif * (dim(.)[1] - 1) + 0.5) / dim(.)[1]
    # ) %>% 
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
  
  ### run models
  
  for (iter in 1:dim(run_mod)[1]){
    data_name=run_mod$data_name[iter]
    response_name=run_mod$response_name[iter]
    response_var=run_mod$response_var[iter]
    fixed_predictor=run_mod$fixed_predictor[iter]
    predictors_list=paste0(run_mod$predictors_list[iter],c("HM", "recruitment","WD"))
    group_var=run_mod$group_var[iter] 
    mod_extension=run_mod$mod_extension[iter]
    mod.type=run_mod$mod.type[iter]
    
    eval(parse(text=paste0("data_mod<-",data_name)))
    data_mod[[response_var]]<-data_mod[[response_name]]
    
    formulas <- generate_formulas(response_var, fixed_predictor, predictors_list, group_var,mod_extension,mod.type)
    model_eval <- data.frame(formulas = unlist(formulas)) |>
      rowwise() |>
      mutate(
        trait = {
          matched_traits <- purrr::map_chr(predictors_list, ~ ifelse(grepl(.x, formulas), .x, NA_character_)) %>% 
            purrr::discard(is.na)
          if (length(matched_traits) == 0) NA_character_ else data.table::first(matched_traits)
        }
      ) |>
      ungroup() |>
      mutate(
        model = paste0("model_", row_number()),
        AIC = NA_real_,
        ncof = NA_real_
      )
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
    best_models<-rbind(best_models,
                       best_model_trait |> mutate(data_name=data_name))
    
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
          mean_effect<-mod_sum[as.logical(grepl("2",rownames(mod_sum))+grepl(predictor,rownames(mod_sum))),
                               1]
          names_effect=rownames(mod_sum)[as.logical(grepl("2",rownames(mod_sum))+grepl(predictor,rownames(mod_sum)))]
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
                         unname(effect_conf[as.logical(!grepl("2",effect_conf$trait)&grepl(predictor,effect_conf$trait)),
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
  
  return(list(best_models=best_models,
              traits_effect=traits_effect,
              predict_traits=predict_traits))
}
