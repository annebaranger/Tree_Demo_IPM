intensity_file= "data/data/Disturbance/jags_dominance.Rdata"


#### Traits ####
# -- Traits data file
traitsNFI_file="data/data/Traits/traitsNFI.csv"
wood.density_file="data/Traits/GlobalWoodDensityDatabase.xls", format = "file"),
tar_target(shade.tolerance_file, "data/Traits/shade_tolerance_FrenchNFI.csv", format = "file"),
tar_target(TRY_file, "data/Traits/TRY_data_request_21092.txt", format = "file"),
# -- Compile traits from different databases
tar_target(traits_compiled, compile_traits(
  wood.density_file, traitsNFI_file, shade.tolerance_file, TRY_file, sp.in.sim)),
# -- Get all species included in the simulations for filtering
tar_target(sp.in.sim, gsub("\\ ", "\\_", unique(species_list$species))),
rdata.file=intensity_file  
# Load rdata
  load(rdata.file)
  
  # Loop on all disturbances
  for(i in 1:length(names(jags.list))){
    
    # Format the table containing parameter value per species, disturbance and iteration
    param.table.i <- ggs(as.mcmc(jags.list[[i]])) %>%
      # Add disturbance
      mutate(disturbance = names(jags.list)[i]) %>%
      # Extract information on parameter, species and country
      mutate(Param = gsub("\\[.+", "", Parameter), 
             sp = as.integer(ifelse(Param == "a0", gsub(".+\\[", "", gsub("\\,.+", "", Parameter)), 
                                    gsub(".+\\[", "", gsub("\\]", "", Parameter)))), 
             co = as.integer(ifelse(Param == "a0", gsub(".+\\,", "", gsub("\\]", "", Parameter)), NA_integer_))) %>%
      # Remove the estimation of intensity and deviance
      filter(Param != "I") %>%
      filter(Param != "deviance") %>%
      # Add name of the country and species, and weight of each species per country
      left_join(corresp.tables[[i]]$country.table, by = "co") %>%
      left_join(corresp.tables[[i]]$species.table, by = "sp") %>%
      left_join(weight.tables[[i]], by = c("species", "country")) %>%
      # No weight for the parameters that do not rely on the country
      mutate(weight = ifelse(Param == "a0", weight, 1)) %>%
      mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
      # Summarize Parameter value per country (only apply to a0)
      group_by(disturbance, Iteration, Chain, Param, species) %>%
      summarize(val = sum(value*weight, na.rm = TRUE)/sum(weight, na.rm = TRUE)) %>%
      # Format to get one column per parameter
      spread(key = Param, value = val) %>%
      # Set a1 to 0 (dominance effect) if disturbance is not storm or snow
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, 0)) %>%
      # Add parameters to scale dbh and logratio
      mutate(dbh.intercept = scale.tables[[i]]$dbh.intercept, 
             dbh.slope = scale.tables[[i]]$dbh.slope, 
             logratio.intercept = scale.tables[[i]]$logratio.intercept, 
             logratio.slope = scale.tables[[i]]$logratio.slope)
    
    # Store table in the final table
    if(i == 1) param.table <- param.table.i
    else param.table <- rbind(param.table, param.table.i)
  }
  
  # return the list
  return(param.table)
  
param=param.table

extend_disturb_coef = function(intensity_file, traits_compiled){
  
  # Get the parameters of the global models
  param = get_param_from_rdata(intensity_file) 
  
  # Extract the parameters for Other broadleaf and other conifer
  param.other = param %>%
    filter(disturbance %in% c("storm", "fire")) %>%
    # filter(species %in% c("Other broadleaf", "Other conifer")) %>%
    group_by(disturbance, species) %>%
    summarize(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c),
              dbh.intercept = mean(dbh.intercept),
              dbh.slope = mean(dbh.slope),
              logratio.intercept = mean(logratio.intercept),
              logratio.slope = mean(logratio.slope))
  
  
  
  # Fit a model with coefficients as a function of bark thickness
  # - Prepare data
  data.bt.fit = param %>%
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    filter(disturbance == "fire") %>%
    ungroup() %>%
    group_by(species)  %>%
    summarise(a0 = mean(a0), b = mean(b), c = mean(c)) %>%
    left_join(traits_compiled$traits_imputed$ShadeDrought %>%
                dplyr::select("species", "BT" = "bark.thickness"), 
              by = "species") %>%
    drop_na()
  # - Fit model
  mod.fire = lm(cbind(a0, b, c) ~ BT, data = data.bt.fit)
  # - Scale parameters
  scale.fire = param %>% ungroup() %>%
    filter(disturbance == "fire") %>%
    dplyr::select(logratio.intercept, logratio.slope, dbh.intercept, dbh.slope) %>%
    distinct()
  
  # Fit a model with coefficients as a function of bark thickness
  # - Prepare data
  data.wd.fit.storm = param %>%
    mutate(species = gsub("\\ ", "\\_", species)) %>%
    filter(disturbance == "storm") %>%
    ungroup() %>%
    group_by(species)  %>%
    summarise(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c)) %>%
    left_join(traits_compiled$traits_imputed$GrSurv %>%
                dplyr::select("species", "WD" = "wood.density"), 
              by = "species") %>%
    drop_na()
  # - Fit model
  mod.storm = lm(cbind(a0, a1, b, c) ~ WD, data = data.wd.fit.storm)
  # - Scale parameters
  scale.storm = param %>% ungroup() %>%
    filter(disturbance == "storm") %>%
    dplyr::select(logratio.intercept, logratio.slope, dbh.intercept, dbh.slope) %>%
    distinct()
  
  # Slightly change disturb_coef to include the right name for Betula
  disturb_coef.2 = matreex::disturb_coef %>%
    filter(!(species %in% c("Betula_pubescens", "Betula_sp"))) %>%
    mutate(species = ifelse(species == "Betula_pendula", "Betula", species))
  # Identify species for which we already have fire or storm parameters
  sp.fire = unique(subset(disturb_coef.2, disturbance == "fire")$species)
  sp.storm = unique(subset(disturb_coef.2, disturbance == "storm")$species)
  # Identify species not in storm and fire
  sp.nostorm = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.storm)))]
  sp.nofire = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.fire)))]
  # Predict parameters based on traits
  out = disturb_coef.2 %>%
    # Add missing storm
    rbind((data.frame(disturbance = "storm", species = sp.nostorm) %>%
             left_join(traits_compiled$traits_imputed$GrSurv %>%
                         dplyr::select("species", "WD" = "wood.density"), 
                       by = "species") %>%
             cbind(predict(mod.storm, newdata = .)) %>%
             drop_na() %>%
             mutate(dbh.intercept = scale.storm$dbh.intercept, 
                    dbh.slope = scale.storm$dbh.slope, 
                    logratio.intercept = scale.storm$logratio.intercept, 
                    logratio.slope = scale.storm$logratio.slope) %>%
             dplyr::select(-WD))) %>%
    # Add missing fire
    rbind((data.frame(disturbance = "fire", species = sp.nofire, a1 = 0) %>%
             left_join(traits_compiled$traits_imputed$ShadeDrought %>%
                         dplyr::select("species", "BT" = "bark.thickness"), 
                       by = "species") %>%
             cbind(predict(mod.fire, newdata = .)) %>%
             drop_na() %>%
             mutate(dbh.intercept = scale.storm$dbh.intercept, 
                    dbh.slope = scale.storm$dbh.slope, 
                    logratio.intercept = scale.storm$logratio.intercept, 
                    logratio.slope = scale.storm$logratio.slope) %>%
             dplyr::select(disturbance, species, a0, a1, b, c, dbh.intercept, 
                           dbh.slope, logratio.intercept, logratio.slope)))
  # Update species for which we don't have parameters
  sp.fire = unique(subset(out, disturbance == "fire")$species)
  sp.storm = unique(subset(out, disturbance == "storm")$species)
  sp.nostorm = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.storm)))]
  sp.nofire = fit_species[which(!(fit_species %in% gsub("\\ ", "\\_", sp.fire)))]
  
  # Add to the final dataset using the parameters of other broadleaf / connifer
  out = out %>%
    rbind(rbind(data.frame(disturbance = "storm", sp = sp.nostorm), 
                data.frame(disturbance = "fire", sp = sp.nofire)) %>%
            mutate(genus = gsub("\\_.+", "", sp)) %>%
            mutate(species = ifelse(genus %in% c(
              "Acer", "Alnus", "Betula", "Carpinus", "Fagus", "Fraxinus", 
              "Populus", "Prunus", "Quercus", "Salix"), 
              "Other broadleaf", "Other conifer")) %>%
            left_join(param.other, by = c("species", "disturbance")) %>%
            dplyr::select(disturbance, species = sp, a0, a1, b, c, dbh.intercept, dbh.slope, 
                          logratio.intercept, logratio.slope))
  
  # Return output
  return(out)
  
}
