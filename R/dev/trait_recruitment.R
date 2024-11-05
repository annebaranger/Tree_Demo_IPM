#' Function to extract recruitment traits from demographic parameters
#' @param fit.list.allspecies demographic parameters of all species
#' @param FUNDIV_data pre-formatted tree data from FUNDIV
#' @param comp.ref character indicating how to calculate competition: 
#'                 "same" = mean competition across all species in dataset
#'                 "specific" = mean competition per species
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

# save(traits_rec,file="trait_rec.RData")
load("inv_30.rdata")
inv_30 |> left_join(traits_rec) |> 
  ggplot(aes(inv_30,recruitment))+
  geom_point()+
  geom_text(aes(label=species))
