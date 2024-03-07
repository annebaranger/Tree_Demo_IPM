#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Anne Baranger, Julien Barrere
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
#' @author Julien Barrere
create_dir_if_needed <- function(file.in){
  
  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}

#' Function to convert a binary code in a vector of species
#' @param code binary code where one indicates presence, 0 absence
#' @param species_vec vector of species, same length as code
decode_species <- function(code, species_vec){
  (data.frame(present = as.numeric(strsplit(code, split = "")[[1]]), 
              species = species_vec) %>%
     filter(present == 1))$species
}

#' Function to load and pre-format data from FUNDIV files
#' @param FUNDIV_tree_file Location of the file containing FUNDIV tree data
#' @param FUNDIV_plot_file Location of the file containing FUNDIV plot data
#' @param FUNDIV_climate_file Location of the file containing FUNDIV climate data
#' @param FUNDIV_species_file Location of the file containing FUNDIV species data
#' @author Julien Barrere
read_FUNDIV = function(FUNDIV_tree_file, FUNDIV_plot_file, 
                       FUNDIV_climate_file, FUNDIV_species_file){
  
  # Read climatic data
  FUNDIV_climate = fread(FUNDIV_climate_file) %>%
    # Remove Belgium
    filter(country != "WA") %>%
    # Only keep the mean sgdd and wai
    dplyr::select(plotcode, sgdd = sgdd_moreno_m, wai = wai_moreno_m) %>%
    # Remove null sgdd (points too close from the sea)
    filter(sgdd != 0)
  
  # Read tree, species and plot data
  FUNDIV_plot = fread(FUNDIV_plot_file)
  FUNDIV_tree = fread(FUNDIV_tree_file)
  FUNDIV_species = fread(FUNDIV_species_file)
  
  # Correct the weight in FUNDIV_tree
  FUNDIV_tree = FUNDIV_tree %>%
    mutate(weight1 = case_when(country == "DE" ~ weight1, 
                               country == "FG" ~ ba_ha1/ba1,
                               country %in% c("ES", "FI", "SW", "WA") ~ 10000/(pi*weight1^2)), 
           weight2 = case_when(country == "DE" ~ weight2, 
                               country == "FG" ~ ba_ha2/ba2,
                               country %in% c("ES", "FG", "FI", "SW", "WA") ~ 10000/(pi*weight2^2))) |> 
    mutate(across(is.numeric,
                  ~na_if(as.numeric(.),Inf)))

  
  # Correct  species anomalies 
  FUNDIV_tree$speciesid[FUNDIV_tree$speciesid %in% c(46,47)] <- 48
  FUNDIV_species$species[FUNDIV_species$id ==277] <- "pubescens"
  FUNDIV_species <- FUNDIV_species |> 
    mutate(sp = paste(genus, species)) |> 
    dplyr::select(c(id, sp))
  FUNDIV_species[FUNDIV_species$id == 277, "sp"] <- "Quercus pubescens"
  FUNDIV_species[FUNDIV_species$id == 48, "sp"] <- "Betula"
  
  ## - Make pca with sgdd and wai
  pca <- prcomp((FUNDIV_climate %>% dplyr::select(sgdd, wai) %>% na.omit()), 
                center = TRUE, scale = TRUE)
  
  ## - Format dataset
  out <- FUNDIV_tree %>%
    # Remove belgium and points in the sea
    filter(plotcode %in% FUNDIV_climate$plotcode) %>%
    # Add species name
    left_join((FUNDIV_species %>% 
                 na.omit() %>%
                 dplyr::select(speciesid = id, species = sp)), 
              by = "speciesid") %>% 
    # Add plot coordinates and dates
    left_join((FUNDIV_plot %>% 
                 dplyr::select(plotcode, longitude, latitude, yearsbetweensurveys)), 
              by = "plotcode") %>%
    # Add climate
    left_join((FUNDIV_climate %>% 
                 dplyr::select(plotcode, sgdd, wai) %>%
                 na.omit() %>%
                 mutate(pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2])), 
              by = "plotcode") %>%
    # Calculate competition
    group_by(plotcode) %>% mutate(BAtot = sum(ba_ha1) - ba_ha1) %>% 
    ungroup %>% group_by(plotcode, species) %>% 
    mutate(BAtotSP = sum(ba_ha1) - ba_ha1, BAtotNONSP = BAtot - BAtotSP) %>%
    ungroup()
  
  # Return output
  return(out)
  
}




#' Function to create a dataset with only climate and presence of species
#' @param FUNDIV_data Pre-formatted FUNDIV data
get_FUNDIV_species_per_climate = function(FUNDIV_data){
  
  FUNDIV_data %>%
    # Keep columns of interest
    dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2, species) %>%
    # Keep only species in the IPM (climate_species is the df of species climate 
    # in matreex package)
    filter(species %in% gsub("\\_", "\\ ", unique(climate_species$sp))) %>%
    # Remove duplicates and na
    distinct() %>%
    na.omit() %>%
    # Transform species presence in a binary variable per column
    mutate(value = 1, species = gsub("\\ ", "\\_", species)) %>%
    tidyr::spread(key = "species", value = "value") %>% 
    replace(is.na(.), 0)
  
}



#' Function to load demographic parameters for several species
#' @param species.names character vector of species ("Genus_species" format)
load_param_demo <- function(species.names){
  
  # Initialize list of fits (one element per species)
  fit.list <- list()
  
  # Loop on all species to load and store demographic parameters
  for(i in 1:length(species.names)){
    eval(parse(text=paste0("fit.list$", species.names[i], " <- fit_", species.names[i])))
  }
  
  # Return the list
  return(fit.list)
  
}


#' Function that creates climate list based on quantiles
#' @param n.clim integer: number of climate to create in the list
#' @param quantile.range numeric vector of length two indicating the climatic range 
create_climate_list = function(n.clim, quantile.range = c(0, 1)){
  
  # Vector that contains a sequence of quantiles value from 0 to 1
  vec.in = seq(from = quantile.range[1], to = quantile.range[2], 
               length.out = n.clim+1)
  
  # Initialize the output list
  list.out = vector(mode = "list", length = n.clim)
  
  # Loop on all climates
  for(i in 1:n.clim){
    # Attribute a name to climate i
    names(list.out)[i] = paste("quantile", vec.in[i], vec.in[i+1], sep = "_")
    # Add vector of quantile for climate i
    list.out[[i]] = c(vec.in[i], vec.in[i+1])
  }
  
  # Return final list
  return(list.out)
}



#' Function to generate a list with climate with species combinations
#' @param FUNDIV_climate_species data with climate and sp presence per plot
#' @param quantiles.in range between 0 and 1 of pca1 value to select
#' @param disturbance.in name of the disturbance we plan to apply to filter 
#'                       species combinations compatible
#' @param nsp_per_richness number of sp combinations to select per sp richness
#' @param exclude.in vector of species to exclude if bad estimation or IPM fit
#' @param method way to select species combination: most frequent ("frequency") or "random" 
#' @param disturb_coef.in disturbance coefficients
make_climate <- function(FUNDIV_climate_species, quantiles.in, disturbance.in, 
                         nsp_per_richness = 10, exclude.in = c("Carpinus_betulus"), 
                         method = "frequency", disturb_coef.in, pc1_per_species){
  
  # Initialize output list
  out = list()
  
  # Get the range of pca_values based on quantiles.in
  climate.in = as.numeric(quantile(FUNDIV_climate_species$pca1, 
                                   probs = quantiles.in))
  
  # climate vector
  out$climate = setNames(
    object = as.numeric(
      FUNDIV_climate_species %>%
        filter(pca1 > climate.in[1] & pca1 < climate.in[2]) %>%
        summarize(sgdd = mean(sgdd), wai = mean(wai), PC1 = mean(pca1), PC2 = mean(pca2)) %>%
        mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
               waib = 1/(1 + wai), N = 2, SDM = 0) %>%
        dplyr::select(sgdd, wai, sgddb, waib, wai2, sgdd2, PC1, PC2, N, SDM)
    ), c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")
  )
  
  # Vector of all species
  species_vec = colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]
  # Remove species with bad estimation or unstable in simulations
  species_vec = species_vec[!(species_vec %in% exclude.in)]
  # Remove species for which we have no trait estimation
  species_vec = species_vec[(species_vec %in% pc1_per_species$species)]
  
  # Adjust species combinations to the disturbance if one is specified
  if(disturbance.in %in% c("storm", "fire", "biotic")){
    # Vector of all species for which we have disturbance parameters
    species_vec_dist = (disturb_coef.in %>%
                          filter(disturbance %in% disturbance.in))$species
    # Restrict the species vector to these species
    species_vec = species_vec[which(species_vec %in% species_vec_dist)]
  }
  
  
  
  # species combinations for this climate based on data (frequency method)
  data.in <- FUNDIV_climate_species |> 
    # -- Paste all species presence absence to create binary code
    mutate(combi = purrr::pmap_chr(across(all_of(species_vec)), paste0, collapse = "")) |> 
    # -- Calculate the number of species per species combination
    mutate(n.sp = rowSums(across(all_of(species_vec)), na.rm = TRUE)) |> 
    # -- Restrict to the climate specified
    filter(pca1 > climate.in[1], pca1 < climate.in[2])
  
  # -- Select the main combinations
  data_codes <- data.in %>%
    group_by(combi, n.sp) %>%
    summarize(n = n()) %>%
    filter(n.sp > 0) %>%
    arrange(desc(n.sp), desc(n))
  codes = c() # initialize list of species combination
  # for the climatic range considered, select the 10 most frequent species combinations
  # at different species richness thresholds
  for(j in 1:length(unique(data_codes$n.sp))){
    codes.j = (data_codes %>%
                 filter(n.sp == unique(data_codes$n.sp)[j]))$combi
    if(length(codes.j) > nsp_per_richness) codes = c(codes, codes.j[c(1:nsp_per_richness)])
    else codes = c(codes, codes.j)
  }
  # -- initialize combinations and species vector
  combinations.in = c(); species.in = c()
  # -- Loop on all codes
  for(i in 1:length(codes)){
    # Decode to get a vector of species
    vec.i = decode_species(codes[i], species_vec)
    # Add combination to the vector
    combinations.in = c(combinations.in, paste(vec.i, collapse = "."))
    # Store all species in species vector
    species.in = c(species.in, vec.i)
  }
  # -- Identify all species for this climate
  species.in = unique(species.in)
  
  
  
  # Make a random selection of species with same number of forest per richness
  # -- Identify the number of combinations per species richness
  combi.per.richness = data_codes %>%
    filter(combi %in% codes) %>%
    group_by(n.sp) %>%
    summarize(n = n()) %>%
    arrange(n.sp)
  # -- Initialize the vector of random combinations
  combinations.random.in = c()
  # -- Loop on all levels of richness
  for(r in 1:dim(combi.per.richness)[1]){
    # All existing combinations for richness r
    combinations.r.all = apply(as.data.frame(t(combn(species.in, r))), 1, 
                               paste, collapse = "." )
    # Randomly sample the same number of forest per richness as with the data approach
    combinations.r = sample(combinations.r.all, size = combi.per.richness$n[r], 
                            replace = FALSE)
    # Add to the vector of random species combinations
    combinations.random.in = c(combinations.random.in, combinations.r)
  }
  
  
  # Add to the final list the combinations and the list of all species
  if(method == "frequency") out$combinations = combinations.in
  if(method == "random") out$combinations = combinations.random.in
  out$species = species.in
  
  # Return output
  return(out)
}


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
    traits_rec$recruitment[i] = exp(sum(vec.rec.i*vec.var.i))
    
  }
  
  # Return the trait dataset generated
  return(traits_rec)
}




#' Get coordinate in first pca traits axis per species
#' @param traits dataframe containing trait value per species
#' @param traits_rec recruitment traits per species
get_pc12_per_species <- function(traits, traits_rec){
  
  # Compile the traits data
  data_traits = left_join(traits, traits_rec, by = "species") %>%
    dplyr::select(-delay) %>%
    drop_na()
  
  # Make the PCA
  pca <- prcomp((data_traits %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # Extract data for individuals
  out = data.frame(species = data_traits$species, 
                   pca1 = get_pca_ind(pca)[[1]][, 1], 
                   pca2 = get_pca_ind(pca)[[1]][, 2])
  # return the output
  return(out)
}







#' Disturbance function
#'
#' @param x population state distribution at time t
#' @param species The species class object of interest to get mesh and RDIcoef
#' values from. RDIcoef is a one line dataframe with RDI coefficient for one
#' species.
#' @param disturb Disturbance parameters. Highly depend on the disturbance
#' impact parameters given to the species.
#' @param ... Not used in this case.
#' \describe{
#' \item{qmd}{Forest Quadratic Mean Diameter}
#' }
#' @author Maxime Jeaunatre
#'
disturb_fun <- function(x, species, disturb = NULL, ...){
  
  dots <- list(...)
  qmd <- dots$qmd 
  size <- species$IPM$mesh
  coef <- species$disturb_coef
  if(any(disturb$type %in% coef$disturbance)){
    coef <- subset(coef, disturbance == disturb$type)
  } else {
    stop(sprintf("The species %s miss this disturbance type (%s) parameters",
                 sp_name(species), disturb$type))
  }
  
  # edits for delay
  size[size == 0] <- min(size[size !=0])
  
  logratio <-  log(size / qmd)
  dbh.scaled = coef$dbh.intercept + size * coef$dbh.slope
  logratio.scaled = coef$logratio.intercept + logratio * coef$logratio.slope
  Pkill <- plogis(coef$a0 + coef$a1 * logratio.scaled + 
                    coef$b * disturb$intensity ^(coef$c * dbh.scaled))
  
  return(x* Pkill) # always return the mortality distribution
}

