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
#' @author Julien Barrere
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
#' @author Julien Barrere

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
#' @author Julien Barrere

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



#' Function to generate a list with climate per species
#' @param FUNDIV_data whole dataset
#' @param species.list.ipm list of species
#' @param n_cat number of categories
make_climate_cat <- function(FUNDIV_data,
                             species.list.ipm,
                             n_cat=3,
                             disturbance.in="storm"){
  FUNDIV_plotcat=FUNDIV_data |>
    # filter species with ipm
    filter(species%in%gsub("_"," ",species.list.ipm)) |>
    dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
                  species,BAtot) |> 
    group_by(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
             species) |> 
    # maximum Ba per plot
    summarize(BA=max(BAtot)) |> 
    distinct() |> 
    # create sgdd/wai category, id and limits
    group_by(species) |> 
    mutate(wai_cat=cut(wai, breaks = n_cat),
           sgdd_cat=cut(sgdd,breaks=n_cat),
           wai_id=cut(wai, breaks = n_cat,labels = 1:n_cat),
           sgdd_id=cut(sgdd, breaks = n_cat,labels = 1:n_cat),
           BA=max(BA,na.rm = TRUE)) |> 
    tidyr::separate_wider_delim(cols="wai_cat",names=c("wai_low","wai_up"),delim=",") |> 
    tidyr::separate_wider_delim(cols="sgdd_cat",names=c("sgdd_low","sgdd_up"),delim=",") |> 
    mutate(across(matches(c("low")),
                  ~as.numeric(stringr::str_sub(.,2,-1))),
           across(matches(c("up")),
                  ~as.numeric(stringr::str_sub(.,1,-2)))) |> 
    ungroup()
  
  # summarise for each category and species combinations
  condi.init<-FUNDIV_plotcat |> 
    group_by(species,BA,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up) |> 
    summarize(n_plot=n(),
              wai=mean(wai,na.rm=TRUE),
              sgdd=mean(sgdd,na.rm=TRUE)) |>
    # filter out category with too few plots
    filter(n_plot>10) |> 
    # create IPM vars
    mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
           waib = 1/(1 + wai), PC1=0, PC2=0, N = 2, SDM = 0) |> 
    ungroup() |> 
    group_by(species) |> 
    mutate(ID.spclim = row_number(),
           ID.species=cur_group_id()) |> 
    ungroup() |> 
    mutate(ID.model=row_number(),
           clim_lab=paste0("wai", wai_id, "sgdd", sgdd_id),
           file = paste0("rds/", disturbance.in, "/", 
                    species,"climate_",ID.spclim,"/species/",species,  ".rds")) 
  
  return(list(FUNDIV_plotcat=FUNDIV_plotcat,
              species.cat=condi.init))
}


#' Function to major species combinations per species and climate category
#' @param FUNDIV_data whole dataset
#' @param FUNDIV_plotcat dataset with plots per species, and corresponding climate category
#' @param species focal species 
#' @param species.list.ipm list of species
#' @param nsp_per_richness maximum of species richness per combination
#' @param prop_threshold threshold of cumulative proportion
make_species_combinations <- function(FUNDIV_data,
                                      FUNDIV_plotcat,
                                      condi.init,
                                      sp_id,
                                      species.list.ipm, 
                                      nsp_per_richness=10,
                                      prop_threshold=0.8){
  s_p=species.list.ipm[sp_id]
  sp=gsub("_"," ",s_p)
  print(sp)
  species.combinations <- FUNDIV_data |>
    left_join(FUNDIV_plotcat |>
                filter(species==sp) |>
                dplyr::select(plotcode,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up),
              by="plotcode") |>
    dplyr::select(treecode,plotcode,species,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up) |>
    filter(!is.na(wai_id)) %>%
    filter(species%in%gsub("_"," ",species.list.ipm)) |> 
    # Group by the climatic condition and plot
    group_by(wai_id, sgdd_id, plotcode) %>%
    # Summarize the species combination in each group, and species richness
    summarise(species_combination = paste(sort(unique(gsub(" ","_",species))), collapse = "."),
              n_species = n_distinct(species), .groups = 'drop') %>%
    # Now count each unique combination's frequency within each wai_cat and sgdd_cat group
    count(wai_id, sgdd_id, species_combination,n_species) |>
    filter(n_species<nsp_per_richness) |>
    group_by(wai_id, sgdd_id,) |>
    mutate(prop=n/sum(n)) |>
    # Optionally, arrange the results for better readability
    arrange(wai_id, sgdd_id, desc(n)) |>
    filter(species_combination!=s_p) |>
    mutate(prop_cum=cumsum(prop)) |>
    # select only frequent combinations (on 3 criteria)
    filter(prop_cum<prop_threshold,
           n>10,
           prop>0.02)

  species.clim.combi <- condi.init |>
    filter(species==sp) |>
    left_join(species.combinations,by=c("wai_id","sgdd_id")) |>
    arrange(wai_id,sgdd_id)
  # species.clim.combi="ok"
  return(species.clim.combi)
}

#' Function to make a species object, save it as rds and return filename
#' @param fit.list.allspecies demographic parameter for all species
#' @param condi.init table with climate condition, species and ID
#' @param ID.model n
#' @author Julien Barrere & Anne Baranger
make_species_rds <- function(fit.list.allspecies, condi.init, ID.model){
  
  climate<-condi.init[ID.model,c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", 
                                 "PC1", "PC2", "N", "SDM")]
  species<-condi.init$species[ID.model]
  s_p<-gsub(" ","_",species)

  # Make IPM
  IPM.in = make_IPM(
    species = s_p, 
    climate = climate, 
    fit =  fit.list.allspecies[[s_p]],
    clim_lab = condi.init$clim_lab[ID.model],
    mesh = c(m = 700, L = 100, U = as.numeric(
      fit.list.allspecies[[s_p]]$info[["max_dbh"]]) * 1.1),
    BA = 0:100, verbose = TRUE, correction = "none"
  )
  
  # Create species object 
  species.in = species(
    IPM.in, init_pop = def_initBA(20), harvest_fun = def_harv, disturb_fun = def_disturb)
  
  # Save species object in a rdata
  create_dir_if_needed(condi.init$file[ID.model])
  saveRDS(species.in, condi.init$file[ID.model])
  
  # Return output list
  return(condi.init$file[ID.model])
  
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

