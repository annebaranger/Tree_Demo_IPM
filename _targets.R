#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE, Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
# lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
source("R/functions_data.R")
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "data.table", 
                 "factoextra", "modi", "sf", "rnaturalearth", "scales", 
                 "cowplot", "multcomp", "piecewiseSEM", "future", "FD", "GGally", 
                 "statmod", "xtable", "car", "modi", "grid", "gridExtra", "semEff")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Disturbance coefficients
  tar_target(disturb_coef.in, fread("data/disturb_coef.csv")),
  
  # Raw data from FUNDIV
  tar_target(FUNDIV_tree_file, "data/FunDiv_trees_Nadja.csv", format = "file"),
  tar_target(FUNDIV_plot_file, "data/FunDiv_plots_Nadja.csv", format = "file"),
  tar_target(FUNDIV_species_file, "data/FunDiv_species_Nadja.csv", format = "file"),
  tar_target(FUNDIV_climate_file, "data/moreno_chelsa_fundiv_clim.csv", format = "file"),
  
  # Read and format FUNDIV data
  tar_target(FUNDIV_data, read_FUNDIV(FUNDIV_tree_file, 
                                      FUNDIV_plot_file, 
                                      FUNDIV_climate_file, 
                                      FUNDIV_species_file)),
  tar_target(FUNDIV_climate_species, get_FUNDIV_species_per_climate(FUNDIV_data)),
  
  # List of species for which we have data
  tar_target(all.species.name, 
             colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]), 
  
  # Get demographic parameters for all species
  tar_target(fit.list.allspecies, load_param_demo(all.species.name)),
  
  # Get all species for which IPM are available, and Fundiv data
  tar_target(species.list.ipm,names(fit.list.allspecies)),
  
  # Mention species to exclude
  tar_target(species.excl.ipm,c("Carpinus_betulus", "Quercus_ilex", "Salix_caprea")),
  
  # Get all species for which disturbance pars are available
  tar_target(species.list.disturbance,species.list.ipm[species.list.ipm %in%
                                                         disturb_coef.in[disturb_coef.in$disturbance=="storm","species"][[1]]]),
  
  # Generate some harvest rules = default? ask Maskimus
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 5, alpha = 1)),
  
  # Data.frame containing storm disturbance to apply
  tar_target(disturbance.df_storm, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))),

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare data for simulations - Monospecific ---
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' for each species define climatic combination of sggd/wai
  
  tar_target(climate.species.condition,
             {FUNDIV_data |> 
                 filter()})


  #' 2. create a vector along which iterate (i.e. number of species)

  #' 3. for each species, create all climate possible
  
  #' 4. create a species list
  
  #' 5. vector of ID per climate /species ?
  
  #' 6. create directory of species object

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make simulations - Monospecific ---
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' 1. vector of sim
  
  #' 2. simulation of equilibrium
  
  #' 3. simulation with perturbation
  
  #' 4. simulation starting from 1m2 (but equil distrib)
  
  #' 5. IBM
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make analysis - Monospecific ---
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' 1. Extract resilience metric
  
  #' 2. extract parameters of each climatic conditions (earlier?)
  
  #' 3. fit gam for each metrics
  
  #' 4. randomisation of contribution and demographic compensation?
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare data for simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Generate some climates
  # -- iterations along all climates that will be created (one iteration per climate)
  tar_target(ID.climate_storm, c(1:10)),
  # -- list of climates
  tar_target(climate_list_storm, create_climate_list(length(ID.climate_storm), 
                                                     quantile.range = c(0.2, 1))),
  # -- generate one climate object per iteration with branching
  tar_target(climate_storm, make_climate(
    FUNDIV_climate_species, quantiles.in = climate_list_storm[[ID.climate_storm]], 
    "storm", 10, exclude.in = c("Carpinus_betulus", "Quercus_ilex", "Salix_caprea"), 
    method = "frequency", disturb_coef.in, traits), 
    pattern = map(ID.climate_storm), iteration = "list"), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare traits data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  # -- Read traits data
  tar_target(traits, fread(traits_file)),
  NULL
)