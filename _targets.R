#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE, Anne Baranger, Maxime Jaunatre, Georges Kunstler
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
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "data.table", "stringr",
                 "factoextra", "modi", "sf", "rnaturalearth", "scales", 
                 "cowplot", "multcomp", "future", "GGally", #"FD",  "piecewiseSEM",
                 "statmod", "xtable", "car", "grid", "gridExtra") #, "semEff"
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
  
  # get parameters names
  tar_target(pars_list,
             get_list_pars(fit.list.allspecies)),
  
  # Mention species to exclude
  tar_target(species.excl.ipm,c("Carpinus_betulus","Juniperus_thurifera" ,"Quercus_ilex", "Salix_caprea")),
  
  # Get all species for which IPM are available, and Fundiv data
  tar_target(species.list.ipm,names(fit.list.allspecies)[!names(fit.list.allspecies) %in%
                                                           species.excl.ipm]),
  
  tar_target(sp_id,seq_along(species.list.ipm)),
  # Get all species for which disturbance pars are available
  tar_target(species.list.disturbance,species.list.ipm[species.list.ipm %in%
                                                         disturb_coef.in[disturb_coef.in$disturbance=="storm","species"][[1]]]),
  
  # Generate some harvest rules = default? ask Maskimus
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1)),
  
  # Data.frame containing storm disturbance to apply
  tar_target(disturbance.df_storm, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))),

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare data for simulations -
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' for each species define climatic combination of sggd/wai

  tar_target(climate.cat,
             make_climate_cat_pca(FUNDIV_data,
                                  species.list.ipm,
                                  n_cat=10)
  ),
  
  # For each cliamte condition, species combinations
  tar_target(species.combination,
             make_species_combinations(FUNDIV_data=FUNDIV_data,
                                       FUNDIV_plotcat=climate.cat$FUNDIV_plotcat,
                                       condi.init=climate.cat$species.cat,
                                       sp_id=sp_id,
                                       species.list.ipm=species.list.ipm, 
                                       nsp_per_richness=10,
                                       prop_threshold=0.8),
             pattern=map(sp_id),iteration = "vector"),
 
  
  # create a list of all species to be computed
  tar_target(species_list,
             make_species_list(species.combination)),

  # subselect species and climate to fit
  tar_target(species.select,
             c("Abies alba","Fagus sylvatica")),# gsub("_"," ",species.list.ipm)),
  tar_target(clim.select,
             1:10), 
  tar_target(species_list.select,
             species_list |> 
               filter(species%in%species.select) |> 
               filter(ID.spclim%in% clim.select) |> 
               mutate(species_combination=factor(species_combination),
                      id.species.obj=row_number(),
                      id.species.mu.obj=as.numeric(species_combination))),
  tar_target(species_list.mu.select,
             levels(species_list.select$species_combination)),
  tar_target(species.combination.select,
             species.combination |> 
               filter(species%in%species.select) |> 
               filter(ID.spclim%in% clim.select)),
  
  # tar_target(species.obj.id, # get Species ID 
  #            species_list.select |> 
  #              pull(id.species.obj)),#species_list$id.species.obj),
  tar_target(species.obj.mu.id,
             seq_along(species_list.mu.select)),

  # create species object, run in parallel (bottleneck step)
  # tar_target(species_object,
  #            make_species_rds(fit.list.allspecies,
  #                             species_list.select,
  #                             species.obj.id),
  #            pattern=map(species.obj.id),
  #            iteration="vector",
  #            format="file"),
  tar_target(species_object_mu,
             make_species_mu(fit.list.allspecies,
                             species_list.mu.select,
                             species.obj.mu.id),
             pattern=map(species.obj.mu.id),
             iteration="vector",
             format="file"),

  # create forest id 
  tar_target(sim_forest_list,
             create_simulation_equil_list(species.combination.select)),


  #%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make mean simulations -
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # simulation until equilibrium
  tar_target(sim_equil.id,
             sim_forest_list$id.simul_eq),#1:dim(species.combination)[1]),
  tar_target(sim_equil,
             make_simulations_equilibrium(sim_forest_list$list.forests,
                                          species_list.select,
                                          species_object_mu,
                                          harv_rules.ref,
                                          sim.type="mu",
                                          id_forest=sim_equil.id),
             pattern=map(sim_equil.id),
             iteration="vector",
             format="file"),

  # simulation for invasion
  tar_target(sim_invasion.id,
             sim_forest_list$id.simul_forest),
  tar_target(sim_invasion,
             make_simulations_invasion(sim_forest_list$list.forests,
                                       species_list.select,
                                       species_object_mu,
                                       harv_rules.ref,
                                       sim_equil,
                                       threshold_pop=200,
                                       id_forest=sim_invasion.id),
             pattern=map(sim_invasion.id),
             iteration="vector",
             format="file"),
  
  # simulation for disturbance
  tar_target(sim_dist.id,
             create_simulation_dist_list(sim_forest_list$list.forests,
                                         species.list.disturbance)),
  tar_target(sim_disturbance,
             make_simulations_disturbance(sim_forest_list$list.forests,
                                          species_list.select,
                                          species_object_mu,
                                          harv_rules.ref,
                                          sim_equil,
                                          disturb_coef.in,
                                          disturbance.df_storm,
                                          id_forest=sim_dist.id),
             pattern=map(sim_dist.id),
             iteration="vector",
             format="file"),
  
  # simulation for ibm
  tar_target(sim_ibm.id,
             sim_forest_list$id.simul_forest),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make sensitivity simulations -
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # create objects
  
  ## parameters selection for species targetted
  tar_target(pars_select,
             unlist(lapply(gsub(" ","_",species.select),
                           function(x){pars_list[grepl(x,pars_list)]})
             )),
  ## create each species object for species targetted in sensitivity analysis
  tar_target(species_object_mu_elast,
             make_species_mu_elast(fit.list.allspecies,
                                   delta=0.01,
                                   pars_select=pars_select),
             pattern=map(pars_select),
             iteration="vector",
             format="file"),
  ## create species list for species selected
  tar_target(species_list_select_elast,
             make_species_select_list_elast(species.combination.select,
                                            pars_list)),
  ## create forest list 
  tar_target(sim_forest_list_elast,
             create_simulation_equil_list_elast(species.combination.select,
                                                sim_forest_list$list.forests,
                                                pars_select)),
  
  # simulation for sensitivity
  
  # # equil
  tar_target(sim_equil_elast.id,
             sim_forest_list_elast$id.simul_forest),
  tar_target(sim_equil_elast,
             make_simulations_equilibrium_elast(sim_forest_list_elast$list.forests,
                                                species_list_select_elast,
                                                species_object_mu,
                                                species_object_mu_elast,
                                                harv_rules.ref,
                                                sim.type="mu",
                                                id_forest=sim_equil_elast.id),
             pattern=map(sim_equil_elast.id),
             iteration="vector",
             format="file"),
  
  ## invasion
  tar_target(sim_invasion_elast.id,
             sim_forest_list_elast$id.simul_forest),
  tar_target(sim_invasion_elast,
             make_simulations_invasion_elast(sim_forest_list_elast,
                                             species_list_select_elast,
                                             species_object_mu,
                                             species_object_mu_elast,
                                             harv_rules.ref,
                                             sim_equil,
                                             sim_equil_elast,
                                             threshold_pop=200,
                                             id_forest=sim_invasion_elast.id),
             pattern=map(sim_invasion_elast.id),
             iteration="vector",
             format="file"),
  
  ## disturbance
  
  #%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make mu matrix sim -
  #%%%%%%%%%%%%%%%%%%%%%%%%
  
  # tar_target(ID.model,
  #            1:4),
  # tar_target(species_clim,
  #            make_species_mu(fit.list.allspecies=fit.list.allspecies,
  #                            species.list.ipm=species.list.ipm, 
  #                            sp_id=sp_id),
  #            pattern=map(sp_id),
  #            iteration="vector",
  #            format="file"),

  #%%%%%%%%%%%%%%%%%%%
  # -- Make analysis -
  #%%%%%%%%%%%%%%%%%%%
  
  #' 1. Compute invasion rate
  tar_target(invasion_metric,
             get_invasion_rate(species.combination=sim_forest_list$list.forests,
                              id.simul_forest=sim_forest_list$id.simul_forest,
                              sim_invasion=sim_invasion,
                              fit.list.allspecies=fit.list.allspecies)),
  #' 2. Extract resilience metric
    tar_target(disturbance_metric,
             get_resilience_metrics(species.combination=sim_forest_list$list.forests,
                                    sim_dist.id,
                                    sim_disturbance,
                                    fit.list.allspecies,
                                    disturbance.df_storm)),
  
  #' 2. extract parameters of each climatic conditions (earlier?)
  
  #' 3. fit gam for each metrics
  
  #' 4. randomisation of contribution and demographic compensation?
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare data for simulations -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  # -- generate one climate object per iteration with branching
  # tar_target(climate_storm, make_climate(
  #   FUNDIV_climate_species, quantiles.in = climate_list_storm[[ID.climate_storm]], 
  #   "storm", 10, exclude.in = c("Carpinus_betulus", "Quercus_ilex", "Salix_caprea"), 
  #   method = "frequency", disturb_coef.in, traits), 
  #   pattern = map(ID.climate_storm), iteration = "list"), 
  # 
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare traits data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # -- Traits data file
  tar_target(traits_file, "data/traits.csv", format = "file"),
  # -- Read traits data
  tar_target(traits, fread(traits_file)),
  NULL
)