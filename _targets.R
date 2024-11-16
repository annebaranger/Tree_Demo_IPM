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
source("R/functions_analysis.R")
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "matreex", "tidyr", "data.table", "stringr",
                 "factoextra", "modi", "sf", "rnaturalearth", "scales", 
                 "cowplot", "multcomp", "future", "GGally", #"FD",  "piecewiseSEM",
                 "statmod", "xtable", "car","mice", "grid", "gridExtra","glmmTMB","lme4") #, "semEff"
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
               # error = "null")
future::plan(future::multisession, workers = 55)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Disturbance coefficients
  tar_target(disturb_coef.raw,read.csv("data/disturb_coef_complete.csv")),
  tar_target(species_meta,get_species_meta(species.list.ipm)),  
  tar_target(disturb_coef.in, get_disturb_coef(disturb_coef.raw,
                                               species_meta,
                                               species.list.ipm)),
  
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
             unique(climate_species$sp[climate_species$sp %in% gsub(" ","_",unique(FUNDIV_data$species))])),
             # colnames(FUNDIV_climate_species)[grep("_", colnames(FUNDIV_climate_species))]), 
  
  # Get demographic parameters for all species
  tar_target(fit.list.allspecies, readRDS("data/new_fit_list.rds")),#load_param_demo(all.species.name)),
  
  # get parameters names
  tar_target(pars_list,
             get_list_pars(fit.list.allspecies)),
  
  # Mention species to exclude
  tar_target(species.excl.ipm,c("Juniperus_thurifera" , "Salix_caprea")), #"Carpinus_betulus","Quercus_ilex",
  
  # Get all species for which IPM are available, and Fundiv data
  tar_target(species.list.ipm,names(fit.list.allspecies)[!names(fit.list.allspecies) %in%
                                                           species.excl.ipm]),
  
  tar_target(sp_id,seq_along(species.list.disturbance)),
  # Get all species for which disturbance pars are available
  tar_target(species.list.disturbance,species.list.ipm[species.list.ipm %in%
                                                         disturb_coef.in[disturb_coef.in$real.coef,"species"]]),
  
  # Generate some harvest rules = default? ask Maskimus
  tar_target(harv_rules.ref, c(Pmax = 0.25, dBAmin = 3, freq = 1, alpha = 1)),
  
  # Data.frame containing storm disturbance to apply
  tar_target(disturbance.df_storm, data.frame(
    type = rep("storm", 3), intensity = rep(0.5, 3), 
    IsSurv = rep(FALSE, 3), t = c(500:502))),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Prepare data for simulations ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #' for each species define climatic combination of sggd/wai
  
  tar_target(climate.cat,
             make_climate_cat_pca(FUNDIV_data,
                                  species.list.ipm,
                                  n_cat=10)
  ),
  
  # For each cliamte condition, species combinations
  tar_target(species_competitors,
             get_species_competitors(FUNDIV_data,
                                     species.list.ipm,
                                     sp_id),
             pattern=map(sp_id),
             iteration="vector"),
  tar_target(species.combination,
             make_species_combinations(FUNDIV_data=FUNDIV_data,
                                       FUNDIV_plotcat=climate.cat$FUNDIV_plotcat,
                                       condi.init=climate.cat$species.cat,
                                       sp_id=sp_id,
                                       species.list.disturbance=species.list.disturbance,
                                       species.list.ipm=species.list.ipm, 
                                       nsp_per_richness=10,
                                       prop_threshold=0.8),
             pattern=map(sp_id),iteration = "vector"),
  
  tar_target(species.combination_2,
             make_species_combinations_2(FUNDIV_data=FUNDIV_data,
                                       FUNDIV_plotcat=climate.cat$FUNDIV_plotcat,
                                       condi.init=climate.cat$species.cat,
                                       sp_id=sp_id,
                                       species.list.disturbance=species.list.disturbance,
                                       species.list.ipm=species.list.ipm, 
                                       nsp_per_richness=10,
                                       prop_threshold=0.8),
             pattern=map(sp_id),iteration = "vector"),
  
  tar_target(clim_bound,
             make_clim_boundaries(species_list.select)),
  # create a list of all species to be computed
  tar_target(species_list,
             make_species_list(species.combination_2)),
  
  # subselect species and climate to fit
  tar_target(species.select,
             gsub("_"," ",species.list.disturbance)),
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
             species.combination_2 |> 
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
                             clim_bound,
                             species.obj.mu.id),
             pattern=map(species.obj.mu.id),
             iteration="vector",
             format="file"),
  
  # create forest id 
  tar_target(sim_forest_list,
             create_simulation_equil_list(species.combination.select)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make mean simulations ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
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
  # Identify simul with exclusion
  tar_target(sim_forest_excl,
             is_species_excluded(sim_forest_list,
                                 sim_equil)
  ),
  
  
  # simulation for invasion
  tar_target(sim_invasion.id,
             sim_forest_list$id.simul_forest),
  tar_target(sim_invasion,
             make_simulations_invasion_2(sim_forest_list$list.forests,
                                         species_list.select,
                                         species_object_mu,
                                         harv_rules.ref,
                                         sim_equil,
                                         BA_target=1,
                                         threshold_pop=200,
                                         id_forest=sim_invasion.id),
             pattern=map(sim_invasion.id),
             iteration="vector",
             format="file"),
  
  tar_target(sim_invasion_2,
             make_simulations_invasion_3(sim_forest_list$list.forests,
                                         species_list.select,
                                         species_object_mu,
                                         harv_rules.ref,
                                         sim_equil,
                                         BA_target=1,
                                         threshold_pop=200,
                                         id_forest=sim_invasion.id),
             pattern=map(sim_invasion.id),
             iteration="vector",
             format="file"),
  
  # simulation for disturbance
  tar_target(sim_dist.id,
             sim_forest_list$id.simul_forest),
             # create_simulation_dist_list(sim_forest_list$list.forests,
             #                             species.list.disturbance)),
  tar_target(sim_disturbance,
             make_simulations_disturbance(sim_forest_excl$list.forests,
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
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Make sensitivity simulations ----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
  tar_target(invasion_metric_elast,
             get_invasion_rate_elast(species.combination=sim_forest_list_elast$list.forests,
                                     id.simul_forest=sim_invasion_elast.id,
                                     sim_invasion=sim_invasion_elast,
                                     fit.list.allspecies=fit.list.allspecies)),
  ## disturbance
  tar_target(sim_dist_elast.id,
             create_simulation_dist_list_elast(sim_forest_list_elast$list.forests,
                                               species.list.disturbance)),
  tar_target(sim_disturbance_elast,
             make_simulations_disturbance_elast(sim_forest_list_elast,
                                                species_list_select_elast,
                                                species_object_mu_elast,
                                                species_object_mu,
                                                harv_rules.ref,
                                                sim_equil_elast,
                                                disturb_coef.in,
                                                disturbance.df_storm,
                                                fit.list.allspecies,
                                                id_forest=sim_dist_elast.id),
             pattern=map(sim_dist_elast.id),
             iteration="vector"),
  tar_target(disturbance_metric_elast,
             bind_rows(sim_disturbance_elast[grep( pattern = "_out", 
                                                   x = names(sim_disturbance_elast))])),
  
  
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
  
  #%%%%%%%%%%%%%%%%%%%%%%
  # -- Make analysis ----
  #%%%%%%%%%%%%%%%%%%%%%%
  
  #' 0. Mean demography of species
  tar_target(mean_demo,
             make_mean_simul(species.list.ipm,
                             climate.cat,
                             species_object=species_object_mu,
                             species_list=species_list.select,
                             fit.list.allspecies,
                             harv_rules.ref,
                             sp_id=species.obj.mu.id),
             pattern=map(species.obj.mu.id),
             iteration="vector"),
  
  #' 1. Compute invasion rate
  tar_target(invasion_metric,
             get_invasion_rate(species.combination=sim_forest_excl$list.forests,
                               id.simul_forest=sim_forest_list$id.simul_forest,
                               sim_invasion=sim_invasion,
                               fit.list.allspecies=fit.list.allspecies)),
  tar_target(invasion_metric_2,
             get_invasion_rate_2(species.combination=sim_forest_excl$list.forests,
                               id.simul_forest=sim_forest_list$id.simul_forest,
                               sim_invasion=sim_invasion,
                               fit.list.allspecies=fit.list.allspecies)),
  tar_target(invasion_metric_3,
             get_invasion_rate_3(species.combination=sim_forest_excl$list.forests,
                                 id.simul_forest=sim_forest_list$id.simul_forest,
                                 sim_invasion=sim_invasion,
                                 fit.list.allspecies=fit.list.allspecies)),
  #' 2. Extract resilience metric
  tar_target(
    forest_iteration,
    sim_dist.id
  ),
  tar_target(disturbance_metric,
             get_resilience_metrics(species.combination=sim_forest_excl$list.forests,
                                    sim_dist.id,
                                    forest_iteration,
                                    sim_disturbance,
                                    fit.list.allspecies,
                                    disturbance.df_storm),
             pattern=map(sim_dist.id),
             iteration = "vector"),
  
  #' 3. Get BA initial
  tar_target(invasion_ba,
             get_bainit_inv(sim_forest_excl,
                            sim_equil,
                            invasion_metric,
                            elast=FALSE,
                            species.list.ipm)),
  tar_target(invasion_ba_2,
             get_bainit_inv(sim_forest_excl,
                            sim_equil,
                            invasion_metric_2,
                            elast=FALSE,
                            species.list.ipm)),
  tar_target(invasion_ba_elast,
             get_bainit_inv(sim_forest_list_elast,
                            sim_equil,
                            invasion_metric_elast,
                            elast=TRUE,
                            species.list.ipm)),
  tar_target(disturbance_ba,
             get_bainit_dist(sim_forest_excl,
                             sim_equil,
                             disturbance_metric,
                             elast=FALSE,
                             species.list.ipm)),
  tar_target(disturbance_ba_elast,
             get_bainit_dist(sim_forest_list_elast,
                             sim_equil_elast,
                             disturbance_metric_elast,
                             elast=TRUE,
                             species.list.ipm)),
  
  #' 4. Get BA final
  tar_target(ba_dif,
             get_dif_ba(disturbance_ba,
                        species.list.ipm)),
  
  
  #'5 Performance final
  tar_target(performance,
             get_performance(species.list.ipm,
                             disturbance_ba,
                             invasion_ba_2,
                             ba_dif,
                             sim_forest_excl,
                             climate.cat,
                             trait_complete,
                             trait_selection=c("HM","recruitment","WD"))),
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
  # -- Traits shade
  tar_target(file.shade,"data/Traits/data_Niinemets&Valladares_2006.csv",format="file"),
  tar_target(trait_shade,
             get_shade(file.shade)),
  # -- Traits Julien (Height/DBH, maxgrowth, WD)
  tar_target(file.traitsJul, "data/Traits/traits_JB_FunEcol.csv", format = "file"),
  tar_target(trait_jul, fread(file.traitsJul)),
  # -- TRY
  tar_target(file.try,"data/Traits/try/try_query.csv",format="file"),
  tar_target(trait_try,get_try(file.try)),
  
  # -- Wood Density
  tar_target(file.wd,"data/Traits/GlobalWoodDensityDatabase.csv",format="file"),
  tar_target(trait_wd,get_WD(file.wd,trait_try)),
  
  # -- Max Heigth
  tar_target(trait_maxH,get_maxH(FUNDIV_data)),
  
  # -- Recruitment
  tar_target(trait_recruitment,
             get_recruitment_traits(fit.list.allspecies, 
                                   FUNDIV_data,
                                   comp.ref="specific")),
  
  # -- Trait tot
  tar_target(trait_raw,
             get_traits(trait_shade,
                        trait_jul,
                        trait_try,
                        trait_wd,
                        trait_maxH,
                        trait_recruitment,
                        species.list.ipm)),
  tar_target(trait_complete,
             complete_trait(trait_raw)),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Run models -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  tar_target(models_traits,
             run_model(performance,
                       climate.cat,
                       trait_complete)),
  NULL
)