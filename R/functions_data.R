#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Anne Baranger, Julien Barrere, Maxime Jaunatre
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


#' Function to generate a list with climate per species
#' @param FUNDIV_data whole dataset
#' @param species.list.ipm list of species
#' @param n_cat number of categories
make_climate_cat <- function(FUNDIV_data,
                             species.list.ipm,
                             max_cat=80,
                             min_cat=10){
  # summarize each plots 
  FUNDIV_plot=FUNDIV_data |> 
    filter(species %in% gsub("_"," ",species.list.ipm)) |> 
    dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
                  species,BAtot) |> 
    group_by(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
             species) |> 
    # maximum Ba per plot
    summarize(BA=max(BAtot)) |> 
    distinct() |> 
    ungroup()
  # 
  # step_cat<- FUNDIV_plot |>
  #   pivot_longer(cols=c("wai","sgdd")) |>
  #   group_by(species,name) |>
  #   summarise(range_clim=quantile(value,probs = 0.95)-quantile(value,probs = 0.05)) |>
  #   arrange(name,range_clim) |>
  #   group_by(name) |>
  #   mutate(extrema=case_when(range_clim==min(range_clim)~"min",
  #                            range_clim==max(range_clim)~"max")) |>
  #   filter(extrema%in%c("min","max")) |>
  #   select(-species) |>
  #   pivot_wider(names_from = extrema,
  #               values_from = range_clim) |>
  #   mutate(step_min=min/6,
  #          step_max=max/15) |>
  #   mutate(step=max(step_min,step_max)) |>
  #   ungroup() |>
  #   select(name,step)
  step_cat<- data.frame(name=c("sgdd","wai"),
                        step=c(95,0.095))
  
  max_obs<-max_cat+1
  while(max_obs>max_cat){
    print(max_obs)
    step_cat$step[1]=step_cat$step[1]+5
    step_cat$step[2]=step_cat$step[2]+0.005
    
    FUNDIV_plot |> 
      pivot_longer(cols=c("wai","sgdd")) |> 
      group_by(species,name) |> 
      summarise(range_clim=quantile(value,probs = 0.95)-quantile(value,probs = 0.05)) |> 
      left_join(step_cat) |> 
      mutate(n_breaks=round(range_clim/step)) |> 
      dplyr::select(species,name, n_breaks) |> 
      pivot_wider(names_from = name,
                  values_from = n_breaks) |> 
      rename(sgdd_breaks=sgdd,wai_breaks=wai)->step_species
    
    
    
    FUNDIV_plotcat <- FUNDIV_plot |> 
      left_join(step_species) |> 
      group_by(species) |> 
      mutate(wai_cat=cut(wai, 
                         breaks = seq(from=min(wai,na.rm=TRUE),
                                      to=max(wai,na.rm=TRUE),
                                      length.out=(wai_breaks+1)),
                         include.lowest = TRUE),
             sgdd_cat=cut(sgdd,
                          breaks=seq(min(sgdd,na.rm=TRUE)
                                     ,max(sgdd,na.rm=TRUE),
                                     length.out=(sgdd_breaks+1)),
                          include.lowest = TRUE)) |> 
      tidyr::separate_wider_delim(cols="wai_cat",names=c("wai_low","wai_up"),delim=",",cols_remove = FALSE) |> 
      tidyr::separate_wider_delim(cols="sgdd_cat",names=c("sgdd_low","sgdd_up"),delim=",",cols_remove = FALSE) |> 
      mutate(across(matches(c("low")),
                    ~as.numeric(stringr::str_sub(.,2,-1))),
             across(matches(c("up")),
                    ~as.numeric(stringr::str_sub(.,1,-2)))
      ) |> 
      ungroup()  |> 
      group_by(species) %>%
      mutate(
        wai_id = as.integer(factor(wai_cat)),
        sgdd_id = as.integer(factor(sgdd_cat))
      ) %>%
      ungroup()
    
    condi.init<- FUNDIV_plotcat |> 
      group_by(species,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up) |> 
      summarize(n_plot=n(),
                wai=mean(wai,na.rm=TRUE),
                sgdd=mean(sgdd,na.rm=TRUE)) |>
      # filter out category with too few plots
      filter(n_plot>min_cat) |> 
      # create IPM vars
      mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
             waib = 1/(1 + wai), PC1=0, PC2=0, N = 2, SDM = 0) |> 
      ungroup() |> 
      group_by(species) |> 
      mutate(ID.spclim = row_number(),
             clim_lab=paste0("wai", wai_id, "sgdd", sgdd_id),
             ID.species=cur_group_id()) |> 
      ungroup()
    condi.init |> 
      group_by(species) |> 
      summarise(n=n()) |> 
      pull(n) |> 
      max() -> max_obs

  }
  

    # mutate(ID.model=row_number(),
    #        clim_lab=paste0("wai", wai_id, "sgdd", sgdd_id),
    #        file = paste0("rds/", disturbance.in, "/", 
    #                 species,"climate_",ID.spclim,"/species/",species,  ".rds")) 
  
  return(list(FUNDIV_plotcat=FUNDIV_plotcat,
              species.cat=condi.init))
}


#' Function to identify main species combinations per species and climate category
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
    filter(species%in%gsub("_"," ",species.list.ipm)) |> 
    # join with climate cat of the targetted species
    left_join(FUNDIV_plotcat |>
                filter(species==sp) |>
                dplyr::select(plotcode,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up),
              by="plotcode") |>
    dplyr::select(treecode,plotcode,species,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up) |>
    filter(!is.na(wai_id)) %>%
    # Group by the climatic condition and plot
    group_by(wai_id, sgdd_id, plotcode) %>%
    # Summarize the species combination in each group, and species richness
    summarise(species_combination = paste(sort(unique(gsub(" ","_",species))), collapse = "."),
              n_species = n_distinct(species), .groups = 'drop')  |> 
    #  count each unique combination's frequency within each wai_cat and sgdd_cat group
    count(wai_id, sgdd_id, species_combination,n_species) |>
    filter(n_species<nsp_per_richness) |>
    group_by(wai_id, sgdd_id) |>
    mutate(prop=n/sum(n)) |>
    # Optionally, arrange the results for better readability
    arrange(wai_id, sgdd_id, desc(n)) 
  
  
  species.target<- FUNDIV_plotcat |>
    filter(species==sp) |> 
    dplyr::select(wai_id,sgdd_id) |> 
    distinct() |> 
    mutate(species_combination=s_p) |> 
    left_join(species.combinations |>
                filter(species_combination==s_p),
              by=c("wai_id","sgdd_id","species_combination")) |> 
    mutate(prop_cum=NA)
  
  species.combinations.other<-species.combinations |> 
    filter(species_combination!=s_p) |>
    mutate(prop_cum=cumsum(prop)) |> 
    filter(prop_cum<prop_threshold,
           n>10,
           prop>0.02) 

  species.clim.combi <- condi.init |>
    filter(species==sp) |>
    left_join(rbind(species.target,species.combinations.other),
              by=c("wai_id","sgdd_id")) |>
    arrange(wai_id,sgdd_id)
  
  # species.clim.combi="ok"
  return(species.clim.combi)
}

#' Function create a list of all species object to be created
#' @param species.combination table of species combination for each climate cat
make_species_list <- function(species.combination){
  species_list<-species.combination |> 
    dplyr::select(species,wai_id,sgdd_id,ID.spclim,clim_lab, # climate info
           wai,sgdd,sgdd2,wai2,sgddb,waib,PC1,PC2,N,SDM, # data for IPM clim
           species_combination) |> # eponyme
    mutate(species_combination=strsplit(species_combination,"\\.")) |> 
    unnest(cols = species_combination) |> 
    unique() |> 
    mutate(file.ipm=paste0("rds/",gsub(" ","_",species),"/clim_",ID.spclim,"/",species_combination,".rds"))
  
  return(species_list)
}


#' Function to make a species object, save it as rds and return filename
#' @param fit.list.allspecies demographic parameter for all species
#' @param condi.init table with climate condition, species and ID
#' @param ID.model n
#' @author Julien Barrere & Anne Baranger
make_species_rds <- function(fit.list.allspecies, species_list, species.obj.id){
  
  climate<-species_list[species.obj.id,c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", 
                                 "PC1", "PC2", "N", "SDM")]
  sp<-species_list$species_combination[species.obj.id]
  s_p<-gsub(" ","_",sp)

  # Make IPM
  IPM.in = make_IPM(
    species = s_p, 
    climate = climate, 
    fit =  fit.list.allspecies[[s_p]],
    clim_lab = species_list$clim_lab[species.obj.id],
    mesh = c(m = 700, L = 100, U = as.numeric(
      fit.list.allspecies[[s_p]]$info[["max_dbh"]]) * 1.1),
    BA = 0:100, verbose = TRUE, correction = "none"
  )
  
  # Create species object 
  species.in = species(
    IPM.in, init_pop = def_initBA(20), harvest_fun = def_harv, disturb_fun = def_disturb)
  
  # Save species object in a rdata
  create_dir_if_needed(species_list$file.ipm[species.obj.id])
  saveRDS(species.in, species_list$file.ipm[species.obj.id])
  
  # Return output list
  return(species_list$file.ipm[species.obj.id])
  
}


#' Function to make a species mu matrix, create species object, and save it 
#' as rds and return filename
#' @param fit.list.allspecies demographic parameter for all species
#' @param ID.model n
#' @author Anne Baranger
make_species_mu <- function(fit.list.allspecies,
                            species.list.ipm, 
                            sp_id){
  s_p=species.list.ipm[sp_id]
  sp=gsub("_"," ",s_p)
  print(sp)

  # Make IPM
  IPM.mu <- make_mu_gr(
    species = s_p, fit = fit.list.allspecies[[s_p]],
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit.list.allspecies[[s_p]]) * 1.1),
    verbose = TRUE, stepMu = 0.0001)
  
    # Create species object 
  species.in = species(IPM=IPM.mu,
                       init_pop = def_initBA(20),
                       harvest_fun = def_harv, 
                       disturb_fun = def_disturb)
  mu.file=paste0("rds/",s_p,"/",s_p,"_mu.rds")
  species.file=paste0("rds/",s_p,"/",s_p,"_species.rds")
  
  
  # Save species object in a rdata
  create_dir_if_needed(mu.file)
  saveRDS(IPM.mu, mu.file)
  create_dir_if_needed(species.file)
  saveRDS(species.in, species.file)
  
  # Return output list
  return(mu.file)
  
}

#' Function to make a list of forest and mean associated climate for mean IPMs
#' @param FUNDIV_data 
#' @param sim_forest_list 
#' @author Anne Baranger
make_mean_forest_list <- function(FUNDIV_data,
                                  sim_forest_list=sim_forest_list[["list.forests"]]){
  list.forest <- sim_forest_list |> 
    pull(species_combination) |> 
    unique() 
  
  
  mean_forest_list <-FUNDIV_data |> 
    group_by(plotcode,wai,sgdd) |> 
    summarise(species_combination = paste(sort(unique(gsub(" ","_",species))), collapse = "."),
              n_species = n_distinct(species), .groups = 'drop') |> 
    ungroup() |> 
    filter(species_combination %in% list.forest) |>
    group_by(species_combination) |> 
    summarise(sgdd=mean(sgdd),
              wai=mean(wai),
              n=n()) |> 
    mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
           waib = 1/(1 + wai), PC1=0, PC2=0, N = 2, SDM = 0,
           file.ipm=paste0("rds/",species_combination,"_meanClim.rds"))
  
  return(mean_forest_list)
}

#' Function to compute IPM for all


#' Function to make a list of simulations till equilibrium
#' @param species.combination.select 
create_simulation_equil_list = function(species.combination.select){
  # create all "partner" forest : forest without the targetted species
  list.forest.bis<-species.combination.select |> 
    dplyr::select(species,species_combination,
                  wai_id,sgdd_id,ID.spclim,clim_lab,
                  wai,sgdd,sgdd2,sgddb,waib,PC1,PC2,N,SDM) |> 
    rowwise() |> 
    mutate(s_p=gsub(" ","_",species),
           species_partner=if_else(sub(paste0(s_p,"\\."),"",species_combination)==species_combination,
                                   if_else(sub(paste0("\\.",s_p),"",species_combination)==species_combination,
                                           NA,
                                           sub(paste0("\\.",s_p),"",species_combination)),
                                   sub(paste0(s_p,"\\."),"",species_combination)),
           species_combination=case_when(species_partner==s_p~NA,
                                         TRUE~species_partner)) |> 
    ungroup() |> 
    dplyr::select(-species_partner,-s_p) |> 
    filter(!is.na(species_combination))
  list.forest<-species.combination.select |> 
    dplyr::select(species,species_combination,
           wai_id,sgdd_id,ID.spclim,clim_lab,
           wai,sgdd,sgdd2,sgddb,waib,PC1,PC2,N,SDM) |> 
    rbind(list.forest.bis) |> 
    unique() |> 
    rowwise() |> 
    mutate(forest.real=grepl(gsub(" ","_",species),species_combination)) |> 
    ungroup() |> 
    arrange(species,ID.spclim) |> 
    mutate(simul_eq=row_number())
  
  id.simul_eq = list.forest$simul_eq
  id.simul_forest = list.forest[list.forest$forest.real,"simul_eq"][[1]]
  return(list(list.forests=list.forest,
              id.simul_eq=id.simul_eq,
              id.simul_forest=id.simul_forest))
  
}

#' Function to make a list of simulations for perturbations
#' @param sim_forest_list
#' @param species.list.disturbance
create_simulation_dist_list = function(sim_forest_list,
                                       species.list.disturbance){
  list.forest=sim_forest_list |> 
    mutate(list.sp=strsplit(species_combination,split="[.]")) |> 
    rowwise() |> 
    mutate(is.dist=(sum(sapply(list.sp,function(x)grepl(x,species.list.disturbance)))==length(list.sp)))
  id.simul_dist=list.forest[list.forest$forest.real&list.forest$is.dist,"simul_eq"][[1]]
  return(id.simul_dist)
  
}





#' Function to make a list of file s outputs of simulations till equilibrium
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param id_forest id of the forest to simulate
make_simulations_equilibrium = function(species.combination,
                                        species_list, 
                                        species_object,
                                        harv_rules.ref,
                                        id_forest){
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  for(i in 1:length(species.in)){
    id.species.obj=species_list[species_list$ID.spclim==clim &
                                  species_list$species==sp &
                                  species_list$species_combination==species.in[i],
                                "id.species.obj"][[1]]
    # Identify the file in species containing species i
    species.file.i = species_object[id.species.obj]
    # Store the file in the list
    list.species[[i]] = readRDS(species.file.i)
    
  }
  
  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
  
  
  sim.in = sim_deter_forest(forest.in, 
                            tlim = 4000,
                            equil_time = 50000, 
                            equil_dist = 2000, 
                            equil_diff = 0.5, 
                            harvest = "default", 
                            SurfEch = 0.03,
                            verbose = TRUE)
  
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_equilibrium/", species.comb, ".rds")
  # Save simulation in a rdata
  create_dir_if_needed(forest.file)
  saveRDS(sim.in, forest.file)
  
  return(forest.file)
  
}


#' Function to make a list of invasion simulations
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param harv_rules.ref rules for harvesting
#' @param sim_equil path of simulation until equilibrium
#' @param threshold_pop dimater threshold of initial population
#' @param id_forest id of the forest to simulate
make_simulations_invasion = function(species.combination,
                                     species_list,
                                     species_object,
                                     harv_rules.ref,
                                     sim_equil,
                                     threshold_pop=0,
                                     id_forest){
  print(id_forest)
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  if(length(species.in)>1){
    partner.comb=if_else(sub(paste0(s_p,"\\."),"",species.comb)==species.comb,
                         if_else(sub(paste0("\\.",s_p),"",species.comb)==species.comb,
                                 "Problem",
                                 sub(paste0("\\.",s_p),"",species.comb)),
                         sub(paste0(s_p,"\\."),"",species.comb))
    simul_eq.partner=species.combination |> 
      filter(species==sp &
               species_combination==partner.comb &
               ID.spclim == clim) |> 
      pull(simul_eq)
    simul_eq.species=species.combination |> 
      filter(species==sp &
               species_combination==s_p &
               ID.spclim == clim) |> 
      pull(simul_eq)
    # Read the simulation at equilibrium
    sim_equilibrium.partner.in = readRDS(sim_equil[simul_eq.partner])
    sim_equilibrium.species.in = readRDS(sim_equil[simul_eq.species])
    
    sim_equilibrium.in=rbind(sim_equilibrium.partner.in,
                             sim_equilibrium.species.in)
    
  }else{
    sim_equilibrium.in= readRDS(sim_equil[id_forest])
  }
  
  
  reached_equil = ifelse(
    is.na(sum((sim_equilibrium.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  
  # Only make the simulation if population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      id.species.obj=species_list[species_list$ID.spclim==clim &
                                    species_list$species==sp &
                                    species_list$species_combination==species.in[i],
                                  "id.species.obj"][[1]]
      # Identify the file in species containing species i
      species.file.i = species_object[id.species.obj]
      # Store the file in the list
      list.species[[i]] = readRDS(species.file.i)
      
      if(species.in[i]!=s_p){
        equil.i = sim_equilibrium.in %>%
          filter(var == "n", equil, species == species.in[i]) %>% 
          pull(value)
      }else{
        # Extract the equilibrium for species i
        equil.i = sim_equilibrium.in %>%
          filter(var == "n", equil, species == species.in[i]) %>%  
          mutate(value=case_when(size>threshold_pop~0,
                                 TRUE~value)) |>
          mutate(BAtot=sum((size/2000)^2*pi*value)) |> 
          pull(value)
      }
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
      
    }
      # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    sim.in = sim_deter_forest(
      forest.in, tlim = 4000, equil_time = 4000, 
      SurfEch = 0.03, verbose = TRUE)
  } else {
    sim.in = matrix()
  }
  
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_invasion/", species.comb, ".rds")
  # Save simulation in a rdata
  create_dir_if_needed(forest.file)
  saveRDS(sim.in, forest.file)
  
  # Return output list
  return(forest.file)
}


#' Function to make a list of simulations with disturbance
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param harv_rules.ref rules for harvesting
#' @param sim_equil path of simulation until equilibrium
#' @param disturb_coef.in table of disturbance coef per species
#' @param disturbance.df_storm disturbance characteristics
#' @param id_forest id of the forest to simulate
make_simulations_disturbance = function(species.combination,
                                     species_list,
                                     species_object,
                                     harv_rules.ref,
                                     sim_equil,
                                     disturb_coef.in,
                                     disturbance.df_storm,
                                     id_forest){
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equil[id_forest])
  
  # Checked that the population reached equilibrium
  reached_equil = ifelse(
    is.na(sum((sim_equilibrium.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  
  
  # Only make the simulation if population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      
      id.species.obj=species_list[species_list$ID.spclim==clim &
                                    species_list$species==sp &
                                    species_list$species_combination==species.in[i],
                                  "id.species.obj"][[1]]
      # Identify the file in species containing species i
      species.file.i = species_object[id.species.obj]
      # Store the file in the list
      list.species[[i]] = readRDS(species.file.i)
      
      
      # Extract the equilibrium for species i
      equil.i =  sim_equilibrium.in %>%
        filter(var == "n", equil, species == species.in[i]) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
      
      
      # Add disturbance coefficients
      list.species[[i]]$disturb_coef <- filter(disturb_coef.in, 
                                               species == species.in[i])
    }
    
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    sim.in = sim_deter_forest(
      forest.in, tlim = 4000, equil_time = 4000, 
      disturbance = disturbance.df_storm,
      verbose = TRUE)
  } else {
    sim.in = matrix()
  }
  
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_disturbance/", species.comb, ".rds")
  # Save simulation in a rdata
  create_dir_if_needed(forest.file)
  saveRDS(sim.in, forest.file)
  
  # Return output list
  return(forest.file)
}


#' Function to make a list of simulations with ibm
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param harv_rules.ref rules for harvesting
#' @param sim_equil path of simulation until equilibrium
#' @param threshold_pop dimater threshold of initial population
#' @param id_forest id of the forest to simulate
make_simulation_ibm <- function(species.combination,
                                species_list,
                                species_object,
                                harv_rules.ref,
                                sim_equil,
                                disturb_coef.in,
                                disturbance.df_storm,
                                id_forest){
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equil[id_forest])
  
  # Checked that the population reached equilibrium
  reached_equil = ifelse(
    is.na(sum((sim_equilibrium.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  
  
  # Only make the simulation if population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      
      id.species.obj=species_list[species_list$ID.spclim==clim &
                                    species_list$species==sp &
                                    species_list$species_combination==species.in[i],
                                  "id.species.obj"][[1]]
      # Identify the file in species containing species i
      species.file.i = species_object[id.species.obj]
      # Store the file in the list
      list.species[[i]] = readRDS(species.file.i)
      
      
      # Extract the equilibrium for species i
      equil.i =  sim_equilibrium.in %>%
        filter(var == "n", equil, species == species.in[i]) %>% 
        pull(value)
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
      
      
      # Add disturbance coefficients
      list.species[[i]]$disturb_coef <- filter(disturb_coef.in, 
                                               species == species.in[i])

      eval(parse(text=paste0("list.species$", species.in[i], "$IPM$fit <- fit_", species.in[i])))
      
      
    }
    
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    sim.in = sim_indiv_forest(
      forest.in, tlim = 500,
      verbose = TRUE)
  } else {
    sim.in = matrix()
  }
  
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_disturbance/", species.comb, ".rds")
  # Save simulation in a rdata
  create_dir_if_needed(forest.file)
  saveRDS(sim.in, forest.file)
  
  # Return output list
  return(forest.file)
  
}
# 
# sim.tot=data.frame(matrix(ncol=8,nrow=0))
# colnames(sim.tot)=c("species","var","time","mesh","size","equil","value","nsim")
# for (i in 1:100){
#   print(i)
#   sim.in = sim_indiv_forest(
#     forest.in, tlim = 1000,
#     verbose = FALSE)
#   sim.tot=rbind(sim.tot,
#                 sim.in |>
#                   filter(species=="Abies_alba") |> 
#                   mutate(nsim=i))
# }
# sim.tot |> 
#   dplyr::filter(var == "N", ! equil) %>%
#   group_by(nsim) |>
#   mutate(is.ext=any(value==2)) |>
#   ungroup() |>
#   filter(is.ext) |>
#   # filter(nsim<6) |> 
#   ggplot(aes(x = time, y = value, color = nsim, group=nsim)) +
#   geom_line(linewidth = .4) + ylab("N") 

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



#' Function to compute invasion rate for a simulation
#' @param species.combination table with all species combi for each climate cat
#' @param sim_invasion table with all species fit info 
#' @param id_forest id of the forest to simulate
#' @param fit.list.allspecies list of species with ipm objects, here to get the delay
get_invasion_rate<-function(species.combination,
                            id.simul_forest,
                            sim_invasion,
                            fit.list.allspecies){
  out=species.combination |> 
    filter(simul_eq %in% id.simul_forest) |> 
    mutate(inv_mean=NA_real_,
           inv_max=NA_real_,
           inv_50=NA_real_)
  for (i in 1:length(sim_invasion)){
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_invasion)))
    
    
    # Read simulation i
    sim.i = readRDS(sim_invasion[i])
    
    if(dim(sim.i)[1]>1){ # check simulation was performed (i.e. equil was reached)
      # get simulation index
      id.sim.i=id.simul_forest[i]
      species.i=species.combination[id.sim.i,"species"][[1]]
      s_p.i=gsub(" ","_",species.i)
      fit.species.i=fit.list.allspecies[[s_p.i]]
      delay.i=as.numeric(fit.species.i$info[["delay"]])
      
      derivative.i <- sim.i |> 
        filter(species==s_p.i,var=="BAsp",time>delay.i) |> 
        mutate(der=(value-lag(value))/(time-lag(time)),
               der2=(der-lag(der))/(time-lag(time)))
      
      first.extr<-derivative.i |>  
        filter(time<1500) |> 
        mutate(sign_eq=(sign(der2)==sign(lag(der2)))) |> 
        filter(sign_eq==FALSE) |> slice(1) |> 
        dplyr::select(time) |> pull(time)
      
      inv.1<-derivative.i |> 
        filter(time<(first.extr-5)) |> 
        summarise(inv_mean=mean(der,na.rm=TRUE),
                  inv_max=max(der,na.rm=TRUE))
      
      inv.2<-derivative.i |> 
        filter(time<delay.i+50) |> 
        summarise(inv_50=mean(der,na.rm=TRUE))
      out[out$simul_eq==id.sim.i,c("inv_mean","inv_max","inv_50")]=
        list(inv.1$inv_mean[[1]],
             inv.1$inv_max[[1]],
             inv.2$inv_50[[1]])
    }

  }
  
  return(out)
}

#' Function to compute invasion rate for a simulation
#' @param species.combination table with all species combi for each climate cat
#' @param sim_disturbance table with all species fit info 
#' @param id_forest id of the forest to simulate
#' @param fit.list.allspecies list of species with ipm objects, here to get the delay
get_resilience_metrics<-function(species.combination,
                                 sim_dist.id,
                                 sim_disturbance,
                                 fit.list.allspecies,
                                 disturbance.df_storm){
  
  out=species.combination |> 
    filter(simul_eq %in% sim_dist.id) |> 
    mutate(resistance = NA_real_, recovery = NA_real_, resilience = NA_real_, 
           t0 = NA_real_, thalf = NA_real_, SD = NA_real_, BA_diff = NA_real_, 
           BA_eq = NA_real_, dbh_mean = NA_real_, dbh_q10 = NA_real_, 
           dbh_q90 = NA_real_, dbh_mean_postdist = NA_real_, 
           dbh_q10_postdist = NA_real_, dbh_q90_postdist = NA_real_)
  
  # Identify disturbance time
  tdist = min(disturbance.df_storm$t)
  
  for (i in 1:length(sim_disturbance)){
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_disturbance)))
    
    # get simulation index
    id.sim.i=sim_dist.id[i]
    species.i=species.combination[id.sim.i,"species"][[1]]
    s_p.i=gsub(" ","_",species.i)
    fit.species.i=fit.list.allspecies[[s_p.i]]
    # delay.i=as.numeric(fit.species.i$info[["delay"]])
    
    # Read simulation i
    sim.i = readRDS(sim_disturbance[i]) |> 
      filter(species==s_p.i)
    
    if(!is.na(sim.i[1, 1])){
      
      # mean dbh at equilibrium and after disturbance
      dbh_i = sim.i %>%
        filter(var == "n") %>%
        filter(time %in% c(1, (max(disturbance.df_storm$t)+1))) %>%
        group_by(size, time) %>%
        summarize(ntot = sum(value)) %>%
        ungroup() %>% group_by(time) %>%
        filter(size > 0) %>%
        mutate(ntot_size = ntot*size) %>%
        summarize(mean_dbh = weighted.mean(size, w = ntot), 
                  q10_dbh = weighted.quantile(size, w = ntot, prob = 0.1), 
                  q90_dbh = weighted.quantile(size, w = ntot, prob = 0.9))
      out$dbh_mean[i] <- subset(dbh_i, time == 1)$mean_dbh
      out$dbh_q10[i] <- subset(dbh_i, time == 1)$q10_dbh
      out$dbh_q90[i] <- subset(dbh_i, time == 1)$q90_dbh
      out$dbh_mean_postdist[i] <- subset(dbh_i, time != 1)$mean_dbh
      out$dbh_q10_postdist[i] <- subset(dbh_i, time != 1)$q10_dbh
      out$dbh_q90_postdist[i] <- subset(dbh_i, time != 1)$q90_dbh
      
      # Format the output
      data.i <- sim.i %>%
        filter(var == "BAsp") %>%
        filter(!equil) %>%
        group_by(time) %>%
        summarize(BA = sum(value))
      
      ## Calculate stability before disturbance (to check equilibrium)
      out$SD[i] = sd(subset(data.i, time %in% c(1:(tdist-1)))$BA)
      out$BA_diff[i] = diff(range(subset(data.i, time %in% c(1:(tdist-1)))$BA))
      
      ## Calculate resistance
      #  - Basal area at equilibrium
      Beq.i = mean((data.i %>% filter(time < min(disturbance.df_storm$t)))$BA)
      out$BA_eq[i] = Beq.i
      # - Basal area after disturbance
      Bdist.i = (data.i %>% filter(time == max(disturbance.df_storm$t)+1))$BA
      # - Resistance : logit of the percentage of basal area that survived 
      #out$resistance[i] = Beq.i/(Beq.i - Bdist.i)
      out$resistance[i] = log((Bdist.i/Beq.i)/(1 - (Bdist.i/Beq.i)))
      
      ## Calculate recovery
      #  - Time at which population recovered fully
      Rec.time.i = min((data.i %>% 
                          filter(time > max(disturbance.df_storm$t)) %>%
                          filter(BA > Beq.i))$time)
      # - Basal area 20 years after disturbance
      Bdist20.i = (data.i %>% filter(time == max(disturbance.df_storm$t)+21))$BA
      # - Recovery = slope of BA increase in teh 20 years after disturbance
      out$recovery[i] = abs(Bdist20.i - Bdist.i)/20
      
      ## Calculate resilience
      out$resilience[i] <- 1/sum((data.i %>%
                                    mutate(BA0 = .[which(.$time == 1), "BA"]) %>%
                                    mutate(diff = abs(BA - BA0)))$diff)
      
      ## Calculate t0
      #  - Time at which population recovered to 5% of the basal area lost
      Rec.0.time.i = min((data.i %>% 
                            filter(time > max(disturbance.df_storm$t)) %>%
                            filter(BA > (Beq.i + 19*Bdist.i)/20))$time)
      # - Recovery = time to recover minus time of disturbance
      out$t0[i] = Rec.0.time.i - max(disturbance.df_storm$t)
      
      ## Calculate thalf
      #  - Time at which population recovered to 50% of the basal area lost
      Rec.half.time.i = min((data.i %>% 
                               filter(time > max(disturbance.df_storm$t)) %>%
                               filter(BA > (Beq.i + Bdist.i)/2))$time)
      # - Recovery = time to recover minus time of disturbance
      out$thalf[i] = Rec.half.time.i - max(disturbance.df_storm$t)
      
      
    }
    
  }
  return(out)
}

#' Function to extract parameter of vital rates
#' @param species.combination table with all species combi for each climate cat
#' @param sim_disturbance table with all species fit info 
#' @param id_forest id of the forest to simulate
#' @param fit.list.allspecies list of species with ipm objects, here to get the delay
get_vitalrates_pars<-function(species_list.select,
                              species_object,
                              fit.list.allspecies){
  pars_list<-  species_list.select |> 
    cbind(file_real=unname(species_object))|> 
    filter(gsub(" ","_",species)==species_combination) 
  
  for (sp in unique(pars_list$species)){
    print(sp)
    s_p=sub(" ","_",sp)
    species.fit=fit.list.allspecies[[s_p]]
    
    for (i in 1:dim(pars_list)[1]){
      print(i)
      species.obj=readRDS(pars_list$file_real[i])
      species.obj$IPM$fit
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}



































