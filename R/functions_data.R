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


#%%%%%%%%%%%%%%%%%%%%%%%%
#### GET FUNDIV DATA ####
#%%%%%%%%%%%%%%%%%%%%%%%%

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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### GET DISTURBANCE DATA ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Function to get species meta data
#' @param species.list.ipm species to simulate
get_species_meta<-function(species.list.ipm){
  species.meta<- data.frame(species=species.list.ipm) |> 
    separate(species,into=c("genus","sp"),remove = FALSE)|> 
    mutate(taxa=case_when(genus%in%c("Pinus","Abies","Picea","Larix")~"conifer",
                          TRUE~"broadleaf")) 
  return(species.meta)
  
}

#' Function to assoacite disturbance parameters to species
#' @param disturb_coef.raw file with param from the model
#' @param species lsit that are simulable from ipm
get_disturb_coef<-function(disturb_coef.raw,
                           species_meta,
                           species.list.ipm){
  disturb_coef.raw<-disturb_coef.raw |> 
    dplyr::select(-X) |> 
    mutate(species=gsub(" ","_",species))|> 
    filter(disturbance=="storm") |> 
    relocate(species,.before=disturbance)
  other.conifer<-subset(disturb_coef.raw,species=="Other_conifer",select=-species)
  other.broadleaf<-subset(disturb_coef.raw,species=="Other_broadleaf",select=-species)
  disturb_coef_out<-data.frame(species=species.list.ipm) |> 
    left_join(disturb_coef.raw) |> 
    left_join(species_meta) |> 
    mutate(real.coef=case_when(is.na(disturbance)~FALSE,
                               TRUE~TRUE))
  
  for(i in 1:dim(disturb_coef_out)[1]){
    if(!disturb_coef_out$real.coef[i]){
      if(disturb_coef_out$taxa[i]=="conifer"){
        disturb_coef_out[i,2:10]<-other.conifer
      }else{
        disturb_coef_out[i,2:10]<-other.broadleaf
      }
    }
  }
  return(disturb_coef_out)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### MAKE SPECIES COMBI ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to generate a list with climate per species
#' @description create climate categories so that numbers of plots in each categories
#' falls between max and min_cat
#' @param FUNDIV_data whole dataset
#' @param species.list.ipm list of species
#' @param max_cat maximum number of categories per species
#' @param min_cat minimum
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

  step_cat<- data.frame(name=c("sgdd","wai"),
                        step=c(95,0.095))
  
  max_obs<-max_cat+1
  while(max_obs>max_cat){
    print(max_obs)
    step_cat$step[1]=step_cat$step[1]+5
    step_cat$step[2]=step_cat$step[2]+0.005
    
    
    # define numbers of breaks for sggd/wai for each species according to step_cat
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
    
    
    # associate each plot x species to a cliamate category of the species
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
  
  return(list(FUNDIV_plotcat=FUNDIV_plotcat,
              species.cat=condi.init))
}


#' Function to generate a list with climate per species over first PCA axis
#' @description create a constant number of climate categories by species, over PCA1
#' @param FUNDIV_data whole dataset
#' @param species.list.ipm list of species
#' @param n_cat number of categories
make_climate_cat_pca <- function(FUNDIV_data,
                                 species.list.ipm,
                                 n_cat=10){
  # Step 1: Filter and group the data by plot
  # Keep data for species in the given list and summarize each plot by selecting relevant columns
  FUNDIV_plot = FUNDIV_data |> 
    filter(species %in% gsub("_", " ", species.list.ipm)) |>  # Filter for species in the list (replace underscores with spaces)
    dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2, species, BAtot) |>  # Select relevant columns
    group_by(plotcode, longitude, latitude, sgdd, wai, pca1, pca2, species) |>  # Group by these variables
    summarize(BA = max(BAtot)) |>  # For each plot, take the maximum basal area (BA)
    distinct() |>  # Remove duplicates
    ungroup()
  
  
  # Step 2: Create climate categories for each species based on PCA1
  FUNDIV_plotcat <- FUNDIV_plot |> 
    filter(species %in% gsub("_", " ", species.list.ipm)) |>  # Filter again for species
    group_by(species) |>  # Group by species
    mutate(clim_cat = cut(pca1,  # Create climate categories based on PCA1
                          breaks = quantile(pca1, probs = seq(0, 1, length.out = n_cat + 1)),  # Use quantiles for equal-sized categories
                          include.lowest = TRUE)) |>  
    # Separate the climate category into lower and upper bounds
    tidyr::separate_wider_delim(cols = "clim_cat",
                                names = c("clim_low", "clim_up"),
                                delim = ",",
                                cols_remove = FALSE) |> 
    # Convert the category bounds from string to numeric by trimming parentheses
    mutate(across(matches(c("low")), ~as.numeric(stringr::str_sub(., 2, -1))),
           across(matches(c("up")), ~as.numeric(stringr::str_sub(., 1, -2)))) |>
    ungroup()  |> 
    group_by(species) %>%
    # Assign an ID to each climate category for each species
    mutate(clim_id = as.integer(factor(clim_cat))) %>%
    ungroup() 
  
  # Step 3: Summarize the climate categories for each species
  species.cat <- FUNDIV_plotcat |> 
    group_by(species, clim_id, clim_low, clim_up) |>  # Group by species and climate categories
    summarize(n_plot = n(),  # Count the number of plots in each climate category
              wai = mean(wai, na.rm = TRUE),  # Calculate the mean WAI (water availability index)
              sgdd = mean(sgdd, na.rm = TRUE)) |>  # Calculate the mean SGDD (growing degree days)
    # Create additional variables used in further models
    mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1 / sgdd, 
           waib = 1 / (1 + wai), PC1 = 0, PC2 = 0, N = 2, SDM = 0) |> 
    ungroup() |> 
    group_by(species) |> 
    # Create unique IDs for each species and climate category combination
    mutate(ID.spclim = row_number(),
           clim_lab = paste0("pca_", clim_id),
           ID.species = cur_group_id()) |> 
    ungroup()
  
  return(list(FUNDIV_plotcat=FUNDIV_plotcat,
              species.cat=species.cat))
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
                                      species.list.disturbance,
                                      species.list.ipm, 
                                      nsp_per_richness=10,
                                      prop_threshold=0.8){
  s_p=species.list.disturbance[sp_id]
  sp=gsub("_"," ",s_p)
  print(sp)
  species.combinations <- FUNDIV_data |>
      # compute basal area ratio of species per plot
      filter(dbh1>0) |>
      group_by(plotcode) |> 
      mutate(BAsum=sum(ba_ha1)) |> 
      group_by(plotcode,species) |> 
      mutate(BAsumsp=sum(ba_ha1), 
             sp_ratio=BAsumsp/BAsum) |> 
      ungroup() |> 
      # join with climate cat of the targetted species
      left_join(FUNDIV_plotcat |>
                  filter(species==sp) |>
                  dplyr::select(plotcode,clim_id,clim_low,clim_up), #wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up
                by="plotcode") |> 
      dplyr::select(treecode,plotcode,species,sp_ratio,clim_id,clim_low,clim_up) |> #wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up
      # filter out competitors that are not present enough
      filter(!(species!=sp&sp_ratio<0.1)) |> 
      # Group by the climatic condition and plot
      group_by(clim_id, plotcode) %>% #wai_id, sgdd_id
      # filter plots with only target species
      filter(sp%in%species)|>
      # Summarize the species combination in each group, and species richness
      mutate(species_combination = paste(sort(unique(gsub(" ","_",species))), collapse = "."),
             n_species = n_distinct(species))  |> 
      ungroup() |>  
      # create one row per plot
      dplyr::select(plotcode,clim_id, species_combination,n_species) |> unique() |> 
      # count each unique combination's frequency within clim cat
      count(clim_id, species_combination,n_species) |>
      # pre filter
      filter(n>2) |> 
      filter(n_species<nsp_per_richness) |>
      # mark wether species combination are simulable
      rowwise() |> 
      mutate(is.sim=(n_species==sum(sapply(species.list.ipm,function(x)grepl(x,species_combination)))),
             is.sp=!(species_combination==s_p)) |> 
      # Compute proportion and cumulated prop
      group_by(clim_id) |>
      arrange(clim_id, desc(n)) |> 
      mutate(prop=(n*is.sp)/sum(n*is.sp), # prop not accounting for species alone
             prop_all=n/sum(n),
             prop_cum=cumsum(prop),
             prop_cum_all=cumsum(prop_all),
             prop_cum_filt=is.sp*cumsum(prop*is.sim),
             max_cum_filt=max(prop_cum_filt))   
    # species.combinations |> 
    #   group_by(species_combination) |> 
    #   summarise(ntot=sum(n)) |> 
    #   ungroup() |>
    #   mutate(prop=ntot/sum(ntot)) |> 
    #   filter(species_combination==s_p) |> pull(prop) |> 
    #   print()
  out<-condi.init |>
    filter(species==sp) |> 
    left_join(species.combinations |> 
                filter((species_combination==s_p|prop_cum_filt<0.8)&is.sim),
              by="clim_id") |> 
    arrange(clim_id,n_species)

  return(out)
}


#' Function that computes climate boundaries for matrix mu fits
#' @param species.combination list of species combinations
make_clim_boundaries <- function(species.combination){
  out<-setNames(data.frame(matrix(ncol=8,nrow = 0)),
                nm =c("species","N","sgdd","sgdd2","sgddb","wai","wai2","waib"))
  for (species in unique(species.combination$species_combination)){
    s_p=gsub(" ","_",species)
    sp.combi<-species.combination %>% 
      filter(grepl(s_p,species_combination)) 
    sgdd_max=max(sp.combi$sgdd)
    wai_max=max(sp.combi$wai)
    sgdd_min=min(sp.combi$sgdd)
    wai_min=min(sp.combi$wai)
    sgdd_med=median(sp.combi$sgdd)
    wai_med=median(sp.combi$wai)
    new_rows <- data.frame(
      species = c(s_p, s_p, s_p),
      N = c(3, 2, 1),
      sgdd = c(sgdd_min, sgdd_med, sgdd_max),
      sgdd2 = c(sgdd_min^2, sgdd_med^2, sgdd_max^2),
      sgddb = c(1/sgdd_min, 1/sgdd_med, 1/sgdd_max),
      wai = c(wai_max, wai_med, wai_min),
      wai2 = c(wai_max^2, wai_med^2, wai_min^2),
      waib = c(1/wai_max, 1/wai_med, 1/wai_min) # errueur ici!!!
    )
    out<-rbind(out,
               new_rows)
  }
  return(out)
}


#' Look for species main competitors
#' @param FUNDIV_data fundiv data
#' @param species.list.ipm list of species in the ipm
#' @param sp_id species id
get_species_competitors<-function(FUNDIV_data,
                                  species.list.ipm,
                                  sp_id){
  s_p=species.list.ipm[sp_id]
  sp=gsub("_"," ",s_p)
  print(sp)
  competitor_sp <-  FUNDIV_data |>
    # compute basal area ratio of species per plot
    filter(dbh1>0) |>
    group_by(plotcode) |> 
    mutate(BAsum=sum(ba_ha1)) |> 
    group_by(plotcode,species) |> 
    mutate(BAsumsp=sum(ba_ha1), 
           sp_ratio=BAsumsp/BAsum) |> 
    ungroup() |>
    dplyr::select(plotcode,species,sp_ratio,BAsumsp) |>
    unique()  |> 
    # filter out competitors that are not present enough
    filter(!(species!=sp&sp_ratio<0.1)) |> 
    # Group by the climatic condition and plot
    group_by(plotcode) %>% 
    # filter plots with only target species
    filter(sp%in%species)|> 
    group_by(species) |> 
    summarise(sp_ratio_mean=mean(sp_ratio),
              BA_tot=sum(BAsumsp)) |> 
    arrange(desc(BA_tot)) |> 
    rowwise() |> 
    mutate(is.sim=gsub(" ","_",species)%in%species.list.ipm,
           species_target=sp) |>
    filter(species!=sp) |> 
    ungroup() |> 
    slice(1:20)
  return(competitor_sp)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### INITIALIZE IPM SIMUL ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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


#' Function create a list of all species object to be created
#' @param species.combination table of species combination for each climate cat
make_species_list <- function(species.combination){
  species_list<-species.combination |> 
    dplyr::select(species,clim_id,ID.spclim,clim_lab, # climate info #wai_id,sgdd_id
           wai,sgdd,sgdd2,wai2,sgddb,waib,PC1,PC2,N,SDM, # data for IPM clim
           species_combination) |> # eponyme
    mutate(species_combination=strsplit(species_combination,"\\.")) |> 
    unnest(cols = species_combination) |> 
    unique() |> 
    mutate(file.ipm=paste0("rds/",gsub(" ","_",species),"/clim_",ID.spclim,"/",species_combination,".rds"))
  
  return(species_list)
}

#' Function create a list of all species object to be created
#' @param species.combination table of species combination for each climate cat
#' @param pars_list
make_species_select_list_elast <- function(species.combination,
                                           pars_list){
  split_pars <- lapply(pars_list, function(x) str_split(x, pattern = "-")[[1]])
  pars_df <- data.frame(
    s_p = sapply(split_pars, `[`, 1),
    gr = sapply(split_pars, `[`, 2),
    param = sapply(split_pars, `[`, 3),
    elast=pars_list
  ) |> 
    mutate(sp=str_replace(s_p,"_"," "))
  
  species_list<-species.combination |> 
    left_join(pars_df,by=c("species"="sp")) |>
    dplyr::select(species,s_p,elast,clim_id,ID.spclim,clim_lab, # climate info #wai_id,sgdd_id
                  wai,sgdd,sgdd2,wai2,sgddb,waib,PC1,PC2,N,SDM, # data for IPM clim
                  species_combination) |> # eponyme
    mutate(species_combination=strsplit(species_combination,"\\.")) |> 
    unnest(cols = species_combination) |> 
    mutate(species_combination=factor(species_combination),
           id.species.mu.obj=as.numeric(species_combination)) |> 
    mutate(species_combination=case_when(s_p==species_combination~elast,
                                         TRUE~species_combination)) |> 
    mutate(elast=factor(elast),
           id.species.elast.obj=as.numeric(elast)) |> 
    unique() 
  
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
                            species.select, 
                            clim_bound,
                            sp_id){
  sp=species.select[sp_id]
  s_p=gsub(" ","_",sp)
  print(sp)

  # Make IPM
  IPM.mu <- make_mu_gr(
    species = s_p, fit = fit.list.allspecies[[s_p]],
    climate = subset(clim_bound,species == s_p, select = -species),
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit.list.allspecies[[s_p]]) * 1.1),
    verbose = TRUE, stepMu = 0.0001)
  
    # Create species object 
  species.in = species(IPM=IPM.mu,
                       init_pop = def_initBA(20),
                       harvest_fun = def_harv, 
                       disturb_fun = def_disturb)
  # mu.file=paste0("rds/",s_p,"/",s_p,"_mu.rds")
  species.file=paste0("rds/",s_p,"/",s_p,"_species.rds")
  
  
  # Save species object in a rdata
  # create_dir_if_needed(mu.file)
  # saveRDS(IPM.mu, mu.file)
  create_dir_if_needed(species.file)
  saveRDS(species.in, species.file)
  
  # Return output list
  return(species.file)
  
}


#' Function to make a species mu matrix with a slightly modification of one pars,
#'  create species object, and save it as rds and return filename
#' @param fit.list.allspecies demographic parameter for all species
#' @param species.select list a species to run
#' @param delta variation to impose on parameter, in percentage
#' @param pars_list list of parameters
#' @param sp_id species id to run
#' @param pars_id parameter id to modify
#' @author Anne Baranger
make_species_mu_elast <- function(fit.list.allspecies,
                                  delta=0.01,
                                  pars_select){
  # get species
  s_p=strsplit(pars_select, split = '-')[[1]][1]
  sp=gsub("_"," ",s_p)
  print(sp)
  # get vital rates and parameter to disturb
  vr=strsplit(pars_select, split = '-')[[1]][2]
  pars=strsplit(pars_select, split = '-')[[1]][3]
  # get fit
  fit_sp=fit.list.allspecies[[s_p]]
  # disturb parameter
  fit_sp[[vr]]$params_m[[pars]]<-fit_sp[[vr]]$params_m[[pars]]*(1+delta)
  
  
  # Make IPM
  IPM.mu <- make_mu_gr(
    species = s_p, fit = fit_sp,
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit.list.allspecies[[s_p]]) * 1.1),
    verbose = TRUE, stepMu = 0.0001)
  
  # Create species object 
  species.in = species(IPM=IPM.mu,
                       init_pop = def_initBA(20),
                       harvest_fun = def_harv, 
                       disturb_fun = def_disturb)
  species.file=paste0("rds/",s_p,"/",
                      gsub(":","x", pars_select), #because ":" are not supported in filenames
                      ".rds")
  
  # Save species object in a rdata
  create_dir_if_needed(species.file)
  saveRDS(species.in, species.file)
  
  # Return output list
  return(species.file)
  
}


#' Function to make a list of forest and mean associated climate for mean IPMs
#' @param FUNDIV_data 
#' @param sim_forest_list 
#' @author Anne Baranger
make_mean_forest_list <- function(FUNDIV_data,
                                  sim_forest_list=sim_forest_list[["list.forests"]]){
  # extract all combination of species
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

#' Function to compute IPM for all mean species
#' @param sim_mean_forest_list
make_mean_species_list <- function(sim_mean_forest_list){
  species_list<-sim_mean_forest_list |> 
    mutate(species_combination=strsplit(species_combination,"\\.")) |> 
    unnest(cols = species_combination) |> 
    unique() |> 
    mutate(file.ipm=paste0("rds/",gsub(" ","_",species),"/clim_",ID.spclim,"/",species_combination,".rds"))
}


#' Function to make a list of simulations till equilibrium
#' @param species.combination.select 
create_simulation_equil_list = function(species.combination.select){
  # create all "partner" forest : forest without the targetted species
  list.forest.bis<-species.combination.select |> 
    dplyr::select(species,species_combination,
                  clim_id,ID.spclim,clim_lab, #wai_id,sgdd_id
                  wai,sgdd,sgdd2,wai2,sgddb,waib,PC1,PC2,N,SDM) |> 
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
           clim_id,ID.spclim,clim_lab, #wai_id,sgdd_id
           wai,sgdd,wai2,sgdd2,sgddb,waib,PC1,PC2,N,SDM) |> 
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

#' Function to make a list of simulations till equilibrium for elasticity analysis
#' @param species.combination.select 
#' @param pars_select list of species and params for elasticity
create_simulation_equil_list_elast <- function(species.combination.select,
                                               sim_forest_list,
                                               pars_select){
  split_pars <- lapply(pars_select, function(x) str_split(x, pattern = "-")[[1]])
  pars_df <- data.frame(
    sp = sapply(split_pars, `[`, 1),
    gr = sapply(split_pars, `[`, 2),
    param = sapply(split_pars, `[`, 3),
    elast=pars_select) |>
    mutate(sp=str_replace(sp,"_"," "))
    
  
  # create all "partner" forest : forest without the targetted species
  list.forest.bis<-species.combination.select |> 
    left_join(pars_df,by=c("species"="sp")) |> 
    dplyr::select(species,elast,species_combination,
                  clim_id,ID.spclim,clim_lab, #wai_id,sgdd_id
                  wai,sgdd,sgdd2,wai2,sgddb,waib,PC1,PC2,N,SDM) |> 
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
    filter(!is.na(species_combination)) |> 
    left_join(sim_forest_list[,c("species","species_combination","clim_id","simul_eq")],
              by=c("species","species_combination","clim_id"))
  list.forest<-species.combination.select |> 
    left_join(pars_df,by=c("species"="sp")) |> 
    dplyr::select(species,elast,species_combination,
                  clim_id,ID.spclim,clim_lab, #wai_id,sgdd_id
                  wai,sgdd,wai2,sgdd2,sgddb,waib,PC1,PC2,N,SDM) |> 
    mutate(simul_eq=NA) |> # because these forest are to be simulated again
    rbind(list.forest.bis) |> 
    rowwise() |> 
    mutate(species_combination=gsub(gsub(" ","_",species), 
                                    elast,
                                    species_combination)) |> 
    unique() |> 
    mutate(forest.real=grepl(gsub(" ","_",species),species_combination)) |> 
    ungroup() |> 
    arrange(species,ID.spclim) |> 
    mutate(simul_eq_elast=row_number())
  
  id.simul_eq = list.forest$simul_eq_elast
  id.simul_forest = list.forest[list.forest$forest.real,"simul_eq_elast"][[1]]
  return(list(list.forests=list.forest,
              id.simul_eq=id.simul_eq,
              id.simul_forest=id.simul_forest))
  
}


#' Function to make a list of simulations for perturbations
#' @description
#' function that checks whether all species present in combinations have disturbance
#' parameters, if not mark the combination as such 
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


#' Function to make a list of simulations for perturbations
#' @description
#' function that checks whether all species present in combinations have disturbance
#' parameters, if not mark the combination as such 
#' @param sim_forest_list
#' @param species.list.disturbance
create_simulation_dist_list_elast = function(sim_forest_list,
                                             species.list.disturbance){
  list.forest=sim_forest_list |> 
    mutate(s_p=gsub(" ","_",species)) %>% 
    rowwise() %>% 
    mutate(species_combination=gsub(elast,
                                    s_p,
                                    species_combination),
           list.sp=strsplit(species_combination,split="[.]")) |> 
    mutate(is.dist=(sum(sapply(list.sp,function(x)grepl(x,species.list.disturbance)))==length(list.sp)))
  id.simul_dist=list.forest[list.forest$forest.real&list.forest$is.dist,"simul_eq_elast"][[1]]
  return(id.simul_dist)
  
}


#' Mark simulation where target species was excluded at equilibrium
#' @param list.forest list of simulated forests
#' @param sim_equil files of equil sim
is_species_excluded<-function(sim_forest_list,
                              sim_equil){
  list.forest<-sim_forest_list$list.forests
  list.forest$excluded<-NA
  list.forest$competexcluded<-NA
  for (i in 1:dim(list.forest)[1]){
    if(list.forest[[i,"forest.real"]]){
      sp=list.forest[[i,"species"]]
      s_p=gsub(" ","_",sp)
      competitors=unlist(strsplit(list.forest[[i,"species_combination"]],"\\."))
      competitors=competitors[competitors!=s_p]

      equil.i=readRDS(sim_equil[[list.forest[[i,"simul_eq"]]]])
      if(equil.i$reached_equil){
        BA_equil<-equil.i$distrib_equil %>% 
          filter(size>0) %>% 
          group_by(species) %>% 
          summarize(BAtot=sum((size/2000)^2*pi*value))
        
        # BA equil of targetted species
        if(BA_equil[BA_equil$species==s_p,"BAtot"]<1){
          list.forest$excluded[i]="excluded"
        }else{
          list.forest$excluded[i]="notExcluded"
        }
        # BA equil of compet species
        compet_excl=""
        for(sp_comp in competitors){
          if(BA_equil[BA_equil$species==sp_comp,"BAtot"]<1){
            compet_excl=paste(compet_excl,sp_comp,sep=".")
          }
        }
        list.forest$competexcluded[i]=compet_excl
      }else{
        list.forest$excluded[i]="noEquil"
      }
    }else{
      list.forest$excluded[i]="notForest"
    }
  }
  sim_forest_list$list.forests<-list.forest
  return(sim_forest_list)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%
#### MAKE SIMULATIONS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%

#' Function to make a list of file s outputs of simulations till equilibrium
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param id_forest id of the forest to simulate
make_simulations_equilibrium = function(species.combination, # table of all real/false forests
                                        species_list, # list of species object id by climate
                                        species_object, # path of species_object
                                        harv_rules.ref,
                                        sim.type="mu",
                                        id_forest){
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  if(sim.type=="mu"){id.obj="id.species.mu.obj"}else{id.obj="id.species.obj"}
  for(i in 1:length(species.in)){
    id.species.obj=species_list[species_list$ID.spclim==clim &
                                  species_list$species==sp &
                                  species_list$species_combination==species.in[i],
                                id.obj][[1]]
    # Identify the file in species containing species i
    species.file.i = species_object[id.species.obj]
    # Store the file in the list
    list.species[[i]] = readRDS(species.file.i)
    
  }
  
  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
  
  if(sim.type=="mu"){
    sim.in = sim_deter_forest(forest.in, 
                              tlim = 4000,
                              climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                      "waib", "wai2", "sgdd2", 
                                                                      "PC1", "PC2", "N", "SDM")],
                              equil_time = 50000, 
                              equil_dist = 2000, 
                              equil_diff = 0.5, 
                              harvest = "default", 
                              SurfEch = 0.03,
                              verbose = TRUE)
  }else{
    sim.in = sim_deter_forest(forest.in, 
                              tlim = 4000,
                              equil_time = 50000, 
                              equil_dist = 2000, 
                              equil_diff = 0.5, 
                              harvest = "default", 
                              SurfEch = 0.03,
                              verbose = TRUE)
  }
  
  reached_equil = ifelse(
    is.na(sum((sim.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  distrib_equil = sim.in %>%
    filter(var == "n", equil)
  
  equil.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_equilibrium/", species.comb, ".rds")
  
  # Save simulation in a rdata
  create_dir_if_needed(equil.file)
  saveRDS(list(reached_equil=reached_equil,
               distrib_equil=distrib_equil), 
          equil.file)

  return(equil.file)
  
}


#' Function to get dsitribution of trees at equilibrium, including elesticity
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param id_forest id of the forest to simulate
make_simulations_equilibrium_elast = function(species.combination,
                                              species_list_select_elast, 
                                              species_object_mu,
                                              species_object_mu_elast,
                                              harv_rules.ref,
                                              sim.type="mu",
                                              id_forest){
  sp=species.combination[id_forest,"elast"][[1]]
  s_p=strsplit(sp, split = '-')[[1]][1]
  
  # s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  if(sim.type=="mu"){id.obj="id.species.mu.obj"}else{id.obj="id.species.obj"}
  for(i in 1:length(species.in)){
    if(species.in[i]==sp){
      id.species.obj=species_list_select_elast[species_list_select_elast$ID.spclim==clim &
                                                 species_list_select_elast$elast==sp &
                                                 species_list_select_elast$species_combination==species.in[i],
                                        "id.species.elast.obj"][[1]]
      species.file.i=species_object_mu_elast[id.species.obj]
    }else{
      id.species.obj=species_list_select_elast[species_list_select_elast$ID.spclim==clim &
                                                 species_list_select_elast$elast==sp &
                                                 species_list_select_elast$species_combination==species.in[i],
                                        "id.species.mu.obj"][[1]]
      species.file.i=species_object_mu[id.species.obj]
    }

    # Store the file in the list
    
    list.species[[i]] = readRDS(species.file.i)

    
  }
  
  # Make forest
  forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
  
  if(sim.type=="mu"){
    sim.in = sim_deter_forest(forest.in, 
                              tlim = 4000,
                              climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                      "waib", "wai2", "sgdd2", 
                                                                      "PC1", "PC2", "N", "SDM")],
                              equil_time = 50000, 
                              equil_dist = 2000, 
                              equil_diff = 0.5, 
                              harvest = "default", 
                              SurfEch = 0.03,
                              verbose = TRUE)
  }else{
    sim.in = sim_deter_forest(forest.in, 
                              tlim = 4000,
                              equil_time = 50000, 
                              equil_dist = 2000, 
                              equil_diff = 0.5, 
                              harvest = "default", 
                              SurfEch = 0.03,
                              verbose = TRUE)
  }
  
  reached_equil = ifelse(
    is.na(sum((sim.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  distrib_equil = sim.in %>%
    filter(var == "n", equil)
  
  equil.file=paste0("rds/", s_p, "/clim_", clim,
                    "/sim_equilibrium/", 
                    gsub(":","x", species.comb),
                    ".rds")
  
  # Save simulation in a rdata
  create_dir_if_needed(equil.file)
  saveRDS(list(reached_equil=reached_equil,
               distrib_equil=distrib_equil), 
          equil.file)
  
  return(equil.file)
  
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
                                     sim.type="mu",
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
    
    sim_equilibrium.in=rbind(sim_equilibrium.partner.in$distrib_equil,
                             sim_equilibrium.species.in$distrib_equil)
    reached_equil=(sim_equilibrium.partner.in$reached_equil&
                     sim_equilibrium.species.in$reached_equil)
    
  }else{
    sim_equilibrium.in= readRDS(sim_equil[id_forest])$distrib_equil
    reached_equil=readRDS(sim_equil[id_forest])$reached_equil
  }
  
  # set the appropriate name for variable selection
  if(sim.type=="mu"){id.obj="id.species.mu.obj"}else{id.obj="id.species.obj"}
  

  # Only make the simulation if population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      id.species.obj=species_list[species_list$ID.spclim==clim &
                                    species_list$species==sp &
                                    species_list$species_combination==species.in[i],
                                  id.obj][[1]]
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
    if(sim.type=="mu"){
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                        "waib", "wai2", "sgdd2", 
                                                                        "PC1", "PC2", "N", "SDM")],
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }else{
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }
    sim.in<-sim.in |> filter(var%in%c("N","BAsp"))
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

#' Function to make a list of invasion simulations
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param harv_rules.ref rules for harvesting
#' @param sim_equil path of simulation until equilibrium
#' @param threshold_pop dimater threshold of initial population
#' @param id_forest id of the forest to simulate
make_simulations_invasion_2 = function(species.combination,
                                     species_list,
                                     species_object,
                                     harv_rules.ref,
                                     sim_equil,
                                     BA_target=1,
                                     threshold_pop=200,
                                     sim.type="mu",
                                     id_forest){
  print(id_forest)
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equil[id_forest])$distrib_equil
  reached_equil = readRDS(sim_equil[id_forest])$reached_equil
  
  if(sim.type=="mu"){id.obj="id.species.mu.obj"}else{id.obj="id.species.obj"}
  
  n_lag<- sim_equilibrium.in %>%
    filter(var == "n", equil, species == s_p,size==0) %>%  #species.in[i]
    summarise(n_tot=sum(value)) |> 
    pull(n_tot)
  
  
  if(reached_equil){
    for(i in 1:length(species.in)){
      
      id.species.obj=species_list[species_list$ID.spclim==clim &
                                    species_list$species==sp &
                                    species_list$species_combination==species.in[i],
                                  id.obj][[1]]
      # Identify the file in species containing species i
      species.file.i = species_object[id.species.obj]
      # Store the file in the list
      list.species[[i]] = readRDS(species.file.i)
      
      lag.i = sim_equilibrium.in %>%
        filter(var == "n", equil, species == species.in[i],size==0) %>%
        mutate(n_tot=sum(value),
               value_2=(value*n_lag)/n_tot,
               n_tot_2=sum(value_2)) %>% 
        pull(value_2)
      rec.i = sim_equilibrium.in %>%
        filter(var == "n", equil, species == species.in[i], size>0) %>%  #species.in[i]
        mutate(value=case_when(size>threshold_pop~0,
                               TRUE~value)) |>
        mutate(BAtot=sum((size/2000)^2*pi*value),
               value_2=(value*BA_target)/BAtot,
               BAtot=sum((size/2000)^2*pi*value_2)) |> 
        pull(value_2)
      equil.i=c(lag.i,rec.i)
      equil.i[is.nan(equil.i)]<-0
   
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
    }
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    if(sim.type=="mu"){
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                        "waib", "wai2", "sgdd2", 
                                                                        "PC1", "PC2", "N", "SDM")],
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }else{
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }
    sim.in<-sim.in |> filter(var%in%c("N","BAsp"))
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

#' Function to make a list of invasion simulations
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param harv_rules.ref rules for harvesting
#' @param sim_equil path of simulation until equilibrium
#' @param threshold_pop dimater threshold of initial population
#' @param id_forest id of the forest to simulate
make_simulations_invasion_3 = function(species.combination,
                                       species_list,
                                       species_object,
                                       harv_rules.ref,
                                       sim_equil,
                                       BA_target=1,
                                       threshold_pop=200,
                                       sim.type="mu",
                                       id_forest){
  print(id_forest)
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in

  if(sim.type=="mu"){id.obj="id.species.mu.obj"}else{id.obj="id.species.obj"}

  for(i in 1:length(species.in)){
    
    id.species.obj=species_list[species_list$ID.spclim==clim &
                                  species_list$species==sp &
                                  species_list$species_combination==species.in[i],
                                id.obj][[1]]
    # Identify the file in species containing species i
    species.file.i = species_object[id.species.obj]
    # Store the file in the list
    list.species[[i]] = readRDS(species.file.i)
    
    mesh.i=list.species[[i]]$IPM$mesh
    
    rank.45 <- findInterval(450, mesh.i)
    
    coef=(0.03*0.5)/(pi*(0.45/2)^2)

    distrib.i=mesh.i
    distrib.i[]<-0
    distrib.i[rank.45]<-1/(pi*(0.45/2)^2) # value to get 0.5 BA/ha
    
    
    
    
    # Initiate the population at equilibrium
    list.species[[i]]$init_pop <- def_init_k(distrib.i)
    
    # Update disturbance function
    list.species[[i]]$disturb_fun <- disturb_fun
  }
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    if(sim.type=="mu"){
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                        "waib", "wai2", "sgdd2", 
                                                                        "PC1", "PC2", "N", "SDM")],
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }else{
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }
    sim.in<-sim.in |> filter(var%in%c("N","BAsp"))
  } 
  
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_invasion/", species.comb, ".rds")
  # Save simulation in a rdata
  create_dir_if_needed(forest.file)
  saveRDS(sim.in, forest.file)
  
  # Return output list
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
make_simulations_invasion_elast = function(sim_forest_list_elast,
                                           species_list_select_elast,
                                           species_object_mu,
                                           species_object_mu_elast,
                                           harv_rules.ref,
                                           sim_equil,
                                           sim_equil_elast,
                                           threshold_pop=0,
                                           sim.type="mu",
                                           id_forest){
  species.combination=sim_forest_list_elast$list.forests
  sim_eq_elast=sim_forest_list_elast$id.simul_forest
  
  print(id_forest)
  sp=species.combination[id_forest,"elast"][[1]]
  s_p=gsub(" ","_",species.combination[id_forest,"species"][[1]])
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  
  ## load equilibrium distributions
  if(length(species.in)>1){
    partner.comb=if_else(sub(paste0(sp,"\\."),"",species.comb)==species.comb,
                         if_else(sub(paste0("\\.",sp),"",species.comb)==species.comb,
                                 "Problem",
                                 sub(paste0("\\.",sp),"",species.comb)),
                         sub(paste0(sp,"\\."),"",species.comb))
    simul_eq.partner=species.combination |> 
      filter(elast==sp &
               species_combination==partner.comb &
               ID.spclim == clim) |> 
      pull(simul_eq)
    simul_eq.species=species.combination |> 
      filter(elast==sp &
               species_combination==sp &
               ID.spclim == clim) |> 
      pull(simul_eq_elast)
    simul_eq.species=match(simul_eq.species,sim_eq_elast) # find simulation in list
    # Read the simulation at equilibrium
    sim_equilibrium.partner.in = readRDS(sim_equil[simul_eq.partner])
    sim_equilibrium.species.in = readRDS(sim_equil_elast[simul_eq.species])
    
    sim_equilibrium.in=rbind(sim_equilibrium.partner.in$distrib_equil,
                             sim_equilibrium.species.in$distrib_equil)
    reached_equil=(sim_equilibrium.partner.in$reached_equil&
                     sim_equilibrium.species.in$reached_equil)
    
  }else{
    simul_eq.species=match(id_forest,sim_eq_elast)
    sim_equilibrium.in= readRDS(sim_equil_elast[simul_eq.species])$distrib_equil
    reached_equil=readRDS(sim_equil_elast[simul_eq.species])$reached_equil
  }

  # Only make the simulation if population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      if(species.in[i]==sp){
        id.species.obj=species_list_select_elast[species_list_select_elast$ID.spclim==clim &
                                                   species_list_select_elast$elast==sp &
                                                   species_list_select_elast$species_combination==species.in[i],
                                                 "id.species.elast.obj"][[1]]
        # Identify the file in species containing species i
        species.file.i = species_object_mu_elast[id.species.obj]
        # Store the file in the list
        list.species[[i]] = readRDS(species.file.i)
        
        # Extract the equilibrium for species i
        equil.i = sim_equilibrium.in %>%
          filter(var == "n", equil, species == s_p) %>%  
          mutate(value=case_when(size>threshold_pop~0,
                                 TRUE~value)) |>
          mutate(BAtot=sum((size/2000)^2*pi*value)) |> 
          pull(value)
        # Initiate the population at equilibrium
        list.species[[i]]$init_pop <- def_init_k(equil.i)
        
        # Update disturbance function
        list.species[[i]]$disturb_fun <- disturb_fun
      }else{
        id.species.obj=species_list_select_elast[species_list_select_elast$ID.spclim==clim &
                                                   species_list_select_elast$elast==sp &
                                                   species_list_select_elast$species_combination==species.in[i],
                                    "id.species.mu.obj"][[1]]
        # Identify the file in species containing species i
        species.file.i = species_object_mu[id.species.obj]
        # Store the file in the list
        list.species[[i]] = readRDS(species.file.i)
        
        # Extract the equilibrium for species i
        equil.i = sim_equilibrium.in %>%
          filter(var == "n", equil, species == species.in[[i]]) %>% 
          pull(value)
        
        # Initiate the population at equilibrium
        list.species[[i]]$init_pop <- def_init_k(equil.i)
        
        # Update disturbance function
        list.species[[i]]$disturb_fun <- disturb_fun
      }

    }
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
    if(sim.type=="mu"){
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                        "waib", "wai2", "sgdd2", 
                                                                        "PC1", "PC2", "N", "SDM")],
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }else{
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 1000,
                                equil_time = 1000, 
                                equil_dist = 50, 
                                equil_diff = 0.5, 
                                harvest = "default", 
                                SurfEch = 0.03,
                                verbose = TRUE)
    }
    sim.in<-sim.in |> filter(var%in%c("N","BAsp"))
  } else {
    sim.in = matrix()
  }
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_invasion/",
                     gsub(":","x", species.comb), #because ":" are not supported in filenames
                     ".rds")
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
                                     sim.type="mu",
                                     id_forest){
  sp=species.combination[id_forest,"species"][[1]]
  s_p=gsub(" ","_",sp)
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equil[id_forest])$distrib_equil
  reached_equil = readRDS(sim_equil[id_forest])$reached_equil
  
  # set the appropriate name for variable selection
  if(sim.type=="mu"){id.obj="id.species.mu.obj"}else{id.obj="id.species.obj"}
  
  sp_excl=(species.combination$excluded[id_forest]!="excluded")
  # Only make the simulation if population reached an equilibrium
  if(reached_equil&sp_excl){
    # Loop on all species
    for(i in 1:length(species.in)){
      
      id.species.obj=species_list[species_list$ID.spclim==clim &
                                    species_list$species==sp &
                                    species_list$species_combination==species.in[i],
                                  id.obj][[1]]
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
    # Run simulation till equilibrium
    if(sim.type=="mu"){
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 4000,
                                climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                        "waib", "wai2", "sgdd2", 
                                                                        "PC1", "PC2", "N", "SDM")],
                                equil_time = 4000, 
                                disturbance = disturbance.df_storm,
                                verbose = TRUE)
    }else{
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 4000,
                                equil_time = 4000, 
                                disturbance = disturbance.df_storm,
                                verbose = TRUE)
    }
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




#' Function to make a list of simulations with disturbance
#' @param species.combination table with all species combi for each climate cat
#' @param species_list table with all species fit info 
#' @param species_object list of files of fitted species
#' @param harv_rules.ref rules for harvesting
#' @param sim_equil path of simulation until equilibrium
#' @param disturb_coef.in table of disturbance coef per species
#' @param disturbance.df_storm disturbance characteristics
#' @param id_forest id of the forest to simulate
# sim_forest_list_elast,
# species_list_select_elast,
# species_object_mu,
# species_object_mu_elast,
# harv_rules.ref,
# sim_equil,
# sim_equil_elast,
# threshold_pop=0,
# sim.type="mu",
# id_forest
make_simulations_disturbance_elast = function(sim_forest_list_elast,
                                              species_list_select_elast,
                                              species_object_mu_elast,
                                              species_object_mu,
                                              harv_rules.ref,
                                              sim_equil_elast,
                                              disturb_coef.in,
                                              disturbance.df_storm,
                                              fit.list.allspecies,
                                              id_forest){
  species.combination=sim_forest_list_elast$list.forests
  sim_eq_elast=sim_forest_list_elast$id.simul_forest
  
  print(id_forest)
  sp=species.combination[id_forest,"elast"][[1]]
  s_p=gsub(" ","_",species.combination[id_forest,"species"][[1]])
  species.comb=species.combination[id_forest,"species_combination"][[1]]
  species.in=unlist(strsplit(species.comb,"\\."))
  clim=species.combination[id_forest,"ID.spclim"][[1]]
  id_simu=match(id_forest,sim_eq_elast)
  
  list.species <- vector("list", length(species.in))
  names(list.species) = species.in
  
  # Read the simulation at equilibrium
  sim_equilibrium.in = readRDS(sim_equil_elast[id_simu])$distrib_equil
  reached_equil = readRDS(sim_equil_elast[id_simu])$reached_equil
  
  # Only make the simulation if population reached an equilibrium
  if(reached_equil){
    # Loop on all species
    for(i in 1:length(species.in)){
      if(species.in[i]==sp){
        id.species.obj=species_list_select_elast[species_list_select_elast$ID.spclim==clim &
                                                   species_list_select_elast$elast==sp &
                                                   species_list_select_elast$species_combination==species.in[i],
                                                 "id.species.elast.obj"][[1]]
        # Identify the file in species containing species i
        species.file.i = species_object_mu_elast[id.species.obj]
        # Store the file in the list
        list.species[[i]] = readRDS(species.file.i)
        
        equil.i = sim_equilibrium.in %>%
          filter(var == "n", equil, species == s_p) %>% 
          pull(value)
        
        # Add disturbance coefficients
        list.species[[i]]$disturb_coef <- filter(disturb_coef.in, 
                                                 species == s_p)
      }else{
        id.species.obj=species_list_select_elast[species_list_select_elast$ID.spclim==clim &
                                                   species_list_select_elast$elast==sp &
                                                   species_list_select_elast$species_combination==species.in[i],
                                                 "id.species.mu.obj"][[1]]
        # Identify the file in species containing species i
        species.file.i = species_object_mu[id.species.obj]
        # Store the file in the list
        list.species[[i]] = readRDS(species.file.i)
        
        # Extract the equilibrium for species i
        equil.i =  sim_equilibrium.in %>%
          filter(var == "n", equil, species == species.in[i]) %>% 
          pull(value)
        # Add disturbance coefficients
        list.species[[i]]$disturb_coef <- filter(disturb_coef.in, 
                                                 species == species.in[i])
      
      }
      
      # Initiate the population at equilibrium
      list.species[[i]]$init_pop <- def_init_k(equil.i)
      
      # Update disturbance function
      list.species[[i]]$disturb_fun <- disturb_fun
      
    }
    
    
    # Make forest
    forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)
    
    # Run simulation till equilibrium
      sim.in = sim_deter_forest(forest.in, 
                                tlim = 4000,
                                climate=species.combination[id_forest,c("sgdd", "wai", "sgddb",
                                                                        "waib", "wai2", "sgdd2", 
                                                                        "PC1", "PC2", "N", "SDM")],
                                equil_time = 4000, 
                                disturbance = disturbance.df_storm,
                                verbose = TRUE)
      sim.short<-sim.in %>% filter(var%in%c("N","BAsp"))
  } else {
    sim.in = matrix()
  }
  
  forest.file=paste0("rds/", s_p, "/clim_", clim,
                     "/sim_disturbance/",
                     gsub(":","x", species.comb), #because ":" are not supported in filenames
                     ".rds")
  # Save simulation in a rdata
  create_dir_if_needed(forest.file)
  saveRDS(sim.short, forest.file)
  
  
  ## Compute resistance indicators
  out=species.combination |> 
    filter(simul_eq_elast == id_forest) |> 
    mutate(resistance = NA_real_, recovery = NA_real_, resilience = NA_real_, 
           t0 = NA_real_, thalf = NA_real_, SD = NA_real_, BA_diff = NA_real_, 
           BA_eq = NA_real_, dbh_mean = NA_real_, dbh_q10 = NA_real_, 
           dbh_q90 = NA_real_, dbh_mean_postdist = NA_real_, 
           dbh_q10_postdist = NA_real_, dbh_q90_postdist = NA_real_)
    
  # Identify disturbance time
  tdist = min(disturbance.df_storm$t)
  fit.species.i=fit.list.allspecies[[s_p]]
      
  if(!is.na(sim.in[1, 1])){
    # Read simulation i
    sim.i = sim.in |> 
      filter(species==s_p)
    
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
    out$dbh_mean<- subset(dbh_i, time == 1)$mean_dbh
    out$dbh_q10 <- subset(dbh_i, time == 1)$q10_dbh
    out$dbh_q90<- subset(dbh_i, time == 1)$q90_dbh
    out$dbh_mean_postdist <- subset(dbh_i, time != 1)$mean_dbh
    out$dbh_q10_postdist <- subset(dbh_i, time != 1)$q10_dbh
    out$dbh_q90_postdist <- subset(dbh_i, time != 1)$q90_dbh
    
    # Format the output
    data.i <- sim.i %>%
      filter(var == "BAsp") %>%
      filter(!equil) %>%
      group_by(time) %>%
      summarize(BA = sum(value))
    
    ## Calculate stability before disturbance (to check equilibrium)
    out$SD= sd(subset(data.i, time %in% c(1:(tdist-1)))$BA)
    out$BA_diff = diff(range(subset(data.i, time %in% c(1:(tdist-1)))$BA))
    
    ## Calculate resistance
    #  - Basal area at equilibrium
    Beq.i = mean((data.i %>% filter(time < min(disturbance.df_storm$t)))$BA)
    out$BA_eq = Beq.i
    # - Basal area after disturbance
    Bdist.i = (data.i %>% filter(time == max(disturbance.df_storm$t)+1))$BA
    # - Resistance : logit of the percentage of basal area that survived 
    #out$resistance[i] = Beq.i/(Beq.i - Bdist.i)
    out$resistance = log((Bdist.i/Beq.i)/(1 - (Bdist.i/Beq.i)))
    
    ## Calculate recovery
    #  - Time at which population recovered fully
    Rec.time.i = min((data.i %>% 
                        filter(time > max(disturbance.df_storm$t)) %>%
                        filter(BA > Beq.i))$time)
    # - Basal area 20 years after disturbance
    Bdist20.i = (data.i %>% filter(time == max(disturbance.df_storm$t)+21))$BA
    # - Recovery = slope of BA increase in teh 20 years after disturbance
    out$recovery = abs(Bdist20.i - Bdist.i)/20
    
    ## Calculate resilience
    out$resilience <- 1/sum((data.i %>%
                                  mutate(BA0 = .[which(.$time == 1), "BA"]) %>%
                                  mutate(diff = abs(BA - BA0)))$diff)
    
    ## Calculate t0
    #  - Time at which population recovered to 5% of the basal area lost
    Rec.0.time.i = min((data.i %>% 
                          filter(time > max(disturbance.df_storm$t)) %>%
                          filter(BA > (Beq.i + 19*Bdist.i)/20))$time)
    # - Recovery = time to recover minus time of disturbance
    out$t0 = Rec.0.time.i - max(disturbance.df_storm$t)
    
    ## Calculate thalf
    #  - Time at which population recovered to 50% of the basal area lost
    Rec.half.time.i = min((data.i %>% 
                             filter(time > max(disturbance.df_storm$t)) %>%
                             filter(BA > (Beq.i + Bdist.i)/2))$time)
    # - Recovery = time to recover minus time of disturbance
    out$thalf = Rec.half.time.i - max(disturbance.df_storm$t)
    
    
  }
      
  # Return output list
  return(list(forest.file=forest.file,
              out=out))
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

# 
# tar_load(c(species.list.ipm,
#            climate.cat))
# species_object=tar_read(species_object_mu)
# species_list=tar_read(species_list.select)
# tar_load(c(fit.list.allspecies,
#            harv_rules.ref))
# sp_id=7
make_mean_simul<-function(species.list.ipm,
                          climate.cat,
                          species_object,
                          species_list,
                          fit.list.allspecies,
                          harv_rules.ref,
                          sp_id){
  print(sp_id)
  s_p=species.list.ipm[sp_id]
  sp=gsub("_"," ",s_p)
  
  # get mean climate
  climate.cat=climate.cat$species.cat
  clim=apply(climate.cat[climate.cat$species==sp&
                     climate.cat$clim_id%in%c(5,6),c("wai","sgdd")],
             MARGIN=2,
             mean)
  clim=data.frame(sgdd=clim[["sgdd"]],
                  wai=clim[["wai"]]) %>% 
    mutate(sgddb=1/sgdd,
           waib=1/(1+wai),
           wai2=wai^2,
           sgdd2=sgdd^2
           )
  # 
  # species.combination[id_forest,c("sgdd", "wai", "sgddb",
  #                                 "waib", "wai2", "sgdd2", 
  #                                 "PC1", "PC2", "N", "SDM")]
  # get species mu
  id_obj=species_list %>% 
    filter(species_combination==s_p) %>% 
    dplyr::select(id.species.mu.obj) %>% 
    unique() %>% pull(id.species.mu.obj)
  
  species_mu<-readRDS(species_object[id_obj])
  list.species<- list(species_mu)
  names(list.species)<-s_p

  # create forest
  forest.in = new_forest(species = list.species,
                         harv_rules = harv_rules.ref)
  
  # launch equil sim
  sim.in = sim_deter_forest(forest.in, 
                            tlim = 4000,
                            climate=clim,
                            equil_time = 50000, 
                            equil_dist = 2000, 
                            equil_diff = 0.5, 
                            harvest = "default", 
                            SurfEch = 0.03,
                            verbose = TRUE)
  
  reached_equil = ifelse(
    is.na(sum((sim.in %>%
                 filter(var == "BAsp") %>%
                 filter(time == max(.$time) - 1))$value)), 
    FALSE, TRUE
  )
  distrib_equil = sim.in %>%
    filter(var == "n", equil)
  ba_equil<-sim.in %>%
    filter(var == "BAsp", equil)
  
  # launch invasion sim
  lag.i = distrib_equil %>%
    filter(size==0) %>%
    pull(value)
  rec.i = distrib_equil %>%
    filter(size>0) %>%
    mutate(value=case_when(size>200~0,
                           TRUE~value)) |>
    mutate(BAtot=sum((size/2000)^2*pi*value),
           value_2=(value*1)/BAtot,
           BAtot=sum((size/2000)^2*pi*value_2)) |> 
    pull(value_2)
  equil.i=c(lag.i,rec.i)
  equil.i[is.nan(equil.i)]<-0
  
  
  ## Initiate the population at equilibrium
  list.species[[1]]$init_pop <- def_init_k(equil.i)
  
  
  forest.inv = new_forest(species = list.species, harv_rules = harv_rules.ref)
  
  sim.inv = sim_deter_forest(forest.inv, 
                            tlim = 100,
                            climate=clim,
                            equil_time = 1000, 
                            equil_dist = 50, 
                            equil_diff = 0.5, 
                            harvest = "default", 
                            SurfEch = 0.03,
                            verbose = TRUE)
  delay.i<-as.numeric(fit.list.allspecies[[s_p]]$info[["delay"]])
  derivative.i <- sim.inv |> 
    filter(var=="BAsp",time>delay.i) |> 
    mutate(der=(value-lag(value))/(time-lag(time)),
           der2=(der-lag(der))/(time-lag(time))) %>% 
    filter(time<delay.i+50) |> 
    summarise(inv_50=mean(der,na.rm=TRUE))
  
  out<-data.frame(species=sp,
             ba_equil=ba_equil$value,
             inv_50=derivative.i$inv_50) %>% 
    cbind(clim)
  return(out)
}
  
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### POST SIMULATION ANALYSIS ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
           inv_50=NA_real_,
           BA_100=NA_real_,
           BA_500=NA_real_,
           BA_1000=NA_real_)
  for (i in 1:length(sim_invasion)){
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_invasion)))
    
    
    # Read simulation i
    sim.i = readRDS(sim_invasion[i])
    id.sim.i=id.simul_forest[i]
    species.i=species.combination[id.sim.i,"species"][[1]]
    s_p.i=gsub(" ","_",species.i)

      if(dim(sim.i)[1]>1){ # check simulation was performed (i.e. equil was reached)
        # get simulation index
        fit.species.i=fit.list.allspecies[[s_p.i]]
        delay.i=as.numeric(fit.species.i$info[["delay"]])
        
        BA_t <- sim.i %>% 
          filter(time%in%c(100,500,1000),!equil,var=="BAsp") %>% 
          filter(species==s_p.i) %>% 
          pull(value)
        derivative.i <- sim.i |> 
          filter(species==s_p.i,var=="BAsp",time>delay.i) |> 
          mutate(der=(value-lag(value))/(time-lag(time)),
                 der2=(der-lag(der))/(time-lag(time)))
        
        first.extr<-derivative.i |>  
          filter(time<1000) |> 
          mutate(sign_eq=(sign(der)==sign(lag(der)))) |> 
          filter(sign_eq==FALSE) |> slice(1) |> 
          dplyr::select(time) |> pull(time)
        
        first.extr=ifelse(length(first.extr)==0,50,first.extr)
        
        inv.1<-derivative.i |> 
          filter(time<(first.extr-5)) |> 
          summarise(inv_mean=mean(der,na.rm=TRUE),
                    inv_max=max(der,na.rm=TRUE))
        
        inv.2<-derivative.i |> 
          filter(time<delay.i+50) |> 
          summarise(inv_50=mean(der,na.rm=TRUE))
        out[out$simul_eq==id.sim.i,c("inv_mean","inv_max","inv_50",
                                     "BA_100","BA_500","BA_1000")]=
          as.list(c(inv.1$inv_mean[[1]],
               inv.1$inv_max[[1]],
               inv.2$inv_50[[1]],
               BA_t))
      }
  }
  
  return(out)
}


#' Function to compute invasion rate for a simulation
#' @param species.combination table with all species combi for each climate cat
#' @param sim_invasion table with all species fit info 
#' @param id_forest id of the forest to simulate
#' @param fit.list.allspecies list of species with ipm objects, here to get the delay
get_invasion_rate_2<-function(species.combination,
                            id.simul_forest,
                            sim_invasion,
                            fit.list.allspecies){
  out=species.combination |> 
    filter(simul_eq %in% id.simul_forest) |> 
    mutate(inv_mean=NA_real_,
           inv_max=NA_real_,
           inv_50=NA_real_,
           BA_100=NA_real_,
           BA_500=NA_real_,
           BA_1000=NA_real_)
  for (i in 1:length(sim_invasion)){
    # Printer
    print(paste0("Reading simulation ", i, "/", length(sim_invasion)))
    
    
    # Read simulation i
    sim.i = readRDS(sim_invasion[i])
    id.sim.i=id.simul_forest[i]
    species.i=species.combination[id.sim.i,"species"][[1]]
    s_p.i=gsub(" ","_",species.i)
    
    if(dim(sim.i)[1]>1){ # check simulation was performed (i.e. equil was reached)
      # get simulation index
      fit.species.i=fit.list.allspecies[[s_p.i]]
      delay.i=as.numeric(fit.species.i$info[["delay"]])
      
      BA_t <- sim.i %>% 
        filter(time%in%c(100,500,1000),!equil,var=="BAsp") %>% 
        filter(species==s_p.i) %>% 
        pull(value)
      derivative.i <- sim.i |> 
        filter(species==s_p.i,var=="BAsp") |> 
        mutate(der=(value-lag(value))/(time-lag(time)),
               der2=(der-lag(der))/(time-lag(time)))
      
      first.extr<-derivative.i |>  
        filter(time<1000) |> 
        mutate(sign_eq=(sign(der)==sign(lag(der)))) |> 
        filter(sign_eq==FALSE) |> slice(1) |> 
        dplyr::select(time) |> pull(time)
      
      first.extr=ifelse(length(first.extr)==0,50,first.extr)
      
      inv.1<-derivative.i |> 
        filter(time<(first.extr-5)) |> 
        summarise(inv_mean=mean(der,na.rm=TRUE),
                  inv_max=max(der,na.rm=TRUE))
      
      inv.2<-derivative.i |> 
        filter(time<delay.i+50) |> 
        summarise(inv_50=mean(der,na.rm=TRUE))
      out[out$simul_eq==id.sim.i,c("inv_mean","inv_max","inv_50",
                                   "BA_100","BA_500","BA_1000")]=
        as.list(c(inv.1$inv_mean[[1]],
                  inv.1$inv_max[[1]],
                  inv.2$inv_50[[1]],
                  BA_t))
    }
  }
  
  return(out)
}


#' Function to compute invasion rate for a simulation
#' @param species.combination table with all species combi for each climate cat
#' @param sim_invasion table with all species fit info 
#' @param id_forest id of the forest to simulate
#' @param fit.list.allspecies list of species with ipm objects, here to get the delay
get_invasion_rate_elast<-function(species.combination,
                                  id.simul_forest,
                                  sim_invasion,
                                  fit.list.allspecies){
  out=species.combination |> 
    filter(simul_eq_elast %in% id.simul_forest) |> 
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
        filter(time<1000) |> 
        mutate(sign_eq=(sign(der2)==sign(lag(der2)))) |> 
        filter(sign_eq==FALSE) |> slice(1) |> 
        dplyr::select(time) |> pull(time)
      
      first.extr=ifelse(length(first.extr)==0,50,first.extr)
      
      inv.1<-derivative.i |> 
        filter(time<(first.extr-5)) |> 
        summarise(inv_mean=mean(der,na.rm=TRUE),
                  inv_max=max(der,na.rm=TRUE))
      
      inv.2<-derivative.i |> 
        filter(time<delay.i+50) |> 
        summarise(inv_50=mean(der,na.rm=TRUE))
      out[out$simul_eq_elast==id.sim.i,c("inv_mean","inv_max","inv_50")]=
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
                                 id.sim.i,
                                 sim_dist.id,
                                 sim_disturbance,
                                 fit.list.allspecies,
                                 disturbance.df_storm){
  
  out=species.combination |> 
    filter(simul_eq == id.sim.i) |> 
    mutate(resistance = NA_real_, recovery = NA_real_, resilience = NA_real_, 
           t0 = NA_real_, thalf = NA_real_, SD = NA_real_, BA_diff = NA_real_, 
           BA_eq = NA_real_, dbh_mean = NA_real_, dbh_q10 = NA_real_, 
           dbh_q90 = NA_real_, dbh_mean_postdist = NA_real_, 
           dbh_q10_postdist = NA_real_, dbh_q90_postdist = NA_real_)
  i=match(id.sim.i,sim_dist.id)
  
  # Identify disturbance time
  tdist = min(disturbance.df_storm$t)
  
  # Printer
  print(paste0("Reading simulation ", i, "/", length(sim_disturbance)))
    
  # get simulation index
  species.i=species.combination[id.sim.i,"species"][[1]]
  s_p.i=gsub(" ","_",species.i)
  fit.species.i=fit.list.allspecies[[s_p.i]]
    # delay.i=as.numeric(fit.species.i$info[["delay"]])
  
  if(species.combination$excluded[id.sim.i]!="excluded"){
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
        out$dbh_mean <- subset(dbh_i, time == 1)$mean_dbh
        out$dbh_q10 <- subset(dbh_i, time == 1)$q10_dbh
        out$dbh_q90 <- subset(dbh_i, time == 1)$q90_dbh
        out$dbh_mean_postdist <- subset(dbh_i, time != 1)$mean_dbh
        out$dbh_q10_postdist <- subset(dbh_i, time != 1)$q10_dbh
        out$dbh_q90_postdist <- subset(dbh_i, time != 1)$q90_dbh
        
        # Format the output
        data.i <- sim.i %>%
          filter(var == "BAsp") %>%
          filter(!equil) %>%
          group_by(time) %>%
          summarize(BA = sum(value))
        
        ## Calculate stability before disturbance (to check equilibrium)
        out$SD = sd(subset(data.i, time %in% c(1:(tdist-1)))$BA)
        out$BA_diff = diff(range(subset(data.i, time %in% c(1:(tdist-1)))$BA))
        
        ## Calculate resistance
        #  - Basal area at equilibrium
        Beq.i = mean((data.i %>% filter(time < min(disturbance.df_storm$t)))$BA)
        out$BA_eq = Beq.i
        # - Basal area after disturbance
        Bdist.i = (data.i %>% filter(time == max(disturbance.df_storm$t)+1))$BA
        # - Resistance : logit of the percentage of basal area that survived 
        #out$resistance[i] = Beq.i/(Beq.i - Bdist.i)
        out$resistance = log((Bdist.i/Beq.i)/(1 - (Bdist.i/Beq.i)))
        
        ## Calculate recovery
        #  - Time at which population recovered fully
        Rec.time.i = min((data.i %>% 
                            filter(time > max(disturbance.df_storm$t)) %>%
                            filter(BA > Beq.i))$time)
        # - Basal area 20 years after disturbance
        Bdist20.i = (data.i %>% filter(time == max(disturbance.df_storm$t)+21))$BA
        # - Recovery = slope of BA increase in teh 20 years after disturbance
        out$recovery = abs(Bdist20.i - Bdist.i)/20
        
        ## Calculate resilience
        out$resilience <- 1/sum((data.i %>%
                                      mutate(BA0 = .[which(.$time == 1), "BA"]) %>%
                                      mutate(diff = abs(BA - BA0)))$diff)
        
        ## Calculate t0
        #  - Time at which population recovered to 5% of the basal area lost
        Rec.0.time.i = min((data.i %>% 
                              filter(time > max(disturbance.df_storm$t)) %>%
                              filter(BA > (Beq.i + 19*Bdist.i)/20))$time)
        # - Recovery = time to recover minus time of disturbance
        out$t0 = Rec.0.time.i - max(disturbance.df_storm$t)
        
        ## Calculate thalf
        #  - Time at which population recovered to 50% of the basal area lost
        Rec.half.time.i = min((data.i %>% 
                                 filter(time > max(disturbance.df_storm$t)) %>%
                                 filter(BA > (Beq.i + Bdist.i)/2))$time)
        # - Recovery = time to recover minus time of disturbance
        out$thalf = Rec.half.time.i - max(disturbance.df_storm$t)
        
        
      }
    }

  return(out)
}


#' Function to compute initial basal area of simulation
#' @description
#' The function takes in argument the list of forest to simulate until equilibrium,
#' and the index of true forests. The latter corresponds to forests
#' 
get_bainit_dist<-function(sim_forest_list,
                          sim_equil,
                          disturbance_metric,
                          elast=TRUE,
                          species.list.ipm){
  forest.list=sim_forest_list$list.forests
  disturbance_metric[species.list.ipm]<-NA
  
  for(id_forest in 1:dim(disturbance_metric)[1]){
    sp=disturbance_metric[id_forest,ifelse(elast,"elast","species")][[1]]
    s_p=gsub(" ","_",sp)
    clim=disturbance_metric[id_forest,"ID.spclim"][[1]]
    species.comb=disturbance_metric[id_forest,"species_combination"][[1]]
    if(elast){
      id_simu=forest.list |>
        filter(elast==sp &
                 species_combination==species.comb &
                 ID.spclim == clim) |> 
          pull(simul_eq_elast)
      simul_eq.partner=match(id_simu,sim_forest_list$id.simul_forest)
      
        }else{
        simul_eq.partner=forest.list |> 
          filter(species==sp &
                   species_combination==species.comb &
                   ID.spclim == clim) |> 
          pull(simul_eq)
      }
      file=sim_equil[simul_eq.partner]
      equil=readRDS(file)
      if(equil$reached_equil){
        sp=disturbance_metric[id_forest,"species"][[1]]
        ba_equil=equil$distrib_equil %>% 
          mutate(basp=pi*(size/(2000))^2) %>% 
          group_by(species) %>% 
          summarise(BAtot=sum(basp*value))
        for(ssp in ba_equil$species){
          disturbance_metric[id_forest,ssp]<-ba_equil[ba_equil$species==ssp,"BAtot"]
        }
      }
  }
  return(disturbance_metric)

}


#' Function to compute initial basal area of simulation for invasion rate
#' @description
#' The function takes in argument the list of forest to simulate until equilibrium,
#' and the index of true forests. The latter corresponds to forests
#' @param sim_forest_list_elast list of forests to simulate
#' @param sim_equil list of sim file until equilibrium
#' @param invasion_metric_elast table with invasion metrics
#' @param elast TRUE/FALSE whether it concerns elasticity analysis or not
#' @param species.list.ipm list of species with models
get_bainit_inv<-function(sim_forest_list_elast,
                         sim_equil,
                         invasion_metric_elast,
                         elast=TRUE,
                         species.list.ipm){
  forest.list=sim_forest_list_elast$list.forests
  invasion_metric_elast[species.list.ipm]<-NA
  
  for(id_forest in 1:dim(invasion_metric_elast)[1]){
    sp=invasion_metric_elast[id_forest,ifelse(elast,"elast","species")][[1]]
    s_p=gsub(" ","_",sp)
    clim=invasion_metric_elast[id_forest,"ID.spclim"][[1]]
    species.comb=invasion_metric_elast[id_forest,"species_combination"][[1]]
    if(length(strsplit(species.comb,"\\.")[[1]])>1){
      partner.comb=if_else(sub(paste0(s_p,"\\."),"",species.comb)==species.comb,
                           if_else(sub(paste0("\\.",s_p),"",species.comb)==species.comb,
                                   "Problem",
                                   sub(paste0("\\.",s_p),"",species.comb)),
                           sub(paste0(s_p,"\\."),"",species.comb))
      if(elast){
        simul_eq.partner=forest.list |> 
          filter(elast==sp &
                   species_combination==partner.comb &
                   ID.spclim == clim) |> 
          pull(simul_eq)
      }else{
        simul_eq.partner=forest.list |> 
          filter(species==sp &
                   species_combination==partner.comb &
                   ID.spclim == clim) |> 
          pull(simul_eq)
      }
      file=sim_equil[simul_eq.partner]
      equil=readRDS(file)
      if(equil$reached_equil){
        sp=invasion_metric_elast[id_forest,"species"][[1]]
        ba_equil=equil$distrib_equil %>% 
          mutate(basp=pi*(size/(2000))^2) %>% 
          group_by(species) %>% 
          summarise(BAtot=sum(basp*value))
        for(ssp in ba_equil$species){
          invasion_metric_elast[id_forest,ssp]<-ba_equil[ba_equil$species==ssp,"BAtot"]
        }
      }
    }
  }
  return(invasion_metric_elast)
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


#' Function to get basal area difference at equil simulations

get_dif_ba<-function(disturbance_ba,
                     species.list.ipm){
  ba_dif<-disturbance_ba %>% 
    pivot_longer(cols=matches(species.list.ipm)) %>% 
    group_by(species,species_combination,clim_id) %>% 
    mutate(ba_multiple=case_when(name==gsub(" ","_",species)~1,
                                 TRUE~0)) %>%
    summarize(ba_partner=sum(value*abs(1-ba_multiple),na.rm = TRUE),
              ba_target=sum(value*ba_multiple,na.rm = TRUE)) %>%
    rowwise() %>% 
    mutate(n_species=length(strsplit(species_combination,"\\.")[[1]])) %>% 
    ungroup() %>% 
    group_by(species,clim_id) %>% 
    arrange(species,clim_id,n_species) %>% 
    mutate(ba_dif=ba_target/ba_target[1]) %>% 
    ungroup()
  return(ba_dif)
}

get_list_pars<-function(fit.list.allspecies){
  list_pars_sp=vector("list",length = length(fit.list.allspecies))
  for(sp in 1:length(fit.list.allspecies)){
    sp_fit=fit.list.allspecies[[sp]]
    names(list_pars_sp)[sp]=names(fit.list.allspecies)[sp]
    n_pars=length(sp_fit$sv$params_m)+length(sp_fit$gr$params_m)+
      length(sp_fit$rec$params_m)
    list_par=c()
    for (vr in c("sv","gr","rec")){
      list_par=c(list_par,paste0(names(list_pars_sp)[sp],"-",vr,"-",names(sp_fit[[vr]]$params_m)))
    } 
    list_pars_sp[[sp]]=list_par
  }
  list_pars_sp=sort(unique(unlist(list_pars_sp)))

  return(list_pars_sp)
}



#' Function to get indices of species in TRY
get_TRY_species<-function(file.try="data/try_species.csv",
                          species.list.ipm){
  try_sp<-read.csv2(file.try) |> 
    mutate(species=paste(AccSpeciesName,ObsNum,sep="_")) |> 
    filter(species%in%c(species.list.ipm,"Betula_pubescens","Betula_pendula")) |> 
    group_by(species) |>
    filter(as.numeric(TraitNum)==max(as.numeric(TraitNum))) 
}




























