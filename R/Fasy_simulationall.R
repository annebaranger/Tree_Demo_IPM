library(matreex)
library(dplyr)
library(tidyr)
library(ggplot2)
library(targets)


tar_load(c(species.list.ipm,harv_rules.ref,disturbance.df_storm,species.combination,climate.cat,species_clim))

# create forest list to be able tu run than in parallel
sp_id=14
s_p=species.list.ipm[sp_id]
sp=gsub("_"," ",s_p)
mu.file=species_clim[[sp_id]]

species.combination |> 
  filter(species==sp) |> 
  mutate(file.sim.equil = paste0("rds/", s_p, "/climate_", ID.spclim,
                                 "/sim_equilibrium/", species_combination, ".rds"),
         file.sim.dist = paste0("rds/", s_p, "/climate_", ID.spclim,
                                "/sim_disturbance/", species_combination, ".rds"),
         file.sim.inv = paste0("rds/", s_p, "/climate_", ID.spclim,
                                "/sim_invastion/", species_combination, ".rds"),
         file.sim.ibm = paste0("rds/", s_p, "/climate_", ID.spclim,
                                "/sim_ibm/", species_combination, ".rds"),
         species_combination=strsplit(species_combination,"\\.")) -> species.combination.list

sim_index=7

clim=species.combination.list[sim_index,c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", 
                                  "PC1", "PC2", "N", "SDM")]

species.in=unlist(species.combination.list[sim_index,
                                    "species_combination"][[1]])
list.species <- vector("list", length(species.in))
names(list.species) = species.in

for(i in 1:length(species.in)){
  
  # Identify the file in species containing species i
  species.file.i = species_clim[grep(species.in[i],species.list.ipm)]
  updated_file_path <- gsub("_mu\\.rds$", "_species.rds", species.file.i)
  # Store the file in the list
  list.species[[i]] = readRDS(updated_file_path)
  
}

# Make forest
forest.in = new_forest(species = list.species, harv_rules = harv_rules.ref)

# Run simulation till equilibrium
sim.in = sim_deter_forest(forest.in,
                          climate=clim,
                          tlim = 3000, equil_time = 50000, equil_dist = 2000, 
                          equil_diff = 0.5, harvest = "default",
                          correction = "cut",verbose = TRUE)

sim.in  %>%
  dplyr::filter(var == "BAsp", ! equil) %>%
  ggplot(aes(x = time, y = value, color = species)) +
  geom_line(linewidth = .4) + ylab("BA") +
  stat_summary(fun = "sum",  aes(col="Total"),
               geom ='line', linetype = "dashed", linewidth = .3)
sim.in = sim_indiv_forest(
  forest.in, tlim = 500,
  verbose = TRUE)

sim.in  %>%
  dplyr::filter(var == "N", ! equil) %>%
  ggplot(aes(x = time, y = value, color = species)) +
  geom_line(linewidth = .4) + ylab("N") 

sim.in  %>%
  filter(species=="Fagus_sylvatica") |> 
  dplyr::filter(var == "BAsp", ! equil) %>%
  filter(time<1000) |> 
  ggplot(aes(x = time, y = value, color = species)) +
  geom_line(linewidth = .4) + ylab("BA") 



# using ipm
species.in
IPM.1 = make_IPM(
  species = species.in[[1]], 
  climate = clim, 
  fit =  fit.list.allspecies[[species.in[[1]]]],
  clim_lab = "clim1",
  mesh = c(m = 700, L = 100, U = as.numeric(
    fit.list.allspecies[[species.in[[1]]]]$info[["max_dbh"]]) * 1.1),
  BA = 0:200, verbose = TRUE, correction = "none"
)
species.1=species(IPM=IPM.1,
                  init_pop = def_initBA(20),
                  harvest_fun = def_harv,
                  disturb_fun = def_disturb)
IPM.2 = make_IPM(
  species = species.in[[2]], 
  climate = clim, 
  fit =  fit.list.allspecies[[species.in[[2]]]],
  clim_lab = "clim1",
  mesh = c(m = 700, L = 100, U = as.numeric(
    fit.list.allspecies[[species.in[[1]]]]$info[["max_dbh"]]) * 1.1),
  BA = 0:200, verbose = TRUE, correction = "none"
)
species.2=species(IPM=IPM.2,
                  init_pop = def_initBA(20),
                  harvest_fun = def_harv,
                  disturb_fun = def_disturb)
list.species=list(species.1,species.2)
names(list.species)=species.in
forest.in.ipm = new_forest(species = list.species, harv_rules = harv_rules.ref)

sim.in.ipm = sim_deter_forest(
  forest.in.ipm, tlim = 4000, equil_time = 50000, equil_dist = 2000, 
  equil_diff = 0.5, harvest = "default", SurfEch = 0.03, verbose = TRUE)


# time required for simulations with 1 IPM per climate
list.species.combi=c()
for(i in 1:max(species.combination.list$ID.spclim)){
  species.combination.list |> 
    filter(ID.spclim==i) |> 
    pull(species_combination) |> unlist() |> unique()->speciesclim
  list.species.combi=c(list.species.combi,length(speciesclim))
}
dim(species.combination.list)[1]
(sum(list.species.combi)*240+dim(species.combination.list)[1]*45*3)


# time required for simulation with mu matrix
dim(species.combination.list)[1]*150*3+180
