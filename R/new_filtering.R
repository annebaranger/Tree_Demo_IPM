library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)

tar_load(c(FUNDIV_data,
           species.list.disturbance,
           species.list.ipm))
condi.init=tar_read(climate.cat)$species.cat
FUNDIV_plotcat=tar_read(climate.cat)$FUNDIV_plotcat


s_p=species.list.ipm[sp_id]
sp=gsub("_"," ",s_p)
print(sp)
species.combinations <- FUNDIV_data |>
  filter(dbh1>0) |>
  group_by(plotcode) |> 
  mutate(BAsum=sum(ba_ha1)) |> 
  group_by(plotcode,species) |> 
  mutate(BAsumsp=sum(ba_ha1), 
         sp_ratio=BAsumsp/BAsum) |> 
  # join with climate cat of the targetted species
  left_join(FUNDIV_plotcat |>
              filter(species==sp) |>
              dplyr::select(plotcode,clim_id,clim_low,clim_up), #wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up
            by="plotcode") |> 
  dplyr::select(treecode,plotcode,species,sp_ratio,clim_id,clim_low,clim_up) |> #wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up
  # filter(!is.na(clim_id)) %>% #wai_id
  # Group by the climatic condition and plot
  group_by(clim_id, plotcode) %>% #wai_id, sgdd_id
  # Summarize the species combination in each group, and species richness
  filter(sp%in%species)|>
  filter(!(species!=sp&sp_ratio<0.05)) |> 
  mutate(species_combination = paste(sort(unique(gsub(" ","_",species))), collapse = "."),
            n_species = n_distinct(species))  |> 
  ungroup() |> 
  #  count each unique combination's frequency within each wai_cat and sgdd_cat group
  count(clim_id, species_combination,n_species) |>  #wai_id, sgdd_id
  filter(n_species<nsp_per_richness) |>
  group_by(clim_id) |> #wai_id, sgdd_id
  mutate(prop=n/sum(n),
         prop_cum=cumsum(prop)) |>
  # Optionally, arrange the results for better readability
  arrange(clim_id, desc(n)) |> # wai_id, sgdd_id
  rowwise() |> 
  mutate(is.sp=sum(sapply(species.list.disturbance,function(x)grepl(x,species_combination)))) |> 
  filter(is.sp==n_species) |> 
  ungroup() |> 
  arrange(clim_id, desc(n)) |>
  mutate(pop_cum_filter=cumsum(prop))


species.target<- FUNDIV_plotcat |>
  filter(species==sp) |> 
  dplyr::select(clim_id) |> # wai_id,sgdd_id
  distinct() |> 
  mutate(species_combination=s_p) |> 
  left_join(species.combinations |>
              filter(species_combination==s_p),
            by=c("clim_id","species_combination")) |> #"wai_id","sgdd_id"
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
            by=c("clim_id")) |> #"wai_id","sgdd_id"
  arrange(clim_id) #wai_id,sgdd_id
