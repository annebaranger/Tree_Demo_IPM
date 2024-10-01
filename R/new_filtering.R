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
out<-setNames(data.frame(matrix(nrow=0,ncol=27)),
              nm=t)
for(s_p in species.list.disturbance){
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
    select(plotcode,clim_id, species_combination,n_species) |> unique() |> 
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
           prop_cum=cumsum(prop),
           prop_cum_filt=is.sp*cumsum(prop*is.sim),
           max_cum_filt=max(prop_cum_filt))   
  
  out_sp<-condi.init |>
    filter(species==sp) |> 
    left_join(species.combinations |> 
                filter((species_combination==s_p|prop_cum_filt<0.8)&is.sim),
              by="clim_id") |> 
    arrange(clim_id,n_species)
  out<-rbind(out,out_sp)
}


  # filter(is.sp==n_species) |> 
  # ungroup() |> 
  # arrange(clim_id, desc(n)) |>
  # mutate(pop_cum_filter=cumsum(prop))

