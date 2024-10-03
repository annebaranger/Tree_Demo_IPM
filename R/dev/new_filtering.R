library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)

tar_load(c(FUNDIV_data,
           species.list.disturbance,
           species.list.ipm))
climate.cat<-make_climate_cat_pca(FUNDIV_data,
                                  species.list.ipm)
condi.init=climate.cat$species.cat
FUNDIV_plotcat=climate.cat$FUNDIV_plotcat
nsp_per_richness=10
names_comp=c("species","clim_id","clim_low","clim_up","n_plot","wai","sgdd","sgdd2",
             "wai2","sgddb","waib","PC1","PC2","N","SDM","ID.spclim","clim_lab","ID.species",
             "species_combination","n_species","n","is.sim","is.sp","prop","prop_all",
             "prop_cum","prop_cum_all","prop_cum_filt","max_cum_filt")
out<-setNames(data.frame(matrix(nrow=0,ncol=29)),
              nm=names_comp)
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
           prop_all=n/sum(n),
           prop_cum=cumsum(prop),
           prop_cum_all=cumsum(prop_all),
           prop_cum_filt=is.sp*cumsum(prop*is.sim),
           max_cum_filt=max(prop_cum_filt))   
  species.combinations |> 
    group_by(species_combination) |> 
    summarise(ntot=sum(n)) |> 
    ungroup() |>
    mutate(prop=ntot/sum(ntot)) |> 
    filter(species_combination==s_p) |> pull(prop) |> 
    print()
  out_sp<-condi.init |>
    filter(species==sp) |> 
    left_join(species.combinations |> 
                filter((species_combination==s_p|prop_cum_filt<0.8)&is.sim),
              by="clim_id") |> 
    arrange(clim_id,n_species)
  out<-rbind(out,out_sp)
}


  names_comp=c("species","clim_id","clim_low","clim_up","n_plot","wai","sgdd","sgdd2",
        "wai2","sgddb","waib","PC1","PC2","N","SDM","ID.spclim","clim_lab","ID.species",
        "species_combination","n_species","n","is.sim","is.sp","is.select","prop",
        "prop_cum","prop_cum_filt","prop_select","max_cum_filt","max_cum_select")
comp<-setNames(data.frame(matrix(nrow=0,ncol=30)),
               nm=names_comp)
for(s_p in species.list.disturbance){
  sp=gsub("_"," ",s_p)
  print(sp)
  preset.combi<-tar_read(species.combination) |> 
    filter(species==sp) |> 
    select(clim_id,species_combination) |> 
    mutate(is.select=TRUE)
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
    left_join(preset.combi,by=c("clim_id","species_combination")) |> 
    mutate(is.select=if_else(is.na(is.select),FALSE,TRUE)) |> 
    mutate(prop=(n*is.sp)/sum(n*is.sp), # prop not accounting for species alone
           prop_cum=cumsum(prop),
           prop_cum_filt=is.sp*cumsum(prop*is.sim),
           prop_select=is.sp*cumsum(prop*is.select),
           max_cum_filt=max(prop_cum_filt),
           max_cum_select=max(prop_select))   
  
  out_sp<-condi.init |>
    filter(species==sp) |> 
    left_join(species.combinations |> 
                filter((species_combination==s_p|prop_cum_filt<0.8)&is.sim),
              by="clim_id") |> 
    arrange(clim_id,n_species)
  comp<-rbind(comp,out_sp)
}

# find principal competitors of species

competitor<-setNames(data.frame(matrix(nrow=0,ncol=5)),
               nm=c("species","sp_ratio_mean","BA_tot","is.sim","sp"))

for(s_p in species.list.ipm){
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
           sp=sp) |> 
    slice(1:20)
  competitor<-rbind(competitor,competitor_sp)
}

