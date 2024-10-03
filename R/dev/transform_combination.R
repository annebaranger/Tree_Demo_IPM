species.combination |> 
  select(species,wai_id,sgdd_id,ID.spclim,clim_lab,wai,sgdd,sgdd2,sgddb,waib,PC1,PC2,N,SDM,species_combination) |> 
  mutate(species_combination=strsplit(species_combination,"\\.")) |> 
  unnest(cols = species_combination) |> 
  unique() |> 
  mutate(file.ipm=paste0("rds/",gsub(" ","_",species),"/clim_",ID.spclim,"/",species_combination,".rds")) |> View()


list.forests<-species.combination |> 
  select(species,species_combination,
         wai_id,sgdd_id,ID.spclim,clim_lab,
         wai,sgdd,sgdd2,sgddb,waib,PC1,PC2,N,SDM) 

list.forests.bis<-species.combination |> 
  select(species,species_combination,
         wai_id,sgdd_id,ID.spclim,clim_lab,
         wai,sgdd,sgdd2,sgddb,waib,PC1,PC2,N,SDM) |> 
  mutate(species_combination=strsplit(species_combination,"\\.")) |> 
  unnest(cols = species_combination) |> 
  unique() 


rbind(list.forests,list.forests.bis) |> unique() |> 
  arrange(species,ID.spclim) |> 
  mutate(forest.real=grepl(gsub(" ","_",species),species_combination)) |> 
  arrange(species,ID.spclim) |> 
  View()

