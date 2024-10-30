# draft 
tar_load(c(species.list.ipm))
species_object=tar_read(species_object_mu)
species_list=tar_read(species_list.select)
tar_load(c(fit.list.allspecies,harv_rules.ref))

sim_tot<-setNames(data.frame(matrix( nrow = 0,ncol=7)),
                  nm=c("species","var","time","mesh","size","equil","value"))
for(sp_id in 1:length(species.list.ipm)){
  print(paste0("Run species ",sp_id,"/",length(species.list.ipm)))
  tar_load(climate.cat)
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

  # get species mu
  id_obj=species_list %>% 
    filter(species_combination==s_p) %>% 
    dplyr::select(id.species.mu.obj) %>% 
    unique() %>% pull(id.species.mu.obj)
  
  species_mu<-readRDS(species_object[id_obj])
  list.species<- list(species_mu)
  names(list.species)<-s_p
  list.species[[1]]$init_pop <-  def_initBA(0.1)
  
  
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
  # delay=as.numeric(fit.list.allspecies[[s_p]]$info[["delay"]])
  sim_tot<-rbind(sim_tot,
                 sim.inv |> filter(var=="BAsp"))
  
}

# sim_fasy=sim.inv
sim_tot|> 
  ggplot(aes(time,value,color=species))+
  geom_line(size=1)+
  xlim(c(0,80))+ylim(c(0,25))+
  scale_y_log10()


inv_18 <- sim_tot |> 
  group_by(species) |> 
  mutate(der=(value-lag(value))/(time-lag(time)),
         der2=(der-lag(der))/(time-lag(time))) |> 
  filter(time<18) |> 
  summarise(inv_50=mean(der,na.rm=TRUE))
save(inv_18,file ="inv_18.RData" )
## get demo traits
tar_load(mean_demo)

