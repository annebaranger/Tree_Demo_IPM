# draft 
library(matreex)
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
  list.species[[1]]$init_pop <-  def_initBA(4)
  
  
  forest.inv = new_forest(species = list.species, harv_rules = harv_rules.ref)
  
  sim.inv = sim_deter_forest(forest.inv, 
                             tlim = 500,
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
  # filter(species%in%c("Pinus_sylvestris","Fagus_sylvatica","Prunus_padus","Quercus_suber")) |> 
  ggplot(aes(time,value,color=species))+
  geom_line(size=1)+
  # xlim(c(0,80))+ylim(c(0,25))+
  scale_y_log10()
delay_sp<-data.frame(species=names(fit.list.allspecies),
           delay=unlist(lapply(names(fit.list.allspecies),
                               function(x)as.numeric(fit.list.allspecies[[x]]$info[["delay"]]))))
sim_tot |> 
  group_by(species) |> 
  # filter(value>5) |> 
  mutate(time=time-min(time)) |> 
  filter(time<20) |> 
  ggplot(aes(time,value,color=species))+
  geom_line(size=1)+
  geom_hline(yintercept = 2)

inv_30<-sim_tot |> 
  group_by(species) |> 
  filter(value>5) |> 
  mutate(time=time-min(time)) |> 
  mutate(der=(value-lag(value))/(time-lag(time)),
         der2=(der-lag(der))/(time-lag(time))) |> 
  filter(time<30) |> 
  summarise(inv_30=mean(der,na.rm=TRUE)) 
quersube<-sim_tot |> 
  filter(species=="Quercus_suber") |> 
  filter(time>150,time<300) |> 
  mutate(der=(value-lag(value))/(time-lag(time)),
         der2=(der-lag(der))/(time-lag(time))) |> 
  summarise(inv_30=mean(der,na.rm=TRUE)) 
inv_30<-bind_rows(inv_30,
                  data.frame(species="Quercus_suber",
                             inv_30=quersube$inv_30[[1]]))
save(inv_30,file="inv_30.RData")
load("inv_30.RData")
load("sim_tot_inv_30.RData")
# save(sim_tot,file="sim_tot_inv_30.RData")
# load("inv_18.RData")

## get demo traits
tar_load(mean_demo)

traits_nfi<-read.csv("data/traitsNFI.csv",sep = " ") 
traits_nfi[traits_nfi$species=="Betula pendula","species"]="Betula"

traits_MASTIF<-read.csv("data/finalFile_maturation.csv") |> 
  filter(speciesN %in% c("Betula_pendula",gsub(" ","_",mean_demo$species))) |> 
  select(speciesN,TSM,sla,wood_density,seed_size,SSP) |> 
  rename(species=speciesN) |> 
  mutate(species=if_else(species=="Betula_pendula","Betula",species))
  
traits<-read.csv("data/traits_complete.csv") |> 
  # left_join(inv_18_2 |> mutate(species=gsub("_"," ",species))) |> 
  left_join(inv_30|> mutate(species=gsub("_"," ",species))) |> 
  left_join(mean_demo) |> 
  left_join(traits_rec|> mutate(species=gsub("_"," ",species)))
  # left_join(traits_MASTIF|> mutate(species=gsub("_"," ",species))) |> 
  # left_join(traits_nfi) |>
  # select(-c(pca1,pca2,taxa,sgdd,wai,sgddb,waib,wai2,sgdd2,bark.thickness_mm)) |> 
  drop_na() 
GGally::ggpairs(traits[,c(3:20)])
# c('shade','SLA','LN','LT','WD','HM','inv_18','ba_equil','inv_50','height.dbh.ratio','growth.max')
#c(shade','SLA','LN','LT','WD','HM','inv_18','inv_30','ba_equil','inv_50','TSM','sla','wood_density','seed_size','SSP','height.dbh.ratio','growth.max','RGR')
pca_traits<-prcomp(traits[,c('recruitment',
                             "HM",
                             'WD')],center=TRUE,scale=TRUE) #c(4:9,12:14)

factoextra::fviz_pca_var(pca_traits)

pca_ind <- as.data.frame(pca_traits$x)
pca_ind$species <- traits$species
pca_ind$taxa=traits$taxa

# Plot individuals and variables on the same biplot, labeling points with the species column
factoextra::fviz_pca_biplot(pca_traits, 
                geom.ind = "point",           # Show individuals as points
                col.var = "black",             # Color of variable arrows
                col.ind = pca_ind$taxa,            # Color of individual points
                label = "var",                # Label only variables
                repel = TRUE) +               # Avoid text overlap for variable labels
  geom_text(data = pca_ind, 
            aes(x = PC1, y = PC2, label = species, color= as.factor(taxa)), 
            vjust = -0.5, 
            hjust = 0.5)    

