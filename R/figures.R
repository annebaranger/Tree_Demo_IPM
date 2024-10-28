library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(viridis)
library(rworldmap)
## maps of plots for one species
tar_load(FUNDIV_data)
sp_ex="Abies alba"


worldmap <- sf::st_as_sf(getMap(resolution = "high"))

fun_other<- FUNDIV_data %>% 
  filter(species!=sp_ex) %>% 
  select(longitude,latitude,plotcode,species,pca1) %>% 
  mutate(pca1=NA) %>% 
  unique()

fun_sp<-FUNDIV_data %>% 
  filter(species==sp_ex) %>% 
  select(longitude,latitude,plotcode,species,pca1) %>% 
  mutate(pca_brk=cut(pca1,  # Create climate categories based on PCA1
                     breaks = quantile(pca1, probs = seq(0, 1, length.out = 10 + 1)),  # Use quantiles for equal-sized categories
                     include.lowest = TRUE),
         clim_id = as.factor(as.integer(factor(pca_brk)))) |> 
  unique()

ggplot()+
  geom_sf(data=worldmap,fill="grey")+
  geom_point(data=fun_other,aes(longitude,latitude),color="grey13",size=0.2)+
  geom_point(data=fun_sp,aes(longitude,latitude,color=pca1),size=0.5)+
  scale_color_gradientn(colors=viridis(n=15))+
  ylim(c(35, 70)) + xlim(c(-10, 35))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="",y="",color="PCA first \naxis")

color_palette <- colorRampPalette(c("red","orange","tan","green","blue"))
colors <- color_palette(10)
names(colors) <- levels(fun_sp$clim_id)

ggplot()+
  geom_sf(data=worldmap,fill="grey")+
  geom_point(data=fun_other,aes(longitude,latitude),color="grey13",size=0.2)+
  geom_point(data=fun_sp,aes(longitude,latitude,color=clim_id),size=0.8)+
  scale_colour_manual(
    name = "Climate \ncategories",
    values = colors
  ) +
  ylim(c(35, 70)) + xlim(c(-10, 35))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="",y="",color="Climate \ncategories")


# distribution climate
pca_breaks=seq(from=min(fun_sp$pca1,na.rm=TRUE),
               to=max(fun_sp$pca1,na.rm=TRUE),
               length.out=(11))
fun_sp %>% ggplot(aes(pca1))+
  geom_density()+
  geom_segment(data = data.frame(br=pca_breaks,col=1:11),
               aes(x=br,xend=br,y=-0.01,yend=0.01,color=col),
               linewidth=2)+
  scale_color_gradientn(colors=turbo(11, direction = -1) )+
  theme_bw( )+
  theme(panel.grid=element_blank(),legend.position = "none")+
  labs(y="Density",x="PCA first axis", title="Abies alba")


# ex simulation equilibrium
library(matreex)
species.combination=tar_read(sim_forest_list)$list.forests
species_list=tar_read(species_list.select)
species_object=tar_read(species_object_mu)
tar_load(harv_rules.ref)
sim.type="mu"
tar_load(sim_equil.id)
id_forest=90


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
sim.in %>% filter(var=="BAsp") %>% 
  # filter(species=="Abies_alba") %>% 
  ggplot(aes(time,value,color=species))+
  geom_line(size=0.81)+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  ylim(c(10,65))+
  labs(y="Basal area (m2/ha)",x="Time (year)",color="Species")


sim.in %>%
  filter(var == "n", equil) %>% 
  filter(size!=0) %>% 
  ggplot(aes(size,value,fill=species,color=species))+geom_col( )+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  facet_wrap(~species,ncol=3)+
  labs(y="Density",x="Size", title="Size distribution at equilibrium")


# species alone
list.sp=list.species[1]
forest.sp= new_forest(species = list.sp, harv_rules = harv_rules.ref)
sim.sp = sim_deter_forest(forest.sp, 
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
sim.sp %>% filter(var=="BAsp") %>% 
  # filter(species=="Abies_alba") %>% 
  ggplot(aes(time,value,color=species))+
  geom_line()+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  ylim(c(10,66))+
  labs(y="Basal area (m2/ha)",x="Time (year)")

rbind(sim.in |>
        filter(species=="Abies_alba") |> 
        mutate(sim="In competition"),
      sim.sp |> 
        mutate(sim="Species alone")) |>
  filter(var=="BAsp") |> 
  ggplot(aes(time,value,color=sim))+
  geom_line(size=1)+
  scale_color_manual(values=c("firebrick2","slateblue3"))+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  labs(y="Basal area (m2/ha)",x="Time (year)",color="")
  


# sim invasion
list.species.inv=list.species[c("Picea_abies","Pinus_sylvestris")]
forest.inv = new_forest(species = list.species.inv, harv_rules = harv_rules.ref)
sim.inv.void = sim_deter_forest(forest.inv, 
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
sim.inv.void %>% filter(var=="BAsp") %>% 
  ggplot(aes(time,value,color=species))+
  geom_line()+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  scale_color_manual(values=c("green4","cornflowerblue"))+
  labs(y="Basal area (m2/ha)",x="Time (year)", title="Equilibrium simulation")


sim.inv.void %>%
  filter(var == "n", equil) %>% 
  filter(size!=0) %>% 
  ggplot(aes(size,value,fill=species,color=species))+geom_col( )+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  scale_color_manual(values=c("green4","cornflowerblue"))+
  scale_fill_manual(values=c("green4","cornflowerblue"))+
  facet_wrap(~species,ncol=2)+
  labs(y="Density",x="Size", title="Size distribution at equilibrium")

sim.inv=readRDS("rds/Abies_alba/clim_6/sim_invasion/Abies_alba.Picea_abies.Pinus_sylvestris.rds")
sim.inv %>% filter(var=="BAsp") %>% 
  ggplot(aes(time,value,color=species))+
  geom_line(size=0.81)+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  labs(y="Basal area (m2/ha)",x="Time (year)")#, title="Equilibrium simulation"


sim.inv=readRDS("rds/Abies_alba/clim_7/sim_invasion/Abies_alba.Fagus_sylvatica.rds")
# sim disturbance
sim.dist=readRDS("rds/Abies_alba/clim_6/sim_disturbance/Abies_alba.Picea_abies.Pinus_sylvestris.rds")
sim.dist %>% filter(var=="BAsp") %>% 
  filter(species=="Abies_alba") %>% 
  ggplot(aes(time,value,color=species))+
  geom_line()+
  theme_bw( )+
  theme(panel.grid=element_blank())+
  labs(y="Basal area (m2/ha)",x="Time (year)", title="Equilibrium simulation")
