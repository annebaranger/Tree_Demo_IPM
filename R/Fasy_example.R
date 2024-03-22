library(matreex)
library(dplyr)
library(tidyr)
library(ggplot2)
library(targets)


tar_load(FUNDIV_data)
worldmap <- sf::st_as_sf(rworldmap::getMap(resolution = "high"))  

disturb_coef.in=data.table::fread("data/disturb_coef.csv") |> 
  filter(disturbance=="storm")
disturbance.df_storm= data.frame(
  type = rep("storm", 3), intensity = rep(0.5, 3), 
  IsSurv = rep(FALSE, 3), t = c(500:502))
sp="Fagus sylvatica"#"Picea abies"#
s_p=gsub(pattern = " ", replacement = "_",  x = sp)
## delimit climate
FUNDIV_data |> 
  filter(species==sp) |> 
  mutate(wai_cat=cut(wai, breaks = 10,labels=1:10),
         sgdd_cat=cut(sgdd,breaks=10,labels=1:10)) |> 
  group_by(wai_cat,sgdd_cat) |> 
  summarise(n=n()) |> 
  ggplot(aes(wai_cat,sgdd_cat,fill=log(n)))+
  geom_tile()


FUNDIV_data |> 
  filter(species==sp) |> 
  mutate(wai_cat=cut(wai, breaks = 10,labels=1:10),
         sgdd_cat=cut(sgdd,breaks=10,labels=1:10)) |> 
  sample_n(50000) |> 
  ggplot()+
  geom_point(aes(x=longitude,y=latitude,color=wai))+
  geom_sf(data=worldmap,fill=NA)+
  xlim(c(-10,25))+ylim(c(35,65))

tar_load(species.list.ipm)
FUNDIV_data |> 
  filter(species %in% gsub("_"," ",species.list.ipm)) |> 
  pivot_longer(cols=c("wai","sgdd")) |> 
  ggplot(aes(value,color=species))+
  geom_density()+
  theme(legend.position = "none")+
  facet_wrap(~name,scales="free")

FUNDIV_data |> 
  filter(species %in% gsub("_"," ",species.list.ipm)) |> 
  pivot_longer(cols=c("wai","sgdd")) |> 
  group_by(species,name) |> 
  summarise(range_clim=quantile(value,probs = 0.97)-quantile(value,probs = 0.03)) |> 
  arrange(name,range_clim) |> 
  group_by(name) |> 
  mutate(extrema=case_when(range_clim==min(range_clim)~"min",
                           range_clim==max(range_clim)~"max")) |> 
  filter(extrema%in%c("min","max")) |> 
  select(-species) |> 
  pivot_wider(names_from = extrema,
              values_from = range_clim) |> 
  mutate(step_min=min/4,
         step_max=max/10) |> 
  mutate(step=max(step_min,step_max)) |> 
  ungroup() |> 
  select(name,step)->step_cat
FUNDIV_data |> 
  filter(species %in% gsub("_"," ",species.list.ipm)) |> 
  pivot_longer(cols=c("wai","sgdd")) |> 
  group_by(species,name) |> 
  summarise(range_clim=quantile(value,probs = 0.97)-quantile(value,probs = 0.03)) |> 
  left_join(step_cat) |> 
  mutate(n_breaks=round(range_clim/step)) |> 
  select(species,name, n_breaks) |> 
  pivot_wider(names_from = name,
              values_from = n_breaks) |> 
  rename(sgdd_breaks=sgdd,wai_breaks=wai)->step_species

FUNDIV_data |> 
  filter(species %in% gsub("_"," ",species.list.ipm)) |> 
  mutate(wai_cat=cut(wai, breaks = step_species[step_species$species==sp,"wai_breaks"][[1]]),
         sgdd_cat=cut(sgdd,breaks=step_species[step_species$species==sp,"sgdd_breaks"][[1]])) |> 
  select(species,wai_cat,sgdd_cat) |> 
  unique() |> 
  group_by(species) |> 
  summarise(n=n()) |> View()

FUNDIV_data |> 
  filter(species %in% gsub("_"," ",species.list.ipm)) |> 
  group_by(species) |> 
  filter(species%in% c("Abies alba","Quercus ilex","Fagus sylvatica")) |> 
  mutate(wai_cat=cut(wai, breaks = step_species[step_species$species==sp,"wai_breaks"][[1]]),
         sgdd_cat=cut(sgdd,breaks=step_species[step_species$species==sp,"sgdd_breaks"][[1]]),
         clim_cat=factor(paste0(wai_cat,sgdd_cat))) |> 
  # sample_n(50000) |> 
  ggplot()+
  geom_point(aes(x=longitude,y=latitude,color=clim_cat))+
  geom_sf(data=worldmap,fill=NA)+
  theme(legend.position = "none")+
  xlim(c(-10,25))+ylim(c(35,65))+
  facet_wrap(~species)


### Climatic condition and species association ###
n_cat=3
FUNDIV_plotcat=FUNDIV_data |>
  filter(species%in%gsub("_"," ",species.list.ipm)) |>
  filter(species=="Abies alba") |>
  dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
                species,BAtot) |> 
  group_by(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
           species) |> 
  summarize(BA=max(BAtot)) |> 
  group_by(species) |> 
  mutate(wai_cat=cut(wai, breaks = n_cat,labels=1:n_cat),
         sgdd_cat=cut(sgdd,breaks=n_cat,labels=1:n_cat),
         BA=max(BA,na.rm = TRUE)) |> 
  ungroup()

condi.init<-FUNDIV_plotcat |> 
  group_by(species,BA,wai_cat,sgdd_cat) |> 
  summarize(wai=mean(wai),
            sgdd=mean(sgdd)) |> 
  mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
         waib = 1/(1 + wai), PC1=0, PC2=0, N = 2, SDM = 0) 


nsp_per_richness=10
prop_threshold=0.8
FUNDIV_data |> 
  left_join(FUNDIV_plotcat |> select(plotcode,wai_cat,sgdd_cat)) |> 
  select(treecode,plotcode,species,wai_cat,sgdd_cat) |> 
  filter(!is.na(wai_cat)) %>%
  filter(species%in%gsub("_"," ",species.list.ipm)) |>
  # Group by the relevant categories
  group_by(wai_cat, sgdd_cat, plotcode) %>%
  # Summarize the species occurrences in each group
  summarise(species_combination = paste(sort(unique(gsub(" ","_",species))), collapse = "."),
            n_species = n_distinct(species), .groups = 'drop') %>%
  # Now count each unique combination's frequency within each wai_cat and sgdd_cat group
  count(wai_cat, sgdd_cat, species_combination,n_species) %>%
  filter(n_species<nsp_per_richness) |> 
  group_by(wai_cat,sgdd_cat) |> 
  mutate(prop=n/sum(n)) |> 
  # Optionally, arrange the results for better readability
  arrange(wai_cat, sgdd_cat, desc(n)) |> 
  mutate(prop_cum=cumsum(prop)) |> 
  filter(prop_cum<prop_threshold) |> 
  View()
### Repeated simulation for different climate ###
#################################################
sym_clim=data.frame(species=as.character(),
                    time=as.numeric(),
                    BAsp=as.numeric(),
                    ba_init=as.numeric(),
                    sgdd_cat=factor(),
                    wai_cat=factor())
for (i in 1:dim(condi.init)[1]){
  climate=condi.init[i,c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")]
  fit_sp=eval(parse(text=paste0("fit_",s_p)))
  sp_ipm <- make_IPM(
    species = s_p, 
    climate = climate, 
    fit = fit_sp,
    clim_lab = paste0("climate_sgdd",condi.init[1,"sgdd_cat"],"_wai",condi.init[1,"wai_cat"]), 
    mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Fagus_sylvatica) * 1.1),
    BA = 0:condi.init[[i,"BA"]], # Default values are 0:200, smaller values speed up this vignette.
    verbose = TRUE
  )
  
  for(ba_init in seq(10,50,by=10)){
    Fasy_sp <- species(IPM = Fasy_ipm, init_pop = def_initBA(ba_init))
    Fasy_for <- forest(species = list(Fasy = Fasy_sp))  
    Fasy_sim <- sim_deter_forest(
      Fasy_for, 
      tlim = 1000, 
      equil_time = 3000, equil_dist = 50, equil_diff = 1,
      SurfEch = 0.03,
      verbose = TRUE
    )
    sym_clim<-sym_clim |> 
      bind_rows(Fasy_sim |> 
                  dplyr::filter(var == "BAsp", ! equil) |> 
                  dplyr::select(species,time,value) |> 
                  rename(BAsp=value) |> 
                  mutate(ba_init=ba_init,
                         wai_cat=condi.init[[i,"wai_cat"]],
                         sgdd_cat=condi.init[[i,"sgdd_cat"]]))
  
  }
}
sym_clim |>
  ggplot(aes(x = time, y = BAsp, color=as.factor(ba_init))) +
  geom_line(linewidth = .4) + ylab("BA")+
   facet_grid(wai_cat~sgdd_cat)

Fasy_sim |> 
  dplyr::filter(var == "n", time == 150,size!=0) %>%
  ggplot(aes(size,value))+
  geom_col()

distrib_t150 <- Fasy_sim |> 
  dplyr::filter(var == "n", time == 150) |> 
  pull(value)

Fasy_sim_init150 <- sim_deter_forest(
  forest(species = list(Fasy =  species(IPM = Fasy_ipm, 
                                        init_pop = def_init_k(distrib_t150)))) , 
  tlim = 1000, 
  equil_time = 3000, equil_dist = 50, equil_diff = 1,
  SurfEch = 0.03,
  verbose = TRUE
)
Fasy_sim_init150 |> 
  filter(var=="BAsp") |> 
  ggplot(aes(time,value))+
  geom_line()


Piab_ipm <- make_IPM(
  species = "Picea abies", 
  climate = climate, 
  fit = fit_Picea_abies,
  clim_lab = paste0("climate_sgdd",condi.init[1,"sgdd_cat"],"_wai",condi.init[1,"wai_cat"]), 
  mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Fagus_sylvatica) * 1.1),
  BA = 0:condi.init[[i,"BA"]], # Default values are 0:200, smaller values speed up this vignette.
  verbose = TRUE
)

### make simulation from equilibrium ###
########################################
climate=condi.init[1,c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")]
fit_sp=eval(parse(text=paste0("fit_",s_p)))
sp_ipm <- make_IPM(
  species = s_p, 
  climate = climate, 
  fit = fit_sp,
  clim_lab = paste0("climate_sgdd",condi.init[1,"sgdd_cat"],"_wai",condi.init[1,"wai_cat"]), 
  mesh = c(m = 700, L = 90, U = get_maxdbh(fit_sp) * 1.1),
  BA = 0:condi.init[[1,"BA"]], # Default values are 0:200, smaller values speed up this vignette.
  verbose = TRUE
)

sp_sp <- species(IPM = sp_ipm, init_pop = def_initBA(20))
species = list(sp_sp)
names(species)=s_p
sp_for <- forest(species = species)  

sim.eq <- sim_deter_forest(sp_for,
                           tlim = 4000, 
                           equil_time = 50000, 
                           equil_dist = 2000, 
                           equil_diff = 0.5, 
                           harvest = "default",
                           SurfEch = 0.03, 
                           verbose = TRUE)
sim.eq |> 
  filter(var=="BAsp") |> 
  ggplot(aes(time,value))+
  geom_line()

# extract equil dist
reached_equil = !is.na(sum((sim.eq %>%
                              filter(var == "BAsp") %>%
                              filter(time == max(.$time) - 1))$value))
if(reached_equil){
  # Extract the equilibrium for species i
  equil.dist = sim.eq %>%
    filter(var == "n", equil) %>% 
    pull(value)
  
  equil.dist.20<-sim.eq %>%
    filter(var == "n", equil) %>% 
    mutate(value=case_when(size>200~0,
                           TRUE~value)) |>
    mutate(BAtot=sum((size/2000)^2*pi*value)) |> 
    pull(value)
  
  ba.eq=sim.eq[sim.eq$time==max(sim.eq$time)-1 & sim.eq$var=="BAsp","value"][[1]]
  
  # Initiate the population at equilibrium
  sp_sp.eq=sp_sp
  sp_sp.eq$init_pop <- def_init_k(equil.dist*1/ba.eq)
  sp_sp.eq$init_pop <- def_init_k(equil.dist.20)
  
  species.eq = list(sp_sp.eq)
  names(species.eq)=s_p
  sp_for.eq <- forest(species = species.eq)  
  sp_for.eq <- forest(species = list(species.in))  
  
  # Run simulation till equilibrium
  sim.in.20 = sim_deter_forest(
    sp_for.eq, tlim = 4000, equil_time = 4000,
    SurfEch = 0.03, verbose = TRUE)
  sim.in.20 |> 
    filter(var=="BAsp") |> 
    ggplot(aes(time,value))+
    geom_line()
} 
bind_rows(sim.eq |> mutate(sim="eq"),
          sim.in |> mutate(sim="low")) |>
  filter(var=="BAsp") |>
  mutate(der=(value-lag(value))/(time-lag(time))) |> 
  ggplot()+
  # geom_line(aes(time,value,color=sim))+
  geom_line(aes(time,der,color=sim),linetype="dashed")+
  NULL

sp_delay<-as.numeric(sp_sp$IPM$info[["delay"]])
data <- bind_rows(sim.in.20 |> mutate(sim="low20"),
                  sim.in |> mutate(sim="low")) |> 
  filter(var=="BAsp") |> 
  group_by(sim) |>
  mutate(der=(value-lag(value))/(time-lag(time)),
         der2=(der-lag(der))/(time-lag(time))) |> 
  ungroup() 
data |> ggplot(aes(time,der,color=sim))+
  geom_line()
# Find a transformation ratio; this is a mock-up as the real relationship needs to be meaningful
# In this case, we assume they are just directly relatable which might not be true in your real data context
# You need to adjust this transformation based on actual interpretive relationships between 'value' and 'der'
max_der <- max(data$der, na.rm = TRUE)
max_value <- max(data$value, na.rm = TRUE)
ratio <- max_value / max_der

# Creating the plot
ggplot(data) +
  geom_line(aes(x = time, y = value, color = sim)) +
  geom_line(aes(x = time, y = der * ratio, color = sim), linetype = "dashed") +
  scale_y_continuous(
    "Value",
    sec.axis = sec_axis(~ . / ratio, name = "Derivative")
  ) +
  labs(x = "Time") +
  geom_vline(xintercept = sp_delay)


first.seq<-data |>  
  filter(time>(sp_delay+5)) |> 
  filter(time<1500) |> 
  group_by(sim) |> 
  mutate(sign_eq=(sign(der)==sign(lag(der)))) |> 
  filter(sign_eq==FALSE) |> slice(1) |> 
  select(sim,time) |> rename(max=time)

data |> left_join(first.seq) |> 
  group_by(sim) |> 
  filter(time<(max-15)) |> 
  summarise(inv_mean=mean(der,na.rm=TRUE),
            inv_max=max(der,na.rm=TRUE))


### make perturbation from equilibrium ###
########################################
equil.dist = sim.eq %>%
  filter(var == "n", equil) %>% 
  pull(value)

ba.eq=sim.eq[sim.eq$time==max(sim.eq$time)-1 & sim.eq$var=="BAsp","value"][[1]]

# Initiate the population at equilibrium
sp_sp.dist=sp_sp
sp_sp.dist$init_pop <- def_init_k(equil.dist)

# Update disturbance function
sp_sp.dist$disturb_fun <- disturb_fun

# Add disturbance coefficients
sp_sp.dist$disturb_coef <- filter(disturb_coef.in,
                                  species == s_p)



species.dist = list(sp_sp.dist)
names(species.dist)=s_p
sp_for.dist <- forest(species = species.dist)  

# Run simulation till equilibrium
sim.dist = sim_deter_forest(
  sp_sp.dist, tlim = 4000, equil_time = 4000,
  disturbance = disturbance.df_storm, 
  SurfEch = 0.03, verbose = TRUE)
sim.dist |> 
  filter(var=="BAsp") |> 
  ggplot(aes(time,value))+
  geom_line()
