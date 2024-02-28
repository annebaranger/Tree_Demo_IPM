tar_load(FUNDIV_data)

## delimit climate
FUNDIV_data |> 
  filter(species=="Fagus sylvatica") |> 
  mutate(wai_cat=cut(wai, breaks = 10,labels=1:10),
         sgdd_cat=cut(sgdd,breaks=10,labels=1:10)) |> 
  group_by(wai_cat,sgdd_cat) |> 
  summarise(n=n()) |> 
  ggplot(aes(wai_cat,sgdd_cat,fill=log(n)))+
  geom_tile()

condi.init=FUNDIV_data |>
  dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
                species,BAtot) |> 
  group_by(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
           species) |> 
  summarize(BA=max(BAtot)) |> 
  ungroup() |> 
  unique() |> 
  filter(species=="Fagus sylvatica") |> 
  mutate(wai_cat=cut(wai, breaks = 10,labels=1:10),
         sgdd_cat=cut(sgdd,breaks=10,labels=1:10)) |> 
  group_by(wai_cat,sgdd_cat) |> 
  summarize(wai=mean(wai),
            sgdd=mean(sgdd),
            BA=max(BA)) |> 
  mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
         waib = 1/(1 + wai), PC1=0, PC2=0, N = 2, SDM = 0)


climate=condi.init[1,c("sgdd", "wai", "sgddb", "waib", "wai2", "sgdd2", "PC1", "PC2", "N", "SDM")]
Fasy_ipm_1 <- make_IPM(
  species = "Fagus_sylvatica", 
  climate = climate, 
  fit = fit_Fagus_sylvatica,
  clim_lab = paste0("climate_sgdd",condi.init[1,"sgdd_cat"],"_wai",condi.init[1,"wai_cat"]), 
  mesh = c(m = 700, L = 90, U = get_maxdbh(fit_Fagus_sylvatica) * 1.1),
  BA = 0:condi.init[[1,"BA"]], # Default values are 0:200, smaller values speed up this vignette.
  verbose = TRUE
)


Fasy_sp <- species(IPM = Fasy_ipm_1, init_pop = def_initBA(10))
Fasy_for <- forest(species = list(Fasy = Fasy_sp))

set.seed(42) # The seed is here for initial population random functions.
Fasy_sim <- sim_deter_forest(
  Fasy_for, 
  tlim = 1000, 
  equil_time = 3000, equil_dist = 50, equil_diff = 1,
  SurfEch = 0.03,
  verbose = TRUE
)


Fasy_sim  %>%
  dplyr::filter(var == "BAsp", ! equil) %>%
  ggplot(aes(x = time, y = value)) +
  geom_line(linewidth = .4) + ylab("BA")
