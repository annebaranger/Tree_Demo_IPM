library(matreex)
library(dplyr)
library(tidyr)
library(ggplot2)
library(targets)


tar_load(FUNDIV_data)
tar_load(species.list.ipm)
tar_load(harv_rules.ref)
tar_load(disturbance.df_storm)

sp="Fagus sylvatica"
s_p="Fagus_sylvatica"
# create climate category
FUNDIV_plot=FUNDIV_data |> 
  filter(species %in% gsub("_"," ",species.list.ipm)) |> 
  dplyr::select(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
                species,BAtot) |> 
  group_by(plotcode, longitude, latitude, sgdd, wai, pca1, pca2,
           species) |> 
  # maximum Ba per plot
  summarize(BA=max(BAtot)) |> 
  distinct() |> 
  ungroup()

FUNDIV_plot |> 
  pivot_longer(cols=c("wai","sgdd")) |> 
  group_by(species,name) |> 
  summarise(range_clim=quantile(value,probs = 0.95)-quantile(value,probs = 0.05)) |> 
  arrange(name,range_clim) |> 
  group_by(name) |> 
  mutate(extrema=case_when(range_clim==min(range_clim)~"min",
                           range_clim==max(range_clim)~"max")) |> 
  filter(extrema%in%c("min","max")) |> 
  select(-species) |> 
  pivot_wider(names_from = extrema,
              values_from = range_clim) |> 
  mutate(step_min=min/6,
         step_max=max/15) |> 
  mutate(step=max(step_min,step_max)) |> 
  ungroup() |> 
  select(name,step)->step_cat_1

FUNDIV_plot |> 
  pivot_longer(cols=c("wai","sgdd")) |> 
  group_by(species,name) |> 
  summarise(range_clim=quantile(value,probs = 0.95)-quantile(value,probs = 0.05)) |> 
  left_join(step_cat) |> 
  mutate(n_breaks=round(range_clim/step)) |> 
  select(species,name, n_breaks) |> 
  pivot_wider(names_from = name,
              values_from = n_breaks) |> 
  rename(sgdd_breaks=sgdd,wai_breaks=wai)->step_species

FUNDIV_plotcat=FUNDIV_plot |> 
  left_join(step_species) |> 
  group_by(species) |> 
  mutate(wai_cat=cut(wai, 
                     breaks = seq(from=min(wai,na.rm=TRUE),
                                  to=max(wai,na.rm=TRUE),
                                  length.out=(wai_breaks+1)),
                     include.lowest = TRUE),
         sgdd_cat=cut(sgdd,
                      breaks=seq(min(sgdd,na.rm=TRUE)
                                 ,max(sgdd,na.rm=TRUE),
                                 length.out=(sgdd_breaks+1)),
                      include.lowest = TRUE)) |> 
  tidyr::separate_wider_delim(cols="wai_cat",names=c("wai_low","wai_up"),delim=",",cols_remove = FALSE) |> 
  tidyr::separate_wider_delim(cols="sgdd_cat",names=c("sgdd_low","sgdd_up"),delim=",",cols_remove = FALSE) |> 
  mutate(across(matches(c("low")),
                ~as.numeric(stringr::str_sub(.,2,-1))),
         across(matches(c("up")),
                ~as.numeric(stringr::str_sub(.,1,-2)))
         ) |> 
  ungroup()  |> 
  group_by(species) %>%
  mutate(
    wai_id = as.integer(factor(wai_cat)),
    sgdd_id = as.integer(factor(sgdd_cat))
  ) %>%
  ungroup()


condi.init<- FUNDIV_plotcat |> 
  group_by(species,wai_id,sgdd_id,wai_low,wai_up,sgdd_low,sgdd_up) |> 
  summarize(n_plot=n(),
            wai=mean(wai,na.rm=TRUE),
            sgdd=mean(sgdd,na.rm=TRUE)) |>
  # filter out category with too few plots
  filter(n_plot>10) |> 
  # create IPM vars
  mutate(sgdd2 = sgdd^2, wai2 = wai^2, sgddb = 1/sgdd, 
         waib = 1/(1 + wai), PC1=0, PC2=0, N = 2, SDM = 0) |> 
  ungroup() |> 
  group_by(species) |> 
  mutate(ID.spclim = row_number(),
         ID.species=cur_group_id()) |> 
  ungroup()

condi.init |> 
  ggplot(aes(wai_id,sgdd_id,fill=log(n_plot)))+
  geom_tile()+
  facet_wrap(~species,scales="free")

condi.init.fasy<-condi.init |> 
  filter(species==sp)

## load param demo
fit.list<-list()
fit.list[[s_p]] <- get(paste0("fit_", s_p))


## make_mu
mu_Fasy <- make_mu_gr(
  species = s_p, fit = fit.list[[s_p]],
  mesh = c(m = 700, L = 90, U = get_maxdbh(fit.list[[s_p]]) * 1.1),
  verbose = TRUE, stepMu = 0.001)

time <- 1500
init_pop_fasy=def_initBA(20)
Fasy <- species(IPM = mu_Fasy, init_pop = init_pop_fasy,
                harvest_fun = def_harv)
forest <- forest(species = list(mu_Fasy = Fasy))

for(clim in 1:dim(condi.init.fasy)[1]){
  set.seed(42)
  memor_mu <- sim_deter_forest(forest,
                               tlim = time,
                               climate = condi.init[clim,c("sgdd", "wai", "sgddb", "waib",
                                                        "wai2", "sgdd2","PC1", "PC2", "N",
                                                        "SDM")],
                               equil_dist = time, equil_time = time,
                               verbose = TRUE, correction = "cut")
  
  }


memor_mu %>%
  filter(var %in% c("BAsp", "H", "N"), ! equil, value != 0) %>%
  ggplot(aes(x = time, y = value)) +
  facet_wrap(~ var, scales = "free_y") +
  geom_line(size = .4) + geom_point(size = .4) +
  NULL


ipm_fasy <- make_IPM(
  s_p,
  climate = condi.init[1,c("sgdd", "wai", "sgddb", "waib",
                           "wai2", "sgdd2","PC1", "PC2", "N",
                           "SDM")], 
  "clim1",
  fit = fit.list[[s_p]],
  mesh = c(m = 700, L = 90, U = get_maxdbh(fit.list[[s_p]]) * 1.1),
  BA = 0:100, verbose = TRUE
)

Fasy_ipm <- species(IPM = ipm_fasy, init_pop = init_pop_fasy,
                           harvest_fun = def_harv)
forest_ipm <- forest(species = list(ipm_Fasy = Fasy_ipm))
set.seed(42)
memor_ipm <- sim_deter_forest(forest_ipm, tlim = time,
                                     equil_dist = time, equil_time = time,
                                     verbose = TRUE, correction = "cut")

e_memor <- dplyr::bind_rows(ipm = memor_ipm, mu = memor_mu, .id = "meth")
e_memor %>%
  filter(var %in% c("BAsp", "N"), ! equil, value != 0) %>%
  ggplot(aes(x = time, y = value, color = meth)) +
  facet_wrap(~ var, scales = "free_y") +
  geom_line(size = .4, linetype = "dotted") + geom_point(size = .4) +
  NULL
e_memor %>%
  filter(var %in% c("BAsp"), ! equil) %>%
  group_by(time) %>% summarise(value = diff(value)) %>%
  ggplot(aes(x = time, y = value)) +
  ylab("BA difference") +
  geom_line(size = .4, linetype = "dotted") + geom_point(size = .4) +
  NULL

e_memor %>%
  filter(var == "n", time ==1) %>%
  ggplot(aes(x = size, y = value, color = meth)) +
  geom_col()
  NULL
