library(targets)
library(dplyr)
library(tidyr)
library(stringr)
tar_load(FUNDIV_data)
tar_load(species.combination)
tar_load(climate.cat)
tar_load(species.list.ipm)
## observed most frequent species.combination

out<-setNames(data.frame(matrix( nrow = 0,ncol=5)),
              nm=c("clim_id","partner","coocc","full","species"))
for (sp in species.list.ipm){
  print(sp)
  climate.cat.sp<-climate.cat$species.cat %>% 
    filter(species==gsub("_"," ",sp))
  clim_breaks=c(-Inf,
               sort(climate.cat.sp$clim_up)[1:9],
               Inf)
  zz<-FUNDIV_data %>% 
    select(plotcode,species,pca1) %>% 
    unique() %>% 
    filter(!is.na(species)) %>% 
    mutate(clim_cat=cut(pca1,
                        breaks=clim_breaks),
           clim_id=as.numeric(clim_cat)) %>% 
    group_by(clim_id,clim_cat,species) %>% 
    count() %>% 
    ungroup() %>% 
    filter(species!=gsub("_"," ",sp)) %>% 
    group_by(clim_id) %>% 
    mutate(prop=100*n/sum(n)) %>% 
    arrange(clim_id,desc(prop)) %>% 
    mutate(cumprop=cumsum(prop)) %>% 
    filter(cumprop<80) %>% 
    select(clim_id,species) %>%
    mutate(full=1,
           species=gsub(" ","_",species)) %>% 
    ungroup()
  
  species.combination.sp=species.combination %>% filter(species==gsub("_"," ",sp)) %>% 
    select(clim_id,species_combination) %>% 
    rename(species=species_combination) %>% 
    rowwise() %>% 
    mutate(species=strsplit(species,"\\.")) |> 
    unnest(cols = species) %>% 
    unique() %>% 
    mutate(coocc=1) %>% 
    filter(species!=sp)
  
  out<- rbind(out,
             full_join(species.combination.sp,zz,by=c("clim_id","species")) %>%
               rename(partner=species) %>%
               mutate(species=sp))
  
}


out %>% 
  filter(partner%in%species.list.ipm) %>% 
  filter(full==1&is.na(coocc)) %>% 
  group_by(species,clim_id) %>%
  summarise(n=n()) %>% 
  ggplot(aes(clim_id,n))+
  geom_col()+
  facet_wrap(~species)
  