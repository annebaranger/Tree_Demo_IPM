library(targets)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(factoextra)
tar_load(fit.list.allspecies)
tar_load(species.list.ipm)

fit.list.allspecies<-fit.list.allspecies[species.list.ipm]

list_pars_sp=data.frame(species=names(fit.list.allspecies))
for(sp in 1:length(fit.list.allspecies)){
  sp_fit=fit.list.allspecies[[sp]]
  sp=names(fit.list.allspecies)[sp]
  for (vr in c("sv","gr","rec")){
    fit_vr=fit.list.allspecies[[sp]][[vr]][["params_m"]]
    for(vvr in 1:length(fit_vr)){
      pars=paste0(vr,"_",names(fit_vr)[vvr])
      list_pars_sp[list_pars_sp$species==sp,pars]=fit_vr[vvr][[1]]
    }
  } 
}
na_count <-data.frame(pars=colnames(list_pars_sp),
                      na_count=sapply(list_pars_sp, function(y) sum(length(which(is.na(y)))))) |> 
  filter(na_count<5)
list_pars_sp<-list_pars_sp |>
  # select(na_count$pars) |> 
  select_if(~ !any(is.na(.))) |>
  # mutate(across(everything(),~replace_na(.,0)))
  # select(species,matches("intercept"),matches("BATOT"),matches("size"),matches("logsize")) |> 
  select(!matches("log"))
pca<-prcomp(list_pars_sp[2:dim(list_pars_sp)[2]],center = TRUE, scale = TRUE)
fviz_pca_ind(pca)
fviz_pca_ind(pca,axes=c(2,3))
fviz_pca_var(pca,col.var="contrib")
fviz_pca_var(pca,axes=c(2,3),col.var="contrib")

tar_load(species.list.ipm)
traits<-tar_read(traits) %>% filter(species%in% species.list.ipm) %>% 
  left_join(list_pars_sp[,c("species","rec_intercept")])
  
pca<-prcomp(traits[,2:5],center=TRUE,scale=TRUE)
library(factoextra)
fviz_pca_var(pca)
