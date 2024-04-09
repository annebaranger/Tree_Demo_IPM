inv_test=invasion_metric |>
  rowwise() |> 
  mutate(ncomb=length(unlist(strsplit(species_combination,"[.]")))) |> 
  ungroup()
inv_test|> 
  filter(ncomb!=1) |> 
  ggplot(aes(wai_id,sgdd_id,color=inv_50))+
  geom_point(size=6,alpha=2)+
  facet_wrap(~ncomb)+
  scale_color_gradientn(colours = viridis(15))


inv_test |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  group_by(ncomb,wai_id,inv) |> 
  summarize(mean.inv=mean(value)) |> 
  ggplot(aes(ncomb,mean.inv,color=inv))+
  geom_point()+
  facet_wrap(~wai_id)

inv_test |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  group_by(ncomb,wai_id,sgdd_id,inv) |> 
  summarize(mean.inv=mean(value)) |>
  filter(inv=="inv_50") |> 
  ggplot(aes(wai_id,mean.inv,color=sgdd_id))+
  geom_point()+
  facet_wrap(~ncomb)+
  scale_color_gradientn(colours = viridis(15))


inv_test |> 
  pivot_longer(cols=c("inv_mean","inv_max","inv_50"),names_to = "inv") |> 
  group_by(ncomb,wai_id,sgdd_id,inv) |> 
  summarize(mean.inv=mean(value)) |> 
  ggplot(aes(ncomb,mean.inv,color=inv))+
  geom_point()+
  facet_grid(sgdd_id~wai_id)


summary(lm(inv_50~ncomb+wai+sgdd,data=inv_test))
car::Anova(lm(inv_50~ncomb+wai+sgdd,data=inv_test))
