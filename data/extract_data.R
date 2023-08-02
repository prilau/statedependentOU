discreteTrait <- read_csv("./bird_brain.csv")
contTrait <- read_csv("./Jess_skel_methods.csv")

discreteTrait <- discreteTrait %>% dplyr::select(2, 4, 5)
contTrait <- contTrait %>%
  filter(sex == "M") %>% 
  mutate(forelimbRatio = humerus/ulna) %>% 
  drop_na() %>% 
  dplyr::select(2, 16) %>% 
  group_by(genus_species) %>% 
  reframe(meanHUratio = mean(forelimbRatio)) %>% 
  rename(Species = 'genus_species')



Traits <- discreteTrait %>% inner_join(contTrait) %>% 
  drop_na()
write_csv(Traits, file = "bird_traits_joined.csv")  











