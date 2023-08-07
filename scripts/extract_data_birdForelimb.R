# Create bird datasets
## Discrete traits: flapper, migrant; continuous trait: humerus-ulna ratio

library(ape)
library(tidyverse)

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


write_csv(Traits, "bird_traits_joined.csv")

#I made discrete traits nexus file by hand lol

#Making joined_discrete trait nexus file
Traits_comb <- Traits %>% 
  mutate(joined = Flapper + 2* Migrant) %>% 
  select(1, 5)

joinedDiscTrait <- list()

for (i in 1:length(Traits_comb$Species)){
  joinedDiscTrait[[Traits_comb$Species[i]]] <- list(Traits_comb$joined[i])
}
write.nexus.data(joinedDiscTrait, "./birdDiscrete_combined.nex", format = "standard")

# Making continuous trait nexus file with R coz RevBayes couldn't read my by-hand version
contTrait <- list()
for (i in 1:length(Traits$Species)){
  contTrait[[Traits$Species[i]]] <- list(Traits$meanHUratio[i])
}
write.nexus.data(contTrait, "./birdForelimbRatio.nex", format = "continuous")



