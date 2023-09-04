# Create bird datasets
## Discrete traits: flapper, migrant; continuous trait: humerus-ulna ratio

library(ape)
library(tidyverse)

discreteTrait <- read_csv("./bird_brain.csv")
contTrait <- read_csv("./Jess_skel_methods.csv")

discreteTrait <- discreteTrait %>% dplyr::select(2, 4, 5)

contTrait <- contTrait %>%
  filter(sex == "M") %>% 
  mutate(limbRatio = humerus/femur) %>% 
  drop_na() %>% 
  dplyr::select(2, 16) %>% 
  group_by(genus_species) %>% 
  reframe(meanLimbRatio = mean(limbRatio)) %>% 
  rename(Species = 'genus_species')

Traits <- discreteTrait %>% inner_join(contTrait) %>% 
  drop_na()


# Making continuous trait nexus file with R coz RevBayes couldn't read my by-hand version
limbRatio <- list()
for (i in 1:length(Traits$Species)){
  limbRatio[[Traits$Species[i]]] <- list(Traits$meanLimbRatio[i])
}
write.nexus.data(limbRatio, "./birdLimbRatio.nex", format = "continuous")

mean(Traits$meanLimbRatio)


#Looking at distributions of the dataset
flapper0 <- Traits %>% filter(Flapper == 0)
flapper1 <- Traits %>% filter(Flapper == 1)

hist(flapper1$meanLimbRatio, ylim = c(0, 100), xlim = c(0, 6))
hist(flapper0$meanLimbRatio, add = TRUE, col = "red")

var(Traits$meanLimbRatio)

migrant0 <- Traits %>% filter(Migrant == 0)
migrant1 <- Traits %>% filter(Migrant == 1)

hist(migrant0$meanLimbRatio, ylim = c(0, 100), xlim = c(0, 6))
hist(migrant1$meanLimbRatio, add = TRUE, col = "red")


