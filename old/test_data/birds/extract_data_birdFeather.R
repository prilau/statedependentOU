# Create bird datasets
## Discrete traits: flapper, migrant; continuous trait: humerus-ulna ratio

library(ape)
library(tidyverse)

discreteTrait <- read_csv("./bird_feather_3.csv")
contTrait <- read_csv("./bird_feather_1.csv")

discreteTrait <- discreteTrait %>% dplyr::select(1, 3, 8) %>% 
  rename(Flight = `Flight type`) %>% 
  mutate(mixFlight = Flight)

for (i in 1:length(discreteTrait$Flight)){
  if (discreteTrait$Flight[i] == "1"){
    discreteTrait$Flight[i] <- 1
    discreteTrait$mixFlight[i] <- 1
  } else if (discreteTrait$Flight[i] == "3") {
    discreteTrait$Flight[i] <- 0
    discreteTrait$mixFlight[i] <- 0
  }
  else {
    if (discreteTrait$Flight[i] == "2a"){
      discreteTrait$mixFlight[i] <- 2
    }
    else if (discreteTrait$Flight[i] == "2b"){
      discreteTrait$mixFlight[i] <- 3
    }
    else if (discreteTrait$Flight[i] == "2c"){
      discreteTrait$mixFlight[i] <- 4
    }
    discreteTrait$Flight[i] <- 2
  }
}

for (i in 1:length(discreteTrait$Habitat)){
  if (discreteTrait$Habitat[i] == "T"){
    discreteTrait$Habitat[i] <- 0
  } else if (discreteTrait$Habitat[i] == "R") {
    discreteTrait$Habitat[i] <- 1
  }
  else {
    discreteTrait$Habitat[i] <- 2
  }
}

for (i in 1:length(discreteTrait$Species)){
  m <- strsplit(discreteTrait$Species[i], split = " ")
  discreteTrait$Species[i] <- paste0(m[[1]][1], "_", m[[1]][2])
}



contTrait_p1 <- contTrait %>%
  filter(p1_p8 == "p1") %>% 
  drop_na() %>% 
  dplyr::select(1:4) %>% 
  group_by(species) %>% 
  rename(Species = 'species',
         RachisW = `Rachis width (mm)`,
         BarbD = `Barb bensity (number/cm)`,
         BarbuleD = `Barbule density (number/mm)`)



#Making discrete trait nexus file
DiscTrait <- list()

for (i in 1:length(discreteTrait$Species)){
  DiscTrait[[discreteTrait$Species[i]]]$Habitat = discreteTrait$Habitat[i]
  DiscTrait[[discreteTrait$Species[i]]]$Flight = discreteTrait$Flight[i]
  DiscTrait[[discreteTrait$Species[i]]]$mixFlight = discreteTrait$mixFlight[i]

}

write.nexus.data(DiscTrait, "./birdFeatherDiscrete.nex", format = "standard")

# Making continuous trait nexus file
contTraitP1 <- list()
for (i in 1:length(contTrait_p1$Species)){
  contTraitP1[[contTrait_p1$Species[i]]]$RachisW = contTrait_p1$RachisW[i]
  contTraitP1[[contTrait_p1$Species[i]]]$BarbD = contTrait_p1$BarbD[i]
  contTraitP1[[contTrait_p1$Species[i]]]$BarbuleD = contTrait_p1$BarbuleD[i]
}

write.nexus.data(contTraitP1, "./birdFeatherCont.nex", format = "continuous")




Traits <- discreteTrait %>% inner_join(contTrait_p1) %>% 
  drop_na()

habitat0 <- Traits %>% filter(Habitat == "T")
habitat1 <- Traits %>% filter(Habitat == "R")
habitat2 <- Traits %>% filter(Habitat == "A")

hist(habitat0$`Rachis width (mm)`, ylim = c(0, 50), xlim = c(0, 8))
hist(habitat2$`Rachis width (mm)`, add = TRUE, col = "red")
hist(habitat1$`Rachis width (mm)`, add = TRUE, col = "blue")

hist(habitat0$`Barb bensity (number/cm)`, ylim = c(0, 50), xlim = c(10, 40))
hist(habitat2$`Barb bensity (number/cm)`, add = TRUE, col = "red")
hist(habitat1$`Barb bensity (number/cm)`, add = TRUE, col = "blue")

hist(habitat0$`Barbule density (number/mm)`, ylim = c(0, 50), xlim = c(10, 60))
hist(habitat2$`Barbule density (number/mm)`, add = TRUE, col = "red")
hist(habitat1$`Barbule density (number/mm)`, add = TRUE, col = "blue")
