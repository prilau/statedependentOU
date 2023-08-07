#Looking at distributions of the dataset

library(ape)
library(tidyverse)

Traits <- read_csv("bird_traits_joined.csv")

flapper0 <- Traits %>% filter(Flapper == 0)
flapper1 <- Traits %>% filter(Flapper == 1)

hist(flapper1$meanHUratio, ylim = c(0, 100), xlim = c(0, 2))
hist(flapper0$meanHUratio, add = TRUE, col = "red")

var(Traits$meanHUratio)

migrant0 <- Traits %>% filter(Migrant == 0)
migrant1 <- Traits %>% filter(Migrant == 1)

hist(migrant1$meanHUratio, ylim = c(0, 40), xlim = c(0.5, 2))
hist(migrant0$meanHUratio, add = TRUE, col = "red")
