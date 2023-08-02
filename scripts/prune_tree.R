library(ape)
tree <- read.nexus("bird_tree.nex")

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
Traits <- Traits %>% 
  select(1, 4)

contTrait <- list()

for (i in 1:length(Traits$Species)){
  contTrait[[Traits$Species[i]]] <- list(Traits$meanHUratio[i])
}
write.nexus.data(contTrait, "./birdForelimbRatio.nex", format = "continuous")

newTree <- keep.tip(tree, Traits$Species)
library(RRphylo)
newNewTree <- fix.poly(newTree,type="resolve",node=NULL,tol=1e-10,random=TRUE)
is.binary(newNewTree)
write.tree(newNewTree, "bird_tree_pruned_resolvedPoly.tre")


flapper0 <- Traits %>% filter(Flapper == 0)
flapper1 <- Traits %>% filter(Flapper == 1)

hist(flapper1$meanHUratio, ylim = c(0, 100), xlim = c(0, 2))
hist(flapper0$meanHUratio, add = TRUE, col = "red")

var(Traits$meanHUratio)


