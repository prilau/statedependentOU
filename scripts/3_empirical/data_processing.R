library(ape)
library(phytools)
library(tidyverse)

# since I obtained the phylogeny and the trait data from different sources
# taxa don't always match
# Let's retain the taxa that are present in both the phylogeny and the trait database.

# read in tree
tree <- read.nexus("data/3_empirical/raw/4705sp_mammal-time_Álvarez-Carretero_2022.tree")
# check root age/tree height -- the raw file has a unit of 100 million-year
max(node.depth.edgelength(tree))
# transform unit to million years because I prefer so
tree$edge.length <- tree$edge.length * 100

trait <- read.csv("data/3_empirical/raw/PHYLACINE_1.2.1_Trait_data.csv")
trait <- trait %>%
  filter(!is.na(Mass.g),
         !is.na(Diet.Plant),
         !is.na(Diet.Vertebrate),
         !is.na(Diet.Invertebrate)) %>% 
  mutate(Binomial.1.2 = tolower(Binomial.1.2),
    herbivore   = ifelse(Diet.Plant >= 100, 1, 0),
    carnivore   = ifelse(Diet.Plant <= 0, 1, 0),
    omnivore    = ifelse(herbivore == 1 | carnivore == 1, 0, 1),
    p_omnivore  = ifelse(Diet.Plant < 100 & Diet.Plant > 50, 1, 0),
    np_omnivore = ifelse(Diet.Plant <= 50 & Diet.Plant > 0, 1, 0),
    diet3 = herbivore + omnivore * 2 + carnivore * 3 - 1,
    diet4 = herbivore + p_omnivore * 2 + np_omnivore * 3 + carnivore * 4 - 1,
    log_mass_kg = log(Mass.g / 1000)) %>% 
  select(Order.1.2, Binomial.1.2, log_mass_kg, diet3, diet4)

# species names found in both trait data and tree
tips_keep <- intersect(tree$tip.label, trait$Binomial.1.2)
# remove species not found in tree
tree <- keep.tip(tree, tips_keep)
# remove species not found in trait
trait <- trait %>% filter(Binomial.1.2 %in% tips_keep)

# save compiled traits as csv
write.csv(trait, "data/3_empirical/mammal_traits.csv")
# save tree
write.tree(tree, file="data/3_empirical/mammal_perMY_n3500.tre")

# save traits as nexus files
mammal_disc3 <- list()
mammal_disc4 <- list()
mammal_cont <- list()
for (i in 1:nrow(trait)){
  sp <- trait$Binomial.1.2[i]
  mammal_disc3[[sp]] <- trait$diet3[i]
  mammal_disc4[[sp]] <- trait$diet4[i]
  mammal_cont[[sp]] <- trait$log_mass_kg[i]
}
write.nexus.data(mammal_disc3, "data/3_empirical/mammal_diet3_n3500_Discrete.nex", format = "standard")
write.nexus.data(mammal_disc4, "data/3_empirical/mammal_diet4_n3500_Discrete.nex", format = "standard")
write.nexus.data(mammal_cont, "data/3_empirical/mammal_log_kg_n3500_Continuous.nex", format = "continuous")

# Subset the data because it's painful to wait for the mcmc to be done.
# But I actually did that in simulation studies already.
# So let's just use that one.

#tree_r500 <- read.tree("data/2_simulation/mammal_diet_perMY_n500.tre")
## always check if tree height is scaled as expected
#max(node.depth.edgelength(tree_r500))
#tips_r500 <- tolower(tree_r500$tip.label)
#trait_r500 <- trait %>% filter(Binomial.1.2 %in% tips_r500)

# Then I realized the tree I subsetted before contains taxa not found in this trait database so...
# I will switch back to plan A and subset another tree :)

# now, I wanna subset a tree that works for both manuscript and thesis
# and the thesis uses a different trait dataset as source
# let's search for taxa that intersect with that dataset as well

# 50-taxa tree for validation
set.seed(1234)
tree_r50 <- keep.tip(tree, sample(trait$Binomial.1.2, 50)) 
trait_r50 <- trait %>% filter(Binomial.1.2 %in% tree_r50$tip.label)
mammal_diet2_r50 <- list()
mammal_cont_r50 <- list()
for (i in 1:nrow(trait_r50)){
  sp <- trait_r50$Binomial.1.2[i]
  mammal_diet2_r50[[sp]] <- ifelse(trait_r50$diet4[i] <= 1, 0, 1)
  mammal_cont_r50[[sp]]  <- trait_r50$log_mass_kg[i]
}
tree_r50$edge.length <- tree_r50$edge.length / max(node.depth.edgelength(tree_r50))
# save characters
write.nexus.data(mammal_diet2_r50, "data/1_validation/mammal/mammal_diet2_r50_Discrete.nex", format = "standard")
write.nexus.data(mammal_cont_r50, "data/1_validation/mammal/mammal_log_kg_r50_Continuous.nex", format = "continuous")
# save tree
write.tree(tree_r50, file="data/1_validation/mammal/mammal_TH1_r50.tre")
simmaps <- make.simmap(tree_r50, x = mammal_diet2_r50 %>% unlist, model = "ARD", nsim=1)
write.simmap(simmaps, file = "data/1_validation/mammal/mammal_r50_simmap.tre")

#tree_r500 <- read.tree("data/3_empirical/mammal_perMY_r500.tre")

tips_r500 <- sample(trait_picky$Binomial.1.2, 500)
tree_r500 <- keep.tip(tree, tips_r500)
max(node.depth.edgelength(tree))
# Keep the sister branch of all other taxa so that the tree height doesn't change
while (max(node.depth.edgelength(tree)) < 202){
  tips_r500 <- sample(tips_keep, 500)
  tree_r500 <- keep.tip(tree, tips_r500)
}
trait_r500 <- trait %>% filter(Binomial.1.2 %in% tips_r500)

# save traits as nexus files
mammal_diet3_r500 <- mammal_diet4_r500 <- list()
mammal_cont_r500 <- list()
for (i in 1:nrow(trait_r500)){
  sp <- trait_r500$Binomial.1.2[i]
  mammal_diet3_r500[[sp]] <- trait_r500$diet3[i]
  mammal_diet4_r500[[sp]] <- trait_r500$diet4[i]
  mammal_cont_r500[[sp]]  <- trait_r500$log_mass_kg[i]
}
write.nexus.data(mammal_diet3_r500, "data/3_empirical/mammal_diet3_r500_Discrete.nex", format = "standard")
write.nexus.data(mammal_diet4_r500, "data/3_empirical/mammal_diet4_r500_Discrete.nex", format = "standard")
write.nexus.data(mammal_cont_r500, "data/3_empirical/mammal_log_kg_r500_Continuous.nex", format = "continuous")

# save tree
write.tree(tree_r500, file="data/3_empirical/mammal_perMY_r500.tre")
