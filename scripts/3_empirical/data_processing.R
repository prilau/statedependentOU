library(ape)
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
    herbivore = ifelse(Diet.Plant == 100, 1, 0),
    carnivore = ifelse(Diet.Vertebrate + Diet.Invertebrate == 100, 1, 0),
    omnivore = ifelse(herbivore == 1 | carnivore == 1, 0, 1),
    diet = herbivore + omnivore * 2 + carnivore * 3 - 1,
    log_mass_kg = log(Mass.g / 1000)) %>% 
  select(Binomial.1.2, log_mass_kg, diet, carnivore)

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
mammal_disc <- list()
mammal_cont <- list()
for (i in 1:nrow(trait)){
  sp <- trait$Binomial.1.2[i]
  mammal_disc[[sp]] <- trait$diet[i]
  mammal_cont[[sp]] <- trait$log_mass_kg[i]
}
write.nexus.data(mammal_disc, "data/3_empirical/mammal_diet_n3500_Discrete.nex", format = "standard")
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
trait_picky <- read.csv("data/3_empirical/raw/COMBINE_imputed_Soria_2021.csv") %>%
  mutate(Binomial.1.2 = tolower(gsub(" ", "_", phylacine_binomial))) %>% 
  filter(!duplicated(phylacine_binomial),
         !is.na(det_fruit),
         Binomial.1.2 %in% tips_keep) %>% 
  mutate(picky = ifelse(det_fruit >=90 | det_inv >=90 | det_vend>=90 | det_vect>=90 | det_vfish >=90 | det_vunk >=90 | det_scav >=90 | det_nect >=90 | det_seed >=90 | det_plantother >=90, 1, 0)) %>% 
  select(Binomial.1.2, picky)

set.seed(1234)
tree_r6 <- keep.tip(tree, sample(trait$Binomial.1.2, 6)) 
tips_r500 <- sample(trait_picky$Binomial.1.2, 500)
tree_r500 <- keep.tip(tree, tips_r500)
max(node.depth.edgelength(tree))
# Keep the sister branch of all other taxa so that the tree height doesn't change
while (max(node.depth.edgelength(tree)) < 202){
  tips_r500 <- sample(tips_keep, 500)
  tree_r500 <- keep.tip(tree, tips_r500)
}
trait_r500 <- trait %>% filter(Binomial.1.2 %in% tips_r500)
trait_picky_r500 <- trait_picky %>% filter(Binomial.1.2 %in% tips_r500)
trait_combined_r500 <- merge(trait_r500, trait_picky_r500, by="Binomial.1.2") %>% 
  mutate(picky_carn = carnivore + picky * 2,
         picky_herb = herbivore + picky * 2,
         picky_diet)

# save traits as nexus files
mammal_diet_r500 <- mammal_carn_r500 <- mammal_pick_r500 <- mammal_pica_r500 <- list()
mammal_cont_r500 <- list()
for (i in 1:nrow(trait_r500)){
  sp <- trait_r500$Binomial.1.2[i]
  mammal_diet_r500[[sp]] <- trait_combined_r500$diet[i]
  mammal_carn_r500[[sp]] <- trait_combined_r500$carnivore[i]
  mammal_pick_r500[[sp]] <- trait_combined_r500$picky[i]
  mammal_cont_r500[[sp]] <- trait_combined_r500$log_mass_kg[i]
  mammal_pica_r500[[sp]] <- trait_combined_r500$picky_carn[i]
}
write.nexus.data(mammal_diet_r500, "data/3_empirical/mammal_diet_r500_Discrete.nex", format = "standard")
write.nexus.data(mammal_carn_r500, "data/3_empirical/mammal_carnivory_r500_Discrete.nex", format = "standard")
write.nexus.data(mammal_pick_r500, "data/3_empirical/mammal_picky_r500_Discrete.nex", format = "standard")
write.nexus.data(mammal_pica_r500, "data/3_empirical/mammal_picky_carnivore_r500_Discrete.nex", format = "standard")
write.nexus.data(mammal_cont_r500, "data/3_empirical/mammal_log_kg_r500_Continuous.nex", format = "continuous")

# save tree
write.tree(tree_r500, file="data/3_empirical/mammal_perMY_r500.tre")
