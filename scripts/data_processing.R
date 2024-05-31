library(ape)
library(tidyverse)
library(phytools)

diet <- read.csv("data/3_empirical/mammal_diet_Discrete.txt", sep = "\t")
diet <- diet %>% filter(Diet != "?") %>% 
  rename(species = Species.Name...Wilson...Reeder.2005)

foot <- read.nexus.data("Desktop/sdOU_local/IRT3/data/kubo_2019/disc_footpostures_wo_footpads.nex")

tree <- read.nexus("data/3_empirical/mammal_2019.trees")[[1]]
old <- strsplit(tree$tip.label, "_")
for (i in 1:length(old)){
  new <- paste0(old[[i]][1], "_", old[[i]][2])
  tree$tip.label[i] <- new
}

#tree <- read.nexus("data/3_empirical/mammal_2022.tree")
#tree$tip.label <- str_to_title(tree$tip.label)

body2009 <- read.csv("data/3_empirical/mammal_bodysize_2009.tsv", sep="\t")
body2009$species <- str_replace(body2009$MSW93_Binomial, " ", "_")
body2009 <- body2009 %>% select(c(species, X5.1_AdultBodyMass_g)) %>% 
  rename(body_mass = X5.1_AdultBodyMass_g) %>% 
  filter(body_mass > 0)
  

#body2022 <- read.csv("data/3_empirical/mammal_bodysize_2022.csv")
#body2022$species <- str_replace(body2022$species, " ", "_")
#body2022 <- body2022 %>% filter(body.mass...units == "kg") %>% 
#  mutate(body_mass = body.mass * 1000) %>% 
#  select(c(species, body_mass))

match_foot <- intersect(intersect(body2009$species, names(foot)), tree$tip.label)
match_diet <- intersect(intersect(body2009$species, diet$species), tree$tip.label)

body2009_foot <- body2009[which(body2009$species %in% match_foot),]
body2009_diet <- body2009[which(body2009$species %in% match_diet),]
body_nex_foot <- list()
body_nex_diet <- list()
for (i in 1:nrow(body2009_foot)){
  sp <- body2009_foot[i,1]
  body_nex_foot[[sp]] <- log(body2009_foot[i,2]/1000)
}
for (i in 1:nrow(body2009_diet)){
  sp <- body2009_diet[i,1]
  body_nex_diet[[sp]] <- log(body2009_diet[i,2]/1000)
}

write.nexus.data(body_nex_diet, "data/3_empirical/mammal_diet_Continuous.nex", format = "continuous")
write.nexus.data(body_nex_foot, "data/3_empirical/mammal_foot_Continuous.nex", format = "continuous")


diet <- diet[which(diet$species %in% match_diet),]
diet_nex <- list()
for (i in 1:nrow(diet)){
  sp <- diet[i,1]
  diet_nex[[sp]] <- diet[i,2]
}
foot <- foot[which(names(foot) %in% match_foot)]

write.nexus.data(diet_nex, "data/3_empirical/mammal_diet_Discrete.nex", format = "standard")
write.nexus.data(foot, "data/3_empirical/mammal_foot_Discrete.nex", format = "standard")


tree_diet <- drop.tip(tree, tree$tip.label[-which(tree$tip.label %in% match_diet)])
tree_foot <- drop.tip(tree, tree$tip.label[-which(tree$tip.label %in% match_foot)])

write.tree(tree_diet, "data/3_empirical/mammal_diet.tre")
write.tree(tree_foot, "data/3_empirical/mammal_foot.tre")

max(node.depth.edgelength(tree_foot))
