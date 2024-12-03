library(ape)
library(tidyverse)

# read in trait data
# remove rows without the continuous trait of interest
# change formatting of species name according to tree tip labels
df <- read.csv("data/3_empirical/raw/COMBINE_imputed_Soria_2021.csv") %>% 
  filter(!is.na(adult_mass_g),
         !duplicated(phylacine_binomial)) %>% 
  mutate(binomial_name = tolower(gsub(" ", "_", phylacine_binomial)), 
         log_body_mass = log(adult_mass_g / 1000))

# plot distribution
ggplot(data=df, aes(x=log_body_mass)) +
  geom_histogram()


# code three binary characters
# herbivore defined as consuming 100% plant materials
# carnivore defined as consuming 100% invertebrate and/or vertebrate materials
# picky eater defined as consuming only one type of material
# aquatic defined as living in freshwater or marine habitat
df_hpa <- df %>% 
  filter(!is.na(dphy_plant),
         !is.na(dphy_invertebrate),
         !is.na(dphy_vertebrate),
         !is.na(det_fruit)) %>% 
  mutate(picky = ifelse(det_fruit >=90 | det_inv >=90 | det_vend>=90 | det_vect>=90 | det_vfish >=90 | det_vunk >=90 | det_scav >=90 | det_nect >=90 | det_seed >=90 | det_plantother >=90, 1, 0),
         herbivore = ifelse(dphy_plant == 100, 1, 0),
         carnivore = ifelse(dphy_invertebrate + dphy_vertebrate == 100, 1, 0),
         picky_herb = picky * 1 + herbivore * 2,
         picky_carn = picky * 1 + carnivore * 2)

# read in tree
tree <- read.nexus("data/3_empirical/raw/4705sp_mammal-time_Álvarez-Carretero_2022.tree")
# check root age
max(node.depth.edgelength(tree))
# transform unit to million years
tree$edge.length <- tree$edge.length * 100
# species names found in both trait data and tree
tips_keep <- intersect(tree$tip.label, df_hpa$binomial_name)
# remove species not found in trait data
tree <- keep.tip(tree, tips_keep)
# remove species not found in tree
df_hpa <- df_hpa %>% filter(binomial_name %in% tips_keep)


# plot characters on tree
# script modified from Boyko et al. (2023)
root_age <- max(node.depth.edgelength(tree))
plot(tree, show.tip.label = FALSE, x.lim = c(0, root_age + 0.2 * root_age))
disc_chars <- c("herbivore", "carnivore", "picky")
begin = root_age + (0.005 * root_age)
color_0 = c('#CCDDAA', '#FFCCCC', '#BBCCEE')
color_1 = c('#225522', '#663333', '#222255')
j=1
for(char in disc_chars){
  end = begin + 5
  for(i in 1:length(tree$tip.label)){
    lines(list(x = c(begin, end), y = c(i,i)), 
               col = ifelse(df_hpa[[char]][which(df_hpa$binomial_name == tree$tip.label[i])] == 0,
                            color_0[j], color_1[j]),
               lwd = 1)
  }
  j=j+1
  begin = end + 1
}
legend("bottomleft", legend = c("Herbivore","Carnivore", "Picky eater"),
       pch=15, col = color_1[1:length(disc_chars)])

# plot distribution of character states (not phylo-corrected)
ggplot(data=df_hpa, aes(x=log_body_mass, fill=factor(herbivore))) +
  geom_histogram()
ggplot(data=df_hpa, aes(x=log_body_mass, fill=factor(carnivore))) +
  geom_histogram()
ggplot(data=df_hpa, aes(x=log_body_mass, fill=factor(picky))) +
  geom_histogram()
ggplot(data=df_hpa, aes(x=log_body_mass, fill=factor(picky_herb))) +
  geom_histogram() +
  scale_fill_discrete()
ggplot(data=df_hpa, aes(x=log_body_mass, fill=factor(picky_carn))) +
  geom_histogram()






# save traits as nexus files
mammal_disc_binary <- list()
mammal_disc_combined <- list()
mammal_log_kg <- list()
for (i in 1:nrow(df_hpa)){
  sp <- df_hpa$binomial_name[i]
  mammal_disc_binary[[sp]] <- c(df_hpa$herbivore[i], df_hpa$carnivore[i], df_hpa$picky[i])
  mammal_disc_combined[[sp]] <- c(df_hpa$picky_herb[i], df_hpa$picky_carn[i])
  mammal_log_kg[[sp]] <- df_hpa$log_body_mass[i]
}
write.nexus.data(mammal_disc_binary, "data/3_empirical/mammal_diet_binary_Discrete.nex", format = "standard")
write.nexus.data(mammal_disc_combined, "data/3_empirical/mammal_diet_combined_Discrete.nex", format = "standard")
write.nexus.data(mammal_log_kg, "data/3_empirical/mammal_log_kg_Continuous.nex", format = "continuous")

# save tree
write.tree(tree, file="data/3_empirical/mammal_perMY_n2955.tre")
