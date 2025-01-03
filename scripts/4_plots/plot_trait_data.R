library(ape)
library(ggplot2)
library(tidyverse)

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
