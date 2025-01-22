library(phytools)
library(RevGadgets)
library(ggplot2)
library(tidyverse)
source("scripts/4_plots/revgadgets_StochMap.R")


#Plotting
tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")
t <- read.tree("data/3_empirical/mammal_perMY_r500.tre")

#################
# 3-state model #
#################
ase <- processAncStates("output/3_empirical/ase_r500_3State/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivorous", "1"="Omnivorous", "2"="Carnivorous"))
# produce the plot object, showing MAP states at nodes.
p1 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivorous"="#44AA99", "Omnivorous"="#ddcc77", "Carnivorous"="#882255"),
                       node_size = c(0.5, 2),
                       tip_states = TRUE,
                       #tip_states_size = c(0.01, 0.01),
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p1
ggsave(paste0("figures/3_empirical/ase_r500_3State/ase_map.pdf"), p1, width = 8, height = 6)


p2 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivorous"="#44AA99", "Omnivorous"="#ddcc77", "Carnivorous"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.2,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p2
ggsave(paste0("figures/3_empirical/ase_r500_3State/ase_pie.pdf"), p2, width = 8, height = 6)


# p3
simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_3State/marginal_character_run_1.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

colnames(stoc_map)[6] = "Herbivorous"
colnames(stoc_map)[7] = "Omnivorous"
colnames(stoc_map)[8] = "Carnivorous"

p3 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivorous"="#44AA99",
                               "Omnivorous"="#ddcc77",
                               "Carnivorous"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p3
ggsave(paste0("figures/3_empirical/ase_r500_3State/stoch_map.pdf"), p3, width = 8, height = 6)

#################
# 4-state model #
#################
ase <- processAncStates("output/3_empirical/ase_r500_4State/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivorous", "1"="Omnivorous (predominantly plants)",
                                       "2"="Omnivorous (predominantly non-plants)", "3"="Carnivorous"))
# produce the plot object, showing MAP states at nodes.
p4 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivorous"                   = "#44AA99",
                                      "Omnivorous (predominantly plants)"     = "#999933",
                                      "Omnivorous (predominantly non-plants)" = "#CC6677",
                                      "Carnivorous"                   = "#882255"),
                       node_size = c(0.5, 2),
                       tip_states = TRUE,
                       #tip_states_size = c(0.01, 0.01),
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p4
ggsave(paste0("figures/3_empirical/ase_r500_4State/ase_map.pdf"), p4, width = 8, height = 6)


p5 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivorous"="#44AA99",
                                      "Omnivorous (predominantly plants)"="#999933",
                                      "Omnivorous (predominantly non-plants)"="#CC6677",
                                      "Carnivorous"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.2,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p5
ggsave(paste0("figures/3_empirical/ase_r500_4State/ase_pie.pdf"), p5, width = 8, height = 6)


simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_4State/marginal_character_run_1.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2", "3"))

colnames(stoc_map)[6] = "Herbivorous"
colnames(stoc_map)[7] = "Omnivorous (predominantly plants)"
colnames(stoc_map)[8] = "Omnivorous (predominantly non-plants)"
colnames(stoc_map)[9] = "Carnivorous"

p6 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivorous"="#44AA99",
                               "Omnivorous (predominantly plants)"="#999933",
                               "Omnivorous (predominantly non-plants)"="#CC6677",
                               "Carnivorous"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p6
ggsave(paste0("figures/3_empirical/ase_r500_4State/stoch_map.pdf"), p6, width = 8, height = 6)



#p_all <- arrangeGrob(p1, p2, p4, p3, nrow = 2)
#ggsave("figures/3_empirical/ase_triState/ase_all.pdf", p_all, width = 16, height = 16)

#########################
# 3-state ordered model #
#########################
ase <- processAncStates("output/3_empirical/ase_r500_3StateOrdered/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivorous", "1"="Omnivorous", "2"="Carnivorous"))
# produce the plot object, showing MAP states at nodes.
p1 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivorous"="#44AA99", "Omnivorous"="#ddcc77", "Carnivorous"="#882255"),
                       node_size = c(0.5, 2),
                       tip_states = TRUE,
                       #tip_states_size = c(0.01, 0.01),
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p1
ggsave(paste0("figures/3_empirical/ase_r500_3StateOrdered/ase_map.pdf"), p1, width = 8, height = 6)


p2 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivorous"="#44AA99", "Omnivorous"="#ddcc77", "Carnivorous"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.2,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p2
ggsave(paste0("figures/3_empirical/ase_r500_3StateOrdered/ase_pie.pdf"), p2, width = 8, height = 6)


# p3
simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_3StateOrdered/marginal_character_run_1.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

colnames(stoc_map)[6] = "Herbivorous"
colnames(stoc_map)[7] = "Omnivorous"
colnames(stoc_map)[8] = "Carnivorous"

p3 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivorous"="#44AA99",
                               "Omnivorous"="#ddcc77",
                               "Carnivorous"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p3
ggsave(paste0("figures/3_empirical/ase_r500_3StateOrdered/stoch_map.pdf"), p3, width = 8, height = 6)









#############
# augmented #
#############
#  3-state  #
#############
ase <- processAncStates("output/3_empirical/aug_r500_3StateOrderedModel/anc_states_run_5.tre",
                        state_labels=c("0"="Herbivorous", "1"="Omnivorous", "2"="Carnivorous"))
# produce the plot object, showing MAP states at nodes.
p7 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivorous"="#44AA99", "Omnivorous"="#ddcc77", "Carnivorous"="#882255"),
                       node_size = c(0.5, 2),
                       tip_states = TRUE,
                       #tip_states_size = c(0.01, 0.01),
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p1
ggsave(paste0("figures/3_empirical/aug_r500_3StateOrderedModel/ase_map.pdf"), p7, width = 8, height = 6)


p8 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivorous"="#44AA99", "Omnivorous"="#ddcc77", "Carnivorous"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.2,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
#p2
ggsave(paste0("figures/3_empirical/aug_r500_3StateOrderedModel/ase_pie.pdf"), p8, width = 8, height = 6)


# p3
simmaps <- read.simmap("output/3_empirical/aug_r500_3StateOrderedModel/augch_run_5.trees", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

colnames(stoc_map)[6] = "Herbivorous"
colnames(stoc_map)[7] = "Omnivorous"
colnames(stoc_map)[8] = "Carnivorous"

p9 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivorous"="#44AA99",
                               "Omnivorous"="#ddcc77",
                               "Carnivorous"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p3
ggsave(paste0("figures/3_empirical/aug_r500_3StateOrderedModel/stoch_map.pdf"), p9, width = 8, height = 6)
