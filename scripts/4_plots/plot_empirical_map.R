library(phytools)
library(RevGadgets)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
library(latex2exp)
source("scripts/4_plots/revgadgets_StochMap.R")


#Plotting
tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")
t <- read.tree("data/3_empirical/mammal_perMY_r500.tre")

######################
# stochastic mapping #
######################
#########################
# 3-state ordered model #
#########################
ase <- processAncStates("output/3_empirical/ase_r500_3StateOrdered/anc_states_run_4.tre",
                        state_labels=c("0"="Herbivore", "1"="Omnivore", "2"="Carnivore"))
# produce the plot object, showing MAP states at nodes.
#p1 <- plotAncStatesMAP(t = ase,
#                       #tip_labels_offset = 0.5,
#                       tip_labels = FALSE,
#                       node_color_as = "state",
#                       node_color = c("Herbivore"="#44AA99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
#                       node_size = c(0.3, 1.2),
#                       tip_states = TRUE,
#                       tip_states_size = 0.3,
#                       #tip_states_shape = 1,
#                       state_transparency = 0.7,
#                       tree_layout = "circular",
#                       #tip_labels_size = 0.5
#                       tree_color = "#bbbbbb",
#                       tree_linewidth = 0.25) +
#  # modify legend location using ggplot2
#  theme(legend.position.inside = c(0.6,0.81))
#  #theme(legend.position = "none")
#p1
#ggsave(paste0("figures/3_empirical/ase_r500_3StateOrdered/run_4_node_map.pdf"), p1, width = 8, height = 6)


p1 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivore"="#44AA99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.3,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  #theme(legend.position.inside = c(0.6,0.81))
  theme(legend.position = "none")

#p1


simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_3StateOrdered/marginal_character_run_4.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

colnames(stoc_map)[6] = "Herbivore"
colnames(stoc_map)[7] = "Omnivore"
colnames(stoc_map)[8] = "Carnivore"

p2 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivore"="#44AA99",
                               "Omnivore"="#ddcc77",
                               "Carnivore"="#882255")) +
  theme(legend.position = "none"
        #legend.position.inside = c(0.6,0.81)
        )
  #theme(legend.position = "none")
#p2


#################
# 4-state model #
#################
ase <- processAncStates("output/3_empirical/ase_r500_4State/anc_states_run_2.tre",
                        state_labels=c("0"="Herbivore", "1"="Omnivore (>50% plants)",
                                       "2"="Omnivore (<=50% plants)", "3"="Carnivore"))
# produce the plot object, showing MAP states at nodes.
#p5 <- plotAncStatesMAP(t = ase,
#                       #tip_labels_offset = 0.5,
#                       tip_labels = FALSE,
#                       node_color_as = "state",
#                       node_color = c("Herbivore"                  = "#44AA99",
#                                      "Omnivore (>50% plants)"     = "#999933",
#                                      "Omnivore (<=50% plants)"    = "#CC6677",
#                                      "Carnivore"                  = "#882255"),
#                       node_size = c(0.5, 2),
#                       tip_states = TRUE,
#                       #tip_states_size = c(0.01, 0.01),
#                       #tip_states_shape = 1,
#                       state_transparency = 0.7,
#                       tree_layout = "circular",
#                       #tip_labels_size = 0.5
#                       tree_color = "#bbbbbb",
#                       tree_linewidth = 0.25) +
#  # modify legend location using ggplot2
#  theme(legend.position = "none"
#    #legend.position.inside = c(0.6,0.81)
#    )
#p4
#ggsave(paste0("figures/3_empirical/ase_r500_4State/run_2_node_map.pdf"), p4, width = 8, height = 6)


p5 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivore"="#44AA99",
                                      "Omnivore (>50% plants)"="#999933",
                                      "Omnivore (<=50% plants)"="#CC6677",
                                      "Carnivore"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.3,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  
  #theme(legend.position.inside = c(0.6,0.81))
  theme(legend.position = "none")

#p5
#ggsave(paste0("figures/3_empirical/ase_r500_4State/run_2_node_pie.pdf"), p5, width = 8, height = 6)


simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_4State/marginal_character_run_2.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2", "3"))

colnames(stoc_map)[6] = "Herbivore"
colnames(stoc_map)[7] = "Omnivore (>50% plants)"
colnames(stoc_map)[8] = "Omnivore (<=50% plants)"
colnames(stoc_map)[9] = "Carnivore"

p6 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivore"="#44AA99",
                               "Omnivore (>50% plants)"="#999933",
                               "Omnivore (<=50% plants)"="#CC6677",
                               "Carnivore"="#882255")) +
  theme(legend.position = "none"
        #legend.position.inside = c(0.6,0.81)
        )
#p6
#ggsave(paste0("figures/3_empirical/ase_r500_4State/run_2_stoch_map.pdf"), p6, width = 8, height = 6)


#p_all <- arrangeGrob(p1, p2, p4, p3, nrow = 2)
#ggsave("figures/3_empirical/ase_triState/ase_all.pdf", p_all, width = 16, height = 16)


#############
# augmented #
#############
#  3-state  #
#############
ase <- processAncStates("output/3_empirical/aug_r500_3StateOrderedModel/anc_states_run_5.tre",
                        state_labels=c("0"="Herbivore", "1"="Omnivore", "2"="Carnivore"))
# produce the plot object, showing MAP states at nodes.
p7 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivore"="#44AA99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
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
                       pie_colors = c("Herbivore"="#44AA99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
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

colnames(stoc_map)[6] = "Herbivore"
colnames(stoc_map)[7] = "Omnivore"
colnames(stoc_map)[8] = "Carnivore"

p9 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivore"="#44AA99",
                               "Omnivore"="#ddcc77",
                               "Carnivore"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p3
ggsave(paste0("figures/3_empirical/aug_r500_3StateOrderedModel/stoch_map.pdf"), p9, width = 8, height = 6)



#############
# missing   #
#############
#  2-state  #
#############
ase <- processAncStates("output/3_empirical/sdOU_r500_missingStateModel/anc_states_nstate_2_run_6.log",
                        state_labels=c("0"="Large", "1"="Small"))

p1 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Large"="#44AA99", "Small"="#882255", "2"="grey0", "3"="grey1"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.5,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  #theme(legend.position.inside = c(0.6,0.81))
  theme(legend.position = "none")
p1


# produce the plot object, showing MAP states at nodes.
#p1 <- plotAncStatesMAP(t = ase,
#                       #tip_labels_offset = 0.5,
#                       tip_labels = FALSE,
#                       node_color_as = "state",
#                       node_color = c("Large"="#44AA99", "Small"="#882255"),
#                       node_size = c(0.5, 2),
#                       tip_states = TRUE,
#                       tip_states_size = 1,
#                       #tip_states_shape = 1,
#                       state_transparency = 0.7,
#                       tree_layout = "circular",
#                       #tip_labels_size = 0.5
#                       tree_color = "#bbbbbb",
#                       tree_linewidth = 0.25) +
#  # modify legend location using ggplot2
#  theme(legend.position = "none"
#        #legend.position.inside = c(0.6,0.81)
#        )

#simmaps <- read.simmap("output/3_empirical/sdOU_r500_missingStateModel/augch_nstate_2_run_6.trees", format="phylip")
#stoc_map <- processStochMaps(tree,
#                             simmap = simmaps,
#                             states=c("0", "1", "2"))
#
#p9 <- plotStochMaps(tree, maps=stoc_map,
#                    tip_labels = FALSE,
#                    tree_layout = "circular",
#                    line_width=0.25,
#                    color_by = "MAP",
#                    colors = c("0"="#44AA99",
#                               "1"="#882255",
#                               "2"="#ddcc77")) +
#  theme(legend.position.inside = c(0.6,0.81))
##p3






################
#  2-state sim #
################
ase <- processAncStates("output/2_simulation/missing_state/anc_states_nstate_2_run_2.log",
                        state_labels=c("0"="0", "1"="1", "2"="2"))
# produce the plot object, showing MAP states at nodes.
p7 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("0"="#44AA99", "1"="#ddcc77", "2"="#882255"),
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
ggsave(paste0("figures/2_simulation/missing_state/nstate_2_ase_map.pdf"), p7, width = 8, height = 6)


p8 <- plotAncStatesPie(t = ase,
                       pie_colors = c("0"="#44AA99", "1"="#ddcc77", "2"="#882255"),
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
ggsave(paste0("figures/2_simulation/missing_state/nstate_2_ase_pie.pdf"), p8, width = 8, height = 6)


# p3
simmaps <- read.simmap("output/2_simulation/missing_state/sim_1/augch_nstate_2_run_1.trees", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

p9 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("0"="#44AA99",
                               "1"="#ddcc77",
                               "2"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p3
ggsave(paste0("figures/2_simulation/missing_state/nstate_2_stoch_map.pdf"), p9, width = 8, height = 6)


#########
# aug   #
#########
#  3-state  #
#############
dir_in = "output/3_empirical/aug_r500_3StateOrderedModel/"
dir_out = "figures/3_empirical/aug_r500_3StateOrderedModel/"
files <- list.files(dir_in, pattern = "anc")

for (file in files){
  run <- strsplit(strsplit(file, "_")[[1]][4], "\\.")[[1]][1]
  ase <- processAncStates(paste0(dir_in, file), state_labels=c("0"="Herbivore",
                                                               "1"="Omnivore",
                                                               "2"="Carnivore"))
  p1 <- plotAncStatesMAP(t = ase,
                         #tip_labels_offset = 0.5,
                         tip_labels = FALSE,
                         node_color_as = "state",
                         node_color = c("Herbivore"="#44AA99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
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
  
  ggsave(paste0(dir_out, "run_", run, "_node_map.pdf"), p1, width = 8, height = 6)
  
  
  p2 <- plotAncStatesPie(t = ase,
                         pie_colors = c("Herbivore"="#44AA99",
                                        "Omnivore"="#ddcc77",
                                        "Carnivore"="#882255",
                                        "4"="#000000"),
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
  ggsave(paste0(dir_out, "run_", run, "_node_pie.pdf"), p2, width = 8, height = 6)
}

augchs <- list.files(dir_in, pattern = ".trees")

for (augch in augchs){
  run <- strsplit(strsplit(augch, "_")[[1]][3], "\\.")[[1]][1]
  simmaps <- read.simmap(paste0(dir_in, augch), format="phylip")
  stoc_map <- processStochMaps(tree,
                               simmap = simmaps,
                               states=c("0", "1", "2"))
  
  colnames(stoc_map)[6] = "Herbivore"
  colnames(stoc_map)[7] = "Omnivore"
  colnames(stoc_map)[8] = "Carnivore"
  
  p3 <- plotStochMaps(tree, maps=stoc_map,
                      tip_labels = FALSE,
                      tree_layout = "circular",
                      line_width=0.25,
                      color_by = "MAP",
                      colors = c("Herbivore"="#44AA99",
                                 "Omnivore"="#ddcc77",
                                 "Carnivore"="#882255")) +
    theme(legend.position.inside = c(0.6,0.81))
  #p3
  ggsave(paste0(dir_out, "run_", run, "_stoch_map.pdf"), p3, width = 8, height = 6)
}

#############
#  4-state  #
#############
dir_in = "output/3_empirical/aug_r500_4StateModel/"
dir_out = "figures/3_empirical/aug_r500_4StateModel/"
files <- list.files(dir_in, pattern = "anc")

for (file in files){
  run <- strsplit(strsplit(file, "_")[[1]][4], "\\.")[[1]][1]
  ase <- processAncStates(paste0(dir_in, file), state_labels=c("0"="Herbivore",
                                                               "1"="Omnivore (>50% plants)",
                                                               "2"="Omnivore (≤50% plants)",
                                                               "3"="Carnivore"))
  p1 <- plotAncStatesMAP(t = ase,
                         #tip_labels_offset = 0.5,
                         tip_labels = FALSE,
                         node_color_as = "state",
                         node_color = c("Herbivore"="#44AA99",
                                        "Omnivore (>50% plants)"="#999933",
                                        "Omnivore (≤50% plants)"="#cc6677",
                                        "Carnivore"="#882255"),
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
  
  ggsave(paste0(dir_out, "run_", run, "_node_map.pdf"), p1, width = 8, height = 6)
  
  
  p2 <- plotAncStatesPie(t = ase,
                         pie_colors = c("Herbivore"="#44AA99",
                                        "Omnivore (>50% plants)"="#999933",
                                        "Omnivore (≤50% plants)"="#cc6677",
                                        "Carnivore"="#882255"),
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
  ggsave(paste0(dir_out, "run_", run, "_node_pie.pdf"), p2, width = 8, height = 6)
}

augchs <- list.files(dir_in, pattern = ".trees")

for (augch in augchs){
  run <- strsplit(strsplit(augch, "_")[[1]][3], "\\.")[[1]][1]
  simmaps <- read.simmap(paste0(dir_in, augch), format="phylip")
  stoc_map <- processStochMaps(tree,
                               simmap = simmaps,
                               states=c("0", "1", "2", "3"))
  
  colnames(stoc_map)[6] = "Herbivore"
  colnames(stoc_map)[7] = "Omnivore (>50% plants)"
  colnames(stoc_map)[8] = "Omnivore (≤50% plants)"
  colnames(stoc_map)[9] = "Carnivore"
  
  p3 <- plotStochMaps(tree, maps=stoc_map,
                      tip_labels = FALSE,
                      tree_layout = "circular",
                      line_width=0.25,
                      color_by = "MAP",
                      colors = c("Herbivore"="#44AA99",
                                 "Omnivore (>50% plants)"="#999933",
                                 "Omnivore (≤50% plants)"="#cc6677",
                                 "Carnivore"="#882255")) +
    theme(legend.position.inside = c(0.6,0.81))
  #p3
  ggsave(paste0(dir_out, "run_", run, "_stoch_map.pdf"), p3, width = 8, height = 6)
}


########
# sdou #
########
tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")
#############
#  3-state  #
#############
dir_in = "output/3_empirical/sdOU_r500_3StateOrderedModel/"
dir_out = "figures/3_empirical/sdOU_r500_3StateOrderedModel/"
file <- list.files(dir_in, pattern = "anc")[1]

ase <- processAncStates(paste0(dir_in, file), state_labels=c("0"="Herbivore",
                                                             "1"="Omnivore",
                                                             "2"="Carnivore"))
p0 <- plotAncStatesMAP(t = ase,
                         #tip_labels_offset = 0.5,
                         tip_labels = FALSE,
                         node_color_as = "state",
                         node_color = c("Herbivore"="#44AA99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
                         node_size = c(0.3, 1.2),
                         tip_states = TRUE,
                         tip_states_size = 0.3,
                         #tip_states_shape = 1,
                         state_transparency = 0.7,
                         tree_layout = "circular",
                         #tip_labels_size = 0.5
                         tree_color = "#bbbbbb",
                         tree_linewidth = 0.25) +
    # modify legend location using ggplot2
    theme(legend.position = "none"
      #legend.position.inside = c(0.6,0.81)
      )
#p1
#ggsave(paste0(dir_out, "run_", run, "_node_map.pdf"), p1, width = 8, height = 6)
  
  
p3 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivore"="#44AA99",
                                      "Omnivore"="#ddcc77",
                                      "Carnivore"="#882255",
                                      "4"="#000000"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.3,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
    # modify legend location using ggplot2
    #theme(legend.position.inside = c(0.6,0.81))
    theme(legend.position = "none")
#p3
#ggsave(paste0(dir_out, "run_", run, "_node_pie.pdf"), p2, width = 8, height = 6)


augch <- list.files(dir_in, pattern = ".trees")[1]

simmaps <- read.simmap(paste0(dir_in, augch), format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

colnames(stoc_map)[6] = "Herbivore"
colnames(stoc_map)[7] = "Omnivore"
colnames(stoc_map)[8] = "Carnivore"

p4 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivore"="#44AA99",
                               "Omnivore"="#ddcc77",
                               "Carnivore"="#882255")) +
    theme(legend.position = "none"
      #legend.position.inside = c(0.6,0.81)
      )
p4
#ggsave(paste0(dir_out, "run_", run, "_stoch_map.pdf"), p3, width = 8, height = 6)

#legend3_color <- get_legend2(p1 + theme(legend.position = "left",
#                                  legend.box.margin = margin(0, 0, 0, 12))
#                                + guides(size = "none",
#                                         color = guide_legend(override.aes = list(size = 3),
#                                                              title='Diet'),
#                                         fill=guide_legend(title='Diet'))
#                                + scale_fill_discrete(name = "Diet"))
#legend3_size <- get_legend2(p1 + theme(legend.position = "right",
#                                        legend.box.margin = margin(0, 0, 0, 12))
#                             + guides(color = "none"))

legend3 <- get_legend2(p0 + theme(legend.position = "right",
                                  legend.box.margin = margin(0, 0, 0, 12))
                                + guides(size = "none",
                                         color = guide_legend(override.aes = list(size = 3),
                                                              title='Diet'),
                                         fill=guide_legend(title='Diet'))
                                + scale_fill_discrete(name = "Diet"))

ase3_sdou3_row1 <- cowplot::plot_grid(p1, p2, ncol=2, labels=c("(a)", "(b)"))
ase3_sdou3_row2 <- cowplot::plot_grid(p3, p4, ncol=2, labels=c("(c)", "(d)"))
ase3_sdou3_row3 <- cowplot::plot_grid(legend3, NULL, rel_widths = c(0.2,0.8))

ase_sdou3_all <- cowplot::plot_grid(ase3_sdou3_row1, ase3_sdou3_row2, ase3_sdou3_row3, nrow=3, rel_heights = c(1,1,0.3))

ggsave("figures/3_empirical/sdOU_r500_3StateOrderedModel/char_hist_ase_sdou.pdf", ase_sdou3_all, width = 7.4, height = 8.3, units = "in")


#plot in two rows
#maps_sdou3_row1 <- cowplot::plot_grid(p1, p2, ncol=2, labels=c("(a)", "(b)"))
#maps_sdou3_row2 <- cowplot::plot_grid(p3, NULL, legend3_size, legend3_color, NULL, ncol=5,
#                                      rel_widths = c(0.5,0.1,0.15,0.15,0.1),
#                                      labels=c("(c)", "", "", "", ""))
#maps_sdou3_tworows <- cowplot::plot_grid(maps_sdou3_row1, maps_sdou3_row2, ncol=1)
#
#
#ggsave("figures/3_empirical/sdOU_r500_3StateOrderedModel/char_hist3.pdf", maps_sdou3_tworows, width = 7.4, height = 7.4, units = "in")



#############
#  4-state  #
#############
dir_in = "output/3_empirical/sdOU_r500_4StateModel/"
dir_out = "figures/3_empirical/sdOU_r500_4StateModel/"

file <- list.files(dir_in, pattern = "anc")[1]

ase <- processAncStates(paste0(dir_in, file), state_labels=c("0"="Herbivore",
                                                             "1"="Omnivore (>50%$ plants)",
                                                             "2"="Omnivore (<=50% plants)",
                                                             "3"="Carnivore"))
p0 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivore"="#44AA99",
                                      "Omnivore (>50%$ plants)"="#999933",
                                      "Omnivore (<=50% plants)"="#cc6677",
                                      "Carnivore"="#882255"),
                       node_size = c(0.3, 1.2),
                       tip_states = TRUE,
                       tip_states_size = 0.3,
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position = "none"
        #legend.position.inside = c(0.6,0.81)
  )
#p1
#ggsave(paste0(dir_out, "run_", run, "_node_map.pdf"), p1, width = 8, height = 6)


p7 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivore"="#44AA99",
                                      "Omnivore (>50%$ plants)"="#999933",
                                      "Omnivore (<=50% plants)"="#cc6677",
                                      "Carnivore"="#882255"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.3,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  #theme(legend.position.inside = c(0.6,0.81))
  theme(legend.position = "none")
#p7
#ggsave(paste0(dir_out, "run_", run, "_node_pie.pdf"), p2, width = 8, height = 6)


augch <- list.files(dir_in, pattern = ".trees")[1]

simmaps <- read.simmap(paste0(dir_in, augch), format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2", "3"))

colnames(stoc_map)[6] = "Herbivore"
colnames(stoc_map)[7] = "Omnivore (>50% plants)"
colnames(stoc_map)[8] = "Omnivore (<=50% plants)"
colnames(stoc_map)[9] = "Carnivore"

p8 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivore"="#44AA99",
                               "Omnivore (>50% plants)"="#999933",
                               "Omnivore (<=50% plants)"="#cc6677",
                               "Carnivore"="#882255")) +
  theme(legend.position = "none"
        #legend.position.inside = c(0.6,0.81)
  )
p8


#legend4_color <- get_legend2(p1 + theme(legend.position = "left",
#                                        legend.box.margin = margin(0, 0, 0, 12))
#                             + scale_color_manual(values=c("#44aa99", "#999933", "#CC6677", "#882255"), 
#                                                name="Diet",
#                                                labels=c("Herbivore", TeX("Omnivore ($>50\\%$ plants)"),
#                                                         TeX("Omnivore ($\\leq 50\\%$ plants)"), "Carnivore"))
#                             + guides(size = "none",
#                                      color = guide_legend(override.aes = list(size = 3),
#                                                           title='Diet'),
#                                      fill=guide_legend(title='Diet'))
#                             + scale_fill_discrete(name = "Diet"))
#legend4_size <- get_legend2(p1 + theme(legend.position = "right",
#                                       legend.box.margin = margin(0, 0, 0, 12))
#                            + guides(color = "none"))
#
##plot in two rows
#maps_sdou4_row1 <- cowplot::plot_grid(p1, p2, ncol=2, labels=c("(a)", "(b)"))
#maps_sdou4_row2 <- cowplot::plot_grid(p3, legend4_size, legend4_color, ncol=3,
#                                      rel_widths = c(0.5,0.2,0.3),
#                                      labels=c("(c)", "", ""))
#maps_sdou4_tworows <- cowplot::plot_grid(maps_sdou4_row1, maps_sdou4_row2, ncol=1)

legend4 <- get_legend2(p0 + theme(legend.position = "left",
                                        legend.box.margin = margin(0, 0, 0, 12))
                             + scale_color_manual(values=c("#44aa99", "#999933", "#CC6677", "#882255"), 
                                                name="Diet",
                                                labels=c("Herbivore", TeX("Omnivore ($>50\\%$ plants)"),
                                                         TeX("Omnivore ($\\leq 50\\%$ plants)"), "Carnivore"))
                             + guides(size = "none",
                                      color = guide_legend(override.aes = list(size = 3),
                                                           title='Diet'),
                                      fill=guide_legend(title='Diet')))

ase4_sdou4_row1 <- cowplot::plot_grid(p5, p6, ncol=2, labels=c("(a)", "(b)"))
ase4_sdou4_row2 <- cowplot::plot_grid(p7, p8, ncol=2, labels=c("(c)", "(d)"))
ase4_sdou4_row3 <- cowplot::plot_grid(legend4, NULL, rel_widths = c(0.3,0.7))

ase4_sdou4_all <- cowplot::plot_grid(ase4_sdou4_row1, ase4_sdou4_row2, ase4_sdou4_row3, nrow=3, rel_heights = c(1,1,0.3))


ggsave("figures/3_empirical/sdOU_r500_4StateModel/char_hist_ase_sdou.pdf", ase4_sdou4_all, width = 7.4, height = 8.3, units = "in")


################
# hidden-state #
################
dir_in = "output/3_empirical/sdOU_r500_hiddenStateModel/"
dir_out = "figures/3_empirical/sdOU_r500_hiddenStateModel/"
files <- list.files(dir_in, pattern = "anc")

for (file in files){
  run <- strsplit(strsplit(file, "_")[[1]][4], "\\.")[[1]][1]
  ase <- processAncStates(paste0(dir_in, file), state_labels=c("0"="Herbivore",
                                                               "1"="Omnivore (>50% plants)",
                                                               "2"="Omnivore (≤50% plants)",
                                                               "3"="Carnivore"))
p1 <- plotAncStatesMAP(t = ase,
                         #tip_labels_offset = 0.5,
                         tip_labels = FALSE,
                         node_color_as = "state",
                         node_color = c("Herbivore"="#44AA99",
                                        "Omnivore (>50% plants)"="#999933",
                                        "Omnivore (≤50% plants)"="#cc6677",
                                        "Carnivore"="#882255"),
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
  
  ggsave(paste0(dir_out, "run_", run, "_node_map.pdf"), p1, width = 8, height = 6)
  
  
  p2 <- plotAncStatesPie(t = ase,
                         pie_colors = c("Herbivore"="#44AA99",
                                        "Omnivore (>50% plants)"="#999933",
                                        "Omnivore (≤50% plants)"="#cc6677",
                                        "Carnivore"="#882255"),
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
  ggsave(paste0(dir_out, "run_", run, "_node_pie.pdf"), p2, width = 8, height = 6)
}

augchs <- list.files(dir_in, pattern = ".trees")

for (augch in augchs){
  run <- strsplit(strsplit(augch, "_")[[1]][3], "\\.")[[1]][1]
  simmaps <- read.simmap(paste0(dir_in, augch), format="phylip")
  stoc_map <- processStochMaps(tree,
                               simmap = simmaps,
                               states=c("0", "1", "2", "3"))
  
  colnames(stoc_map)[6] = "Herbivore"
  colnames(stoc_map)[7] = "Omnivore (>50% plants)"
  colnames(stoc_map)[7] = "Omnivore (≤50% plants)"
  colnames(stoc_map)[8] = "Carnivore"
  
  p3 <- plotStochMaps(tree, maps=stoc_map,
                      tip_labels = FALSE,
                      tree_layout = "circular",
                      line_width=0.25,
                      color_by = "MAP",
                      colors = c("Herbivore"="#44AA99",
                                 "Omnivore (>50% plants)"="#999933",
                                 "Omnivore (≤50% plants)"="#cc6677",
                                 "Carnivore"="#882255")) +
    theme(legend.position.inside = c(0.6,0.81))
  #p3
  ggsave(paste0(dir_out, "run_", run, "_stoch_map.pdf"), p3, width = 8, height = 6)
}



