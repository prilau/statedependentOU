library(phytools)
library(RevGadgets)
library(ggplot2)
source("scripts/4_plots/revgadgets_StochMap.R")


#tree <- read.tree("data/3_empirical/dummy_r6.tre")
#index_to_rev <- matchNodes(tree)
#
#dir_in = "output/3_empirical/aug_tipMissingDiscrete/"
#
#for (file in list.files(dir_in)){
#  path = paste0(dir_in, file)
#  log <- read_tsv(path)
#  log <- as.data.frame(log$char_hist)
#  write_tsv(log, file=path, col_names = FALSE)
#}
#
##simmap_to_ancStates <- function(input_path, output_path){
#input_path <- paste0(dir_in, list.files(dir_in))
#simmaps <- read.simmap(input_path, format="phylip")
#
#df_rev <- data.frame()
#
#for (row_num in 1:length(simmaps)){
#  simmap <- simmaps[[row_num]]
#  
#  # Iteration column
#  df_rev[row_num, 1] <- row_num-1
#  
#  for (i in 1:(length(simmap$maps))){
#    ape_node <- which(index_to_rev[,2]==i)
#    ape_edge <- which(simmap$edge[,2]==ape_node)
#    map <- simmap$maps[[ape_edge]]
#    df_rev[row_num, i+1] <- names(map[length(map)])
#  }
#  df_rev[row_num, length(simmap$maps)+2] <- names(map[1])
#}
#header <- paste0("end_", as.character(1:(length(simmap$maps)+1)))
#colnames(df_rev) <- c("Iteration", header)
#write_tsv(df_rev, output_path)
#
#simmap_to_ancStates("output/3_empirical/3.2a_sdOU/augch.trees", "output/3_empirical/3.2a_sdOU/anc_states.log")


#Do this in RevBayes!
#tree <- readTrees("data/3_empirical/dummy_r6.tre")[1]
#
#for (file in listFiles("output/3_empirical/aug_tipMissingDiscrete/augch/")){
#  anc_states = readAncestralStateTrace(file)
#  anc_tree = ancestralStateTree(
#    tree=tree,
#    ancestral_state_trace_vector=anc_states,
#    include_start_states=false,
#    file=file + "ancStates.tre",
#    summary_statistic="MAP",
#    reconstruction="marginal",
#    burnin=0.1,
#    nStates=3,
#    site=1)
#}



#Plotting
tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")
t <- read.tree("data/3_empirical/mammal_perMY_r500.tre")


#################
# 3-state model #
#################
ase <- processAncStates("output/3_empirical/ase_r500_3State/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivore", "1"="Omnivore", "2"="Carnivore"))
# produce the plot object, showing MAP states at nodes.
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
p1
ggsave(paste0("figures/3_empirical/ase_r500_3State/ase_map.pdf"), p1, width = 8, height = 6)


p2 <- plotAncStatesPie(t = ase,
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
p2
ggsave(paste0("figures/3_empirical/ase_r500_3State/ase_pie.pdf"), p2, width = 8, height = 6)


# p3
simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_3State/marginal_character_run_1.tree", format="phylip")
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
p3
ggsave(paste0("figures/3_empirical/ase_r500_3State/stoch_map.pdf"), p3, width = 8, height = 6)

#################
# 4-state model #
#################
ase <- processAncStates("output/3_empirical/ase_r500_4State/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivore", "1"="Plant-dominant omnivore",
                                       "2"="Non-plant-dominant omnivore", "3"="Carnivore"))
# produce the plot object, showing MAP states at nodes.
p4 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivore"                   = "#44AA99",
                                      "Plant-dominant omnivore"     = "#999933",
                                      "Non-plant-dominant omnivore" = "#CC6677",
                                      "Carnivore"                   = "#882255"),
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
                       pie_colors = c("Herbivore"="#44AA99",
                                      "Plant-dominant omnivore"="#999933",
                                      "Non-plant-dominant omnivore"="#CC6677",
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
#p5
ggsave(paste0("figures/3_empirical/ase_r500_4State/ase_pie.pdf"), p5, width = 8, height = 6)


simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_4State/marginal_character_run_1.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2", "3"))

colnames(stoc_map)[6] = "Herbivore"
colnames(stoc_map)[7] = "Plant-dominant omnivore"
colnames(stoc_map)[8] = "Non-plant-dominant omnivore"
colnames(stoc_map)[9] = "Carnivore"

p6 <- plotStochMaps(tree, maps=stoc_map,
                    tip_labels = FALSE,
                    tree_layout = "circular",
                    line_width=0.25,
                    color_by = "MAP",
                    colors = c("Herbivore"="#44AA99",
                               "Plant-dominant omnivore"="#999933",
                               "Non-plant-dominant omnivore"="#CC6677",
                               "Carnivore"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
#p6
ggsave(paste0("figures/3_empirical/ase_r500_4State/stoch_map.pdf"), p6, width = 8, height = 6)

######################
# hidden-state model #
######################
ase <- processAncStates("output/3_empirical/ase_r500_hiddenState/anc_states_run_2.tre",
                        state_labels=c("0"="Herbivore_0", "1"="Omnivore_0",
                                       "2"="Carnivore_0", "3"="Herbivore_1",
                                       "4"="Omnivore_1", "5"="Carnivore_1"))
# produce the plot object, showing MAP states at nodes.
p7 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Herbivore_0"="#44AA99", "Omnivore_0"="#DDCC77", "Carnivore_0"="#CC6677",
                                      "Herbivore_1"="#009988", "Omnivore_1"="#997700", "Carnivore_1"="#994455"),
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
p7
ggsave(paste0("figures/3_empirical/ase_r500_hiddenState/ase_map.pdf"), p, width = 8, height = 6)


p8 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Herbivore_0"="#44AA99", "Omnivore_0"="#DDCC77", "Carnivore_0"="#CC6677",
                                      "Herbivore_1"="#009988", "Omnivore_1"="#997700", "Carnivore_1"="#994455"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 2,
                       tip_pie_size = 1,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
p8
ggsave(paste0("figures/3_empirical/ase_r500_hiddenState/ase_pie.pdf"), p8, width = 8, height = 6)


simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_r500_hiddenState/marginal_character.tree", format="phylip")
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
                               "Omnivore"="#999933",
                               "Non-plant-dominant omnivore"="#CC6677",
                               "Carnivore"="#882255")) +
  theme(legend.position.inside = c(0.6,0.81))
p9
ggsave(paste0("figures/3_empirical/ase_r500_hiddenState/stoch_map.pdf"), p, width = 8, height = 6)

#p_all <- arrangeGrob(p1, p2, p4, p3, nrow = 2)
#ggsave("figures/3_empirical/ase_triState/ase_all.pdf", p_all, width = 16, height = 16)

