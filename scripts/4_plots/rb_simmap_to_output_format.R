library(cowplot)
library(ggplot2)
library(gridExtra)
library(phytools)
library(readr)
library(RevGadgets)
source("scripts/4_plots/revgadgets_StochMap.R")


#tree <- read.tree("data/3_empirical/mammal_diet_perMY.tre")
#tree <- read.tree("data/2_simulation/mammal_diet_perMY_n500.tre")
#index_to_rev <- matchNodes(tree)
#
#dir_in = "output/2_simulation/convergence/augch_trees/"
#for (file in list.files(dir_in)){
#  path = paste0(dir_in, file)
#  log <- read_tsv(path)
#  log <- as.data.frame(log$char_hist)
#  write_tsv(log, file=path, col_names = FALSE)
#}
#
##simmap_to_ancStates <- function(input_path, output_path){
#simmaps <- read.simmap(input_path, format="phylip")
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

#simmap_to_ancStates("output/3_empirical/3.2a_sdOU/augch.trees", "output/3_empirical/3.2a_sdOU/anc_states.log")

# Do this in RevBayes!
# tree <- readTrees("data/2_simulation/mammal_diet_perMY_n500.tre")[1]
# 
# for (file in listFiles("output/2_simulation/convergence/augch_trees/")){
#   anc_states = readAncestralStateTrace(file)
#   anc_tree = ancestralStateTree(
#     tree=tree,
#     ancestral_state_trace_vector=anc_states,
#     include_start_states=false,
#     file=file + "ancStates.tre",
#     summary_statistic="MAP",
#     reconstruction="marginal",
#     burnin=0.1,
#     nStates=3,
#     site=1)
# }



#Plotting
tree <- readTrees("data/3_empirical/mammal_perMY_n2955.tre")
t <- read.tree("data/3_empirical/mammal_perMY_n2955.tre")

##############
# tip states #
##############
tip_states <- list()
tip_states[["Picky"]] = df_hpa %>% filter(picky == 1) %>% select(binomial_name) %>% unlist() %>% unname()
tip_states[["Not picky"]] = df_hpa %>% filter(picky == 0) %>% select(binomial_name) %>% unlist() %>% unname()
t <- tidytree::groupOTU(tree, tip_states)

colors <- c("Not picky"="#bbccee", "Picky"="#222255")

p <- ggtree(t, layout = "circular", color = "grey", linewidth = 0.1) +
  ggtree::geom_tippoint(ggplot2::aes(colour = group), size = 0.1, shape = 1, alpha = 0.4) +
  #geom_tiplab(aes(color=group), size = 0.4, offset = 0.5) +
  scale_color_manual(values=colors) +
  guides(color = guide_legend( 
    override.aes=list(size = 4, shape = 20, alpha = 0.8)))
ggsave(paste0("figures/3_empirical/ase_picky/tip_states.pdf"), p, width = 8, height = 6)



###################
# ancestral state #
###################
ase <- processAncStates("output/3_empirical/ase_picky/anc_states.tre",
                        state_labels=c("0"="Not picky", "1"="Picky"))
# produce the plot object, showing MAP states at nodes.
p <- plotAncStatesMAP(t = ase,
                      #tip_labels_offset = 0.5,
                      tip_labels = FALSE,
                      node_color_as = "state",
                      node_color = c("Not picky"="#bbccee", "Picky"="#222255"),
                      node_size = c(0.05, 1.2),
                      tip_states = FALSE,
                      #tip_states_size = c(0.01, 0.01),
                      #tip_states_shape = 1,
                      state_transparency = 0.6,
                      tree_layout = "circular",
                      #tip_labels_size = 0.5
                      tree_color = "#D3D3D3",
                      tree_linewidth = 0.1) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))

ggsave(paste0("figures/3_empirical/ase_picky/ase_map.pdf"), p, width = 8, height = 6)


p <- plotAncStatesPie(t = ase,
                      pie_colors = c("Not picky"="#bbccee", "Picky"="#222255"),
                      #tip_labels_size = 1,
                      tip_pies = FALSE,
                      #tip_labels_offset = 0.5,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.2,
                      tree_layout = "circular",
                      tip_labels = FALSE,
                      state_transparency = 0.6,
                      tree_color = "#D3D3D3",
                      tree_linewidth = 0.1) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))

ggsave(paste0("figures/3_empirical/ase_picky/ase_pie.pdf"), p, width = 8, height = 6)

#############
# carnivory #
#############
tip_states <- list()
tip_states[["Carnivore"]] = df_hpa %>% filter(carnivore == 1) %>% select(binomial_name) %>% unlist() %>% unname()
tip_states[["Not carnivore"]] = df_hpa %>% filter(carnivore == 0) %>% select(binomial_name) %>% unlist() %>% unname()
t <- tidytree::groupOTU(tree, tip_states)

colors <- c("Not carnivore"="#ffcccc", "Carnivore"="#663333")
p1 <- ggtree(t, layout = "circular", color = "grey", linewidth = 0.1) +
  ggtree::geom_tippoint(ggplot2::aes(colour = group), size = 0.1, shape = 1, alpha = 0.4) +
  #geom_tiplab(aes(color=group), size = 0.4, offset = 0.5) +
  scale_color_manual(values=colors) +
  theme(legend.position = "none")
  #guides(color = guide_legend( 
  #  override.aes=list(size = 4, shape = 20, alpha = 0.8)))
ggsave(paste0("figures/3_empirical/ase_carnivory/tip_states.pdf"), p, width = 8, height = 6)


ase <- processAncStates("output/3_empirical/ase_carnivory/anc_states.tre",
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        state_labels=c("0"="Not carnivore", "1"="Carnivore"))

# produce the plot object, showing MAP states at nodes.
# color corresponds to state, size to the state's posterior probability
p2 <- plotAncStatesMAP(t = ase,
                      #tip_labels_offset = 0.5,
                      node_color = c("Not carnivore"="#ffcccc", "Carnivore"="#663333"),
                      tip_labels = FALSE,
                      tip_states = FALSE,
                      node_color_as = "state",
                      node_size = c(0.05, 1.2),
                      #tip_states_size = c(0.01, 0.01),
                      #tip_states_shape = 1,
                      state_transparency = 0.6,
                      tree_layout = "circular",
                      #tip_labels_size = 0.5
                      tree_color = "#D3D3D3",
                      tree_linewidth = 0.1) +
  # modify legend location using ggplot2
  guides(color = "none")
legend_size <- get_legend(p2)
p2 <- p2 + guides(color = "none", size = "none")

ggsave(paste0("figures/3_empirical/ase_carnivory/ase_map.pdf"), p, width = 16, height = 9)

p3 <- plotAncStatesPie(t = ase,
                      pie_colors = c("Not carnivore"="#ffcccc", "Carnivore"="#663333"),
                      #tip_labels_size = 1,
                      tip_pies = FALSE,
                      #tip_labels_offset = 0.5,
                      node_pie_size = 0.4,
                      tip_pie_size = 0.2,
                      tree_layout = "circular",
                      tip_labels = FALSE,
                      state_transparency = 0.6,
                      tree_color = "#D3D3D3",
                      tree_linewidth = 0.1)

ggsave(paste0("figures/3_empirical/ase_carnivory/ase_pie.pdf"), p, width = 16, height = 9)

simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_carnivory/marginal_character.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1"))


colnames(stoc_map)[6] = "Not carnivore"
colnames(stoc_map)[7] = "Carnivore"
p4 <- plotStochMaps(tree, maps=stoc_map,
                   tip_labels = FALSE,
                   tree_layout = "circular",
                   line_width=0.1,
                   color_by = "MAP",
                   colors = c("Not carnivore"="#ffcccc", "Carnivore"="#663333")) +
  guides(color = "none")
ggsave(paste0("figures/3_empirical/ase_carnivory/stoch_map.pdf"), p4, width = 8, height = 6)


p_all <- arrangeGrob(p1, p2, p4, p3, nrow = 2)
ggsave("figures/3_empirical/ase_carnivory/ase_all.pdf", p_all, width = 16, height = 16)

##################
# stochastic map #
##################

# picky eater #
simmaps <- list()
simmaps[[1]] <- read.simmap("output/3_empirical/ase_picky/marginal_character.tree", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1"))

colnames(stoc_map)[6] = "Not picky"
colnames(stoc_map)[7] = "Picky"

p <- plotStochMaps(tree, maps=stoc_map,
                   tip_labels = FALSE,
                   tree_layout = "circular",
                   line_width=0.1,
                   color_by = "MAP",
                   colors = c("Not picky"="#BBCCEE",
                              "Picky"="#222255")) +
  theme(legend.position.inside = c(0.6,0.81))
ggsave(paste0("figures/3_empirical/ase_picky/stoch_map.pdf"), p, width = 8, height = 6)
