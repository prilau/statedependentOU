library(ggplot2)
library(phytools)
library(readr)
library(RevGadgets)
source("scripts/revgadgets_StochMap.R")

matchNodes = function(phy) {
  
  # get some useful info
  num_tips = length(phy$tip.label)
  num_nodes = phy$Nnode
  tip_indexes = 1:num_tips
  node_indexes = num_tips + num_nodes:1
  
  node_map = data.frame(R=1:(num_tips + num_nodes), Rev=NA, visits=0)
  current_node = phy$Nnode + 2
  k = 1
  t = 1
  
  while(TRUE) {
    
    if ( current_node <= num_tips ) {
      node_map$Rev[node_map$R == current_node] = t
      current_node = phy$edge[phy$edge[,2] == current_node,1]
      t = t + 1
    } else {
      
      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1
      
      if ( node_map$visits[node_map$R == current_node] == 1 ) {
        # go right
        current_node = phy$edge[phy$edge[,1] == current_node,2][2]
      } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
        # go left
        current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      } else if ( node_map$visits[node_map$R == current_node] == 3 ) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node = phy$edge[phy$edge[,2] == current_node,1]
        }
      }
    }
  }
  return(node_map[,1:2])
}
#tree <- read.tree("data/3_empirical/mammal_diet.tre")

tree <- read.tree("data/3_empirical/mammal_diet_perMY.tre")
index_to_rev <- matchNodes(tree)

log <- read_tsv("output/2_simulation/false_positive_stateless/sim_1/augch_blackPrior_run_1.trees")

log <- as.data.frame(log$char_hist)
write_tsv(log, "output/2_simulation/false_positive_stateless/sim_1/augch_blackPrior_run_1_modified.trees", col_names = FALSE)

simmap_to_ancStates <- function(input_path, output_path){
  simmaps <- read.simmap(input_path, format="phylip")
  df_rev <- data.frame()
  
  for (row_num in 1:length(simmaps)){
    simmap <- simmaps[[row_num]]
    
    # Iteration column
    df_rev[row_num, 1] <- row_num-1
    
    for (i in 1:(length(simmap$maps))){
      ape_node <- which(index_to_rev[,2]==i)
      ape_edge <- which(simmap$edge[,2]==ape_node)
      map <- simmap$maps[[ape_edge]]
      df_rev[row_num, i+1] <- names(map[length(map)])
    }
    df_rev[row_num, length(simmap$maps)+2] <- names(map[1])
  }
  header <- paste0("end_", as.character(1:(length(simmap$maps)+1)))
  colnames(df_rev) <- c("Iteration", header)
  write_tsv(df_rev, output_path)
}

simmap_to_ancStates("output/3_empirical/3.2a_sdOU/augch.trees", "output/3_empirical/3.2a_sdOU/anc_states.log")

# in rb
#tree <- readTrees("data/3_empirical/mammal_diet.tre")[1]
#anc_states = readAncestralStateTrace("output/3_empirical/3c_statedependentOU/anc_states_run2.log")
#anc_tree = ancestralStateTree(
#  tree=tree,
#  ancestral_state_trace_vector=anc_states,
#  include_start_states=false,
#  file="output/3_empirical/3c_statedependentOU/anc_states_run2.tre",
#  summary_statistic="MAP",
#  reconstruction="marginal",
#  burnin=0.0,
#  nStates=3,
#  site=1)


ase <- processAncStates("output/3_empirical/3c_statedependentOU/anc_states_run2.tre",
                        # Specify state labels.
                        # These numbers correspond to
                        # your input data file.
                        state_labels = c("0" = "Herbivore", "1" = "Omnivore", "2" = "Carnivore"))

# produce the plot object, showing MAP states at nodes.
# color corresponds to state, size to the state's posterior probability
p <- plotAncStatesMAP(t = ase,
                      tip_labels_offset = 0.00,
                      #tip_labels = FALSE,
                      node_color_as = "state",
                      node_color = c("Herbivore"="#44aa99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
                      node_size = c(0.1, 2.5),
                      tip_states_size = c(0.01, 0.01),
                      tip_states_shape = 1,
                      state_transparency = 0.70,
                      tree_layout = "rectangular",
                      tip_labels_size = 0.5) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.92,0.81))


ggsave(paste0("figures/3_empirical/3c_statedependentOU/anc_run2_wTipLabel.pdf"), p, width = 9, height = 20)


p <- plotAncStatesPie(t = ase,
                      pie_colors = c("Herbivore"="#44aa99", "Omnivore"="#ddcc77", "Carnivore"="#882255"),
                      tip_labels_size = 1,
                      tip_pies = TRUE,
                      tip_labels_offset = 0.02,
                      #tip_labels = FALSE,
                      state_transparency = 0.55,
                      tree_layout = "circular") +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.92,0.81))

ggsave(paste0("figures/3_empirical/3c_statedependentOU/pie_run2_wTipLabel.pdf"), p, width = 9, height = 12)


tree <- readTrees("data/3_empirical/mammal_foot.tre")
simmaps <- read.simmap("output/3_empirical/3.2a_sdOU/augch.trees", format="phylip")
stoc_map <- processStochMaps(tree,
                             simmap = simmaps,
                             states=c("0", "1", "2"))

colnames(stoc_map)[6] = "Digitigrady"
colnames(stoc_map)[7] = "Plantigrady"
colnames(stoc_map)[8] = "Unguligrady"

p <- plotStochMaps(tree, maps=stoc_map,
                   tip_labels = TRUE,
                   tree_layout = "circular",
                   line_width=0.25,
                   color_by = "MAP",
                   colors = c("Herbivore"="#44AA99",
                              "Omnivore"="#ddcc77",
                              "Carnivore"="#882255"),
                   tip_labels_size=1) +
  theme(legend.position.inside = c(0.92,0.81))
ggsave(paste0("figures/3_empirical/3c_statedependentOU/augch_run2_circ_wTipLabels.pdf"), p, width = 9, height = 12)

