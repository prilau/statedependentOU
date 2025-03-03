library(phytools)
library(ggtree)
library(gridExtra)
source("scripts/5_miscellaneous/functions.R")

tree <- read.tree("data/1_validation/artiodactyla/artiodactyla.tree")
t <- readTrees("data/1_validation/artiodactyla/artiodactyla.tree")
dir_in = "output/2_simulation/missing_state/"
files <- list.files(dir_in, recursive = TRUE)[grepl("trees", list.files(dir_in, recursive = TRUE))]



sim <- 1:16

p_true <- list()
p_inf <- list()

for ( i in 1:length(sim) ){
  simmaps <- list()
  load(paste0("data/2_simulation/missing_state/sim_",
              sim[i], "/nstate_2_history.Rda"))
  simmaps[[1]] <- history2
  stoc_map <- processStochMaps(t,
                               simmap = simmaps,
                               states=c("1", "2", "3", "4"))
  colnames(stoc_map)[6] = "0 (large optimum)"
  colnames(stoc_map)[7] = "1 (small optimum)"
  
  p_true[[i]] <- plotStochMaps(t, maps=stoc_map,
                               tip_labels = FALSE,
                               #tree_layout = "circular",
                               line_width=0.5,
                               color_by = "MAP",
                               colors = c("0 (large optimum)"="#364B9a",
                                          "1 (small optimum)"="#c2e4ef")) +
    theme(legend.position = "none") +
    scale_x_reverse()
  
  ase <- processAncStates(paste0("output/2_simulation/missing_state/sim_", sim[i], "/anc_states_nstate_2_run_2.log"),
                          state_labels=c("0"="0 (large optimum)",
                                         "1"="1 (small optimum)",
                                         "2"="2",
                                         "3"="3"))
  p_inf[[i]] <- plotAncStatesPie(t = ase,
                     pie_colors = c("0 (large optimum)"="#364B9a", "1 (small optimum)"="#c2e4ef", "2"="#bbccbb", "3"="#000000"),
                     #tip_labels_size = 1,
                     tip_pies = TRUE,
                     #tip_labels_offset = 0.5,
                     node_pie_size = 3,
                     tip_pie_size = 3,
                     #tree_layout = "circular",
                     tip_labels = FALSE,
                     state_transparency = 0.85,
                     tree_color = "#aaaaaa",
                     tree_linewidth = 0.5) +
    # modify legend location using ggplot2
    theme(legend.position = "none")
}



row1 <- cowplot::plot_grid(p_inf[[1]], p_true[[1]],
                           p_inf[[2]], p_true[[2]],
                           p_inf[[3]], p_true[[3]],
                           p_inf[[4]], p_true[[4]],
                           p_inf[[5]], p_true[[5]],
                           ncol=10, labels=c("#1", "",
                                             "#2", "",
                                             "#3", "",
                                             "#4", "",
                                             "#5", ""))
row2 <- cowplot::plot_grid(p_inf[[6]], p_true[[6]],
                           p_inf[[7]], p_true[[7]],
                           p_inf[[8]], p_true[[8]],
                           p_inf[[9]], p_true[[9]],
                           p_inf[[10]], p_true[[10]],
                           ncol=10, labels=c("#6", "",
                                             "#7", "",
                                             "#8", "",
                                             "#9", "",
                                             "#10", ""))

p_all <- cowplot::plot_grid(row1, row2, ncol=1)

ggsave(file="figures/2_simulation/missing_state/compareTipStates.pdf", p_all, width = 8, height = 5.5)
