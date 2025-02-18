library(ape)
library(collapse)
library(ggtree)
library(phytools)
library(tidyverse)

tipMAP <- function(input_path, t){
  ntip <- length(t$tip.label)
  run <- read_tsv(input_path)
  map <- run %>% summarise(across(end_1:colnames(run[ntip+1]),
                                  collapse::fmode)) %>% as.numeric()
  return(map)
}

compareMAP <- function(map1, map2, rm.na=TRUE){
  if(isTRUE(rm.na)){
    map2 <- map2[!is.na(map1)]
    map1 <- map1[!is.na(map1)]
    map1 <- map1[!is.na(map2)]
    map2 <- map2[!is.na(map2)]
  }
  
  num_nodes = length(map1)
  p_same_state <- length(which(map1 == map2)) / num_nodes
  return(p_same_state)
}

markCompareMAP <- function(ref, map){
  mark_map <- map
  mark_map[which(map==ref)] = "match"
  mark_map[which(map!=ref)] = "mismatch"
  mark_map[which(is.na(ref))] = "NA"
  return(mark_map)
}

t <- read.tree("data/3_empirical/mammal_perMY_r500.tre")
dir_in = "output/3_empirical/sdOU_r500_missingStateModel/"
run6_sum <- tipMAP(paste0(dir_in, "augch_nstate_2_run_6.log"), t)


dir_in = "output/3_empirical/ase_r500_3StateOrdered/"
ase3_run4_sum <- tipMAP(paste0(dir_in, "states_run_4.log"), t)

ase3_run4_sum[which(ase3_run4_sum==1)] = NA
ase3_run4_sum[which(ase3_run4_sum==2)] = 1
compareMAP(ase3_run4_sum, run6_sum)
run6_mark <- markCompareMAP(ase3_run4_sum, run6_sum)


ase <- processAncStates("output/3_empirical/ase_r500_3StateOrdered/anc_states_run_4.tre",
                        state_labels=c("0"="match", "1"="mismatch", "2"="NA"))

ase@data[1:500,2] <- "1"
ase@data[1:500,4] <- "0"
ase@data[1:500,6] <- "0"


for (i in 1:500){
  ase@data[i,1] <- run6_mark[501-i]
}

p1 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = ,
                       node_color = c("match"="green", "mismatch"="red", "NA"="grey"),
                       node_size = c(0.0, 0.000001),
                       tip_states = TRUE,
                       tip_states_size = 1,
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  theme(legend.position.inside = c(0.6,0.81))
