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

t <- read.tree("data/1_validation/artiodactyla/artiodactyla.tree")
dir_in = "output/2_simulation/missing_state/"
dir_in2 = "data/2_simulation/missing_state/"

num_match = 0
for (i in 1:16){
  subdir <- paste0(dir_in, "sim_", i)
  sum <- tipMAP(paste0(subdir, "/augch_nstate_2_run_2.log"), t)
  
  subdir2 <- paste0(dir_in2, "sim_", i) 
  sum2 <- read_tsv(paste0(subdir2, "/nstate_2_true_history.log")) %>%
    select(2:44) %>% 
    pivot_longer(cols=1:43, names_to = "node", values_to = "state") %>% 
    mutate(state = as.numeric(substring(state, 2,2))-1) %>% 
    select(state) %>% unlist %>% unname
  num_match = num_match + sum(sum==sum2)
}
perc_match = num_match / (43 * 16)
