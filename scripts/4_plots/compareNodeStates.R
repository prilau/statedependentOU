library(ape)
library(collapse)
library(phytools)
library(tidyverse)

MAP <- function(input_path, t){
  ntip <- length(t$tip.label)
  run <- read_tsv(input_path)
  map <- run %>% summarise(across(colnames(run[ntip+1]):colnames(run[length(run)]),
                                   collapse::fmode)) %>% as.numeric()
  return(map)
}

compareMAP <- function(map1, map2){
  num_nodes = length(map1)
  p_same_state <- length(which(map1 == map2)) / num_nodes
  return(p_same_state)
}


dir_in = "output/3_empirical/aug_r500_3StateOrderedModel/"
aug3_run1_sum <- MAP(paste0(dir_in, "augch_run_1.log"), t)
aug3_run2_sum <- MAP(paste0(dir_in, "augch_run_2.log"), t)
aug3_run3_sum <- MAP(paste0(dir_in, "augch_run_3.log"), t)
aug3_run4_sum <- MAP(paste0(dir_in, "augch_run_4.log"), t)
aug3_run5_sum <- MAP(paste0(dir_in, "augch_run_5.log"), t)

dir_in = "output/3_empirical/ase_r500_3StateOrdered/"
ase3_run1_sum <- MAP(paste0(dir_in, "states_run_1.log"), t)
ase3_run4_sum <- MAP(paste0(dir_in, "states_run_4.log"), t)
ase3_run5_sum <- MAP(paste0(dir_in, "states_run_5.log"), t)

compareMAP(ase3_run4_sum, aug3_run5_sum)
