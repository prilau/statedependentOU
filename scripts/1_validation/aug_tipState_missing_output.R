library(ape)
library(phytools)
library(ggplot2)
library(tidyverse)

dir_in = "output/1_validation/aug_tipState_missing"

# Uncomment the block below if outputs are freshly from RevBayes

#for (file in list.files(dir_in)){
#  path = paste0(dir_in, file)
#  log <- read_tsv(path, show_col_types = FALSE)
#  log <- as.data.frame(log$char_hist)
#  write_tsv(log, file=path, col_names = FALSE)
#}

# validate each tip expected sample state equally
# validate each tree expected sample state == 3
path <- paste0(dir_in, list.files(dir_in))
simmaps <- read.simmap(path, format="phylip")

# obtain indices of tip branches
tip_edges <- c(which.edge(simmaps[[1]], "t1"),
               which.edge(simmaps[[1]], "t2"),
               which.edge(simmaps[[1]], "t3"),
               which.edge(simmaps[[1]], "t4"),
               which.edge(simmaps[[1]], "t5"),
               which.edge(simmaps[[1]], "t6"))

# initiate number of state 1 sampled per tip
state_count_per_tip <- rep(0,6)

bar = txtProgressBar(style=3, width=40)
for (i in 1:length(simmaps)){
  for (j in 1:length(tip_edges)){ # add the tip state (0 or 1) for each tip
    num_states <- length(simmaps[[i]]$maps[[tip_edges[j]]])
    state <- names(simmaps[[i]]$maps[[tip_edges[j]]][num_states]) # obtain last (aka most recent) state on the tip branch
    state_count_per_tip[j] <- state_count_per_tip[j] + as.numeric(state)
  }
  setTxtProgressBar(bar, i / length(simmaps))
}

p_state1_per_tip = state_count_per_tip / length(simmaps) # expected == 0.5 for each tip
p_state1_per_tree <- sum(state_count_per_tip) / length(simmaps) # expected == 3