library(ape)
library(phytools)
library(ggplot2)
library(tidyverse)

dir_in = "output/1_validation/aug_tipState_hidden/"

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
state_count_per_tip <- tibble(t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0)

bar = txtProgressBar(style=3, width=40)
for (i in 1:length(simmaps)){
  for (j in 1:length(tip_edges)){ # add the tip state (0 or 1) for each tip
    num_states <- length(simmaps[[i]]$maps[[tip_edges[j]]])
    state <- names(simmaps[[i]]$maps[[tip_edges[j]]][num_states]) # obtain last (aka most recent) state on the tip branch
    state_count_per_tip[i,j] <- as.numeric(state)
  }
  setTxtProgressBar(bar, i / length(simmaps))
}

# t1
c(length(which(state_count_per_tip$t1 == 0)) > 0,
  length(which(state_count_per_tip$t1 == 1)) > 0,
  length(which(state_count_per_tip$t1 == 2)) > 0,
  length(which(state_count_per_tip$t1 == 3)) > 0)
c(length(which(state_count_per_tip$t1 == 0)) / length(which(state_count_per_tip$t1 == 1)),
  length(which(state_count_per_tip$t1 == 2)) / length(which(state_count_per_tip$t1 == 3)),
  length(which(state_count_per_tip$t1 == 0)) / length(which(state_count_per_tip$t1 == 2)),
  length(which(state_count_per_tip$t1 == 1)) / length(which(state_count_per_tip$t1 == 3)))

# t2
c(length(which(state_count_per_tip$t2 == 0)) > 0,
  length(which(state_count_per_tip$t2 == 1)) > 0,
  length(which(state_count_per_tip$t2 == 2)) > 0,
  length(which(state_count_per_tip$t2 == 3)) > 0)
c(length(which(state_count_per_tip$t2 == 0)) / length(which(state_count_per_tip$t2 == 1)),
  length(which(state_count_per_tip$t2 == 2)) / length(which(state_count_per_tip$t2 == 3)),
  length(which(state_count_per_tip$t2 == 0)) / length(which(state_count_per_tip$t2 == 2)),
  length(which(state_count_per_tip$t2 == 1)) / length(which(state_count_per_tip$t2 == 3)))

# t3
c(length(which(state_count_per_tip$t3 == 0)) > 0,
  length(which(state_count_per_tip$t3 == 1)) > 0,
  length(which(state_count_per_tip$t3 == 2)) > 0,
  length(which(state_count_per_tip$t3 == 3)) > 0)
c(length(which(state_count_per_tip$t3 == 0)) / length(which(state_count_per_tip$t3 == 1)),
  length(which(state_count_per_tip$t3 == 2)) / length(which(state_count_per_tip$t3 == 3)),
  length(which(state_count_per_tip$t3 == 0)) / length(which(state_count_per_tip$t3 == 2)),
  length(which(state_count_per_tip$t3 == 1)) / length(which(state_count_per_tip$t3 == 3)))


# t4
c(length(which(state_count_per_tip$t4 == 0)) > 0,
  length(which(state_count_per_tip$t4 == 1)) > 0,
  length(which(state_count_per_tip$t4 == 2)) > 0,
  length(which(state_count_per_tip$t4 == 3)) > 0)
c(length(which(state_count_per_tip$t4 == 0)) / length(which(state_count_per_tip$t4 == 1)),
  length(which(state_count_per_tip$t4 == 2)) / length(which(state_count_per_tip$t4 == 3)),
  length(which(state_count_per_tip$t4 == 0)) / length(which(state_count_per_tip$t4 == 2)),
  length(which(state_count_per_tip$t4 == 1)) / length(which(state_count_per_tip$t4 == 3)))

# t5
c(length(which(state_count_per_tip$t5 == 0)) > 0,
  length(which(state_count_per_tip$t5 == 1)) > 0,
  length(which(state_count_per_tip$t5 == 2)) > 0,
  length(which(state_count_per_tip$t5 == 3)) > 0)
c(length(which(state_count_per_tip$t5 == 0)) / length(which(state_count_per_tip$t5 == 1)),
  length(which(state_count_per_tip$t5 == 2)) / length(which(state_count_per_tip$t5 == 3)),
  length(which(state_count_per_tip$t5 == 0)) / length(which(state_count_per_tip$t5 == 2)),
  length(which(state_count_per_tip$t5 == 1)) / length(which(state_count_per_tip$t5 == 3)))

# t6
c(length(which(state_count_per_tip$t6 == 0)) > 0,
  length(which(state_count_per_tip$t6 == 1)) > 0,
  length(which(state_count_per_tip$t6 == 2)) > 0,
  length(which(state_count_per_tip$t6 == 3)) > 0)
c(length(which(state_count_per_tip$t6 == 0)) / length(which(state_count_per_tip$t6 == 1)),
  length(which(state_count_per_tip$t6 == 2)) / length(which(state_count_per_tip$t6 == 3)),
  length(which(state_count_per_tip$t6 == 0)) / length(which(state_count_per_tip$t6 == 2)),
  length(which(state_count_per_tip$t6 == 1)) / length(which(state_count_per_tip$t6 == 3)))

# t7
#length(which(state_count_per_tip$t7 == 0))
#length(which(state_count_per_tip$t7 == 1))
#length(which(state_count_per_tip$t7 == 2))
#length(which(state_count_per_tip$t7 == 3))
#
## t8
#length(which(state_count_per_tip$t8 == 0))
#length(which(state_count_per_tip$t8 == 1))
#length(which(state_count_per_tip$t8 == 2))
#length(which(state_count_per_tip$t8 == 3))
#
## t9
#length(which(state_count_per_tip$t9 == 0))
#length(which(state_count_per_tip$t9 == 1))
#length(which(state_count_per_tip$t9 == 2))
#length(which(state_count_per_tip$t9 == 3))
#
## t10
#length(which(state_count_per_tip$t10 == 0))
#length(which(state_count_per_tip$t10 == 1))
#length(which(state_count_per_tip$t10 == 2))
#length(which(state_count_per_tip$t10 == 3))







# To subset tree to a smaller size

tree <- read.nexus("data/3_empirical/raw/4705sp_mammal-time_Álvarez-Carretero_2022.tree")

n=7
m=3
tree_small <- keep.tip(tree, sample(tree$tip.label, n))
tree_small$tip.label <- paste0("t", 1:n)
max(node.depth.edgelength(tree_small))
plot(tree_small)
write.tree(tree_small, paste0("data/1_validation/dummy/dummy_r", n, "_", m, ".tre"))

