library(TESS)
library(phytools)
source("scripts/readWriteCharacterData.R")

cat("simulating discrete characters.\n")

# simulation parameters
num_tips = 500
exp_change = c(5, 10, 20)
num_traits = 4

grid = expand.grid(exp_change=exp_change,
                   traits=1:num_traits, stringsAsFactors=FALSE)

# first, we need to compute the tree length
this_dir = paste0("data/simulation_3_sameTree_diffRate/n500/t1")
tree <- read.tree(paste0(this_dir, "/tree.tre"))
tree_length = sum(tree$edge.length)


# specify the Mk2 rate matrix
Q = matrix(1, 2, 2)
diag(Q) = -1
rownames(Q) = colnames(Q) = 1:2 - 1

# simulate the discrete characters
# track the number of rejected simulations based on proportional
# representation
colors = c("0"="blue","1"="red")
num_rejections = numeric(length(num_tips))
num_simulations = numeric(length(num_tips))
names(num_rejections) = names(num_simulations) = num_tips

bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  
  this_row        = grid[i,]
  this_exp_change = this_row[[1]]
  this_num_traits = this_row[[2]]
  
  
  # read the tree
  this_dir = paste0("data/simulation_3_sameTree_diffRate/n500/t1")
  
  tree = read.tree(paste0(this_dir, "/tree.tre"))
  
  # get the rate
  # specify rates so that the expected number of changes is 10
  this_rate = this_exp_change / tree_length
  
  # simulate the character
  history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
  
  # make sure at least 20% of the tips are in either state
  while (! (mean(history$states == "0") > 0.2 & (mean(history$states == "1") > 0.2) ) ) {
    history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
    num_rejections[as.character(500)] = num_rejections[as.character(500)] + 1
  }
  num_simulations[as.character(500)] = num_simulations[as.character(500)] + 1
  
  # make state-0 and state-1 trees
  # the branch length of state-i trees is the proportional
  # amount of time each branch spends in state i
  
  maps = history$mapped.edge[,c("0","1")]
  
  state_0_tree = tree
  state_0_tree$edge.length = maps[,1] / tree$edge.length
  
  state_1_tree = tree
  state_1_tree$edge.length = maps[,2] / tree$edge.length
  
  # save these trees
  this_sub_dir = paste0("data/simulation_3_sameTree_diffRate/n500/t1", "/r",this_exp_change, "/d", this_num_traits)
  if ( !dir.exists(this_sub_dir) ) {
    dir.create(this_sub_dir, recursive=TRUE, showWarnings=FALSE)
  }
  
  write.tree(state_0_tree, file=paste0(this_sub_dir, "/n500t1r", this_exp_change, "d", this_num_traits, "_State0.tre"))
  write.tree(state_1_tree, file=paste0(this_sub_dir, "/n500t1r", this_exp_change, "d", this_num_traits, "_State1.tre"))
  
  # save the discrete trait as a nexus file
  writeCharacterData(t(t(history$states)), file=paste0(this_sub_dir, "/n500t1r", this_exp_change, "d", this_num_traits, "_Discrete.nex"), type="Standard")
  
  # save the character history
  save(history, file=paste0(this_sub_dir, "/n500t1r", this_exp_change, "d", this_num_traits, "_History.Rda"))
  
  # write a pdf
  pdf(paste0(this_sub_dir, "/n500t1r", this_exp_change, "d", this_num_traits, "_History.pdf"))
  plot(history, col=colors)
  dev.off()
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
  
}

cat("\n")
