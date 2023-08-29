library(TESS)
library(phytools)
source("scripts/readWriteCharacterData.R")

cat("simulating discrete characters.\n")

# simulation parameters
num_tips   = c(50, 100, 250)
num_traits = 1
reps       = 5

grid = expand.grid(num_tips=num_tips, num_traits=num_traits,
                   tree=1:reps, stringsAsFactors=FALSE)

# first, we need to compute the tree length of each replicate
tree_lengths = vector("list", length(num_tips))
for(i in 1:length(num_tips)) {
  this_dir = paste0("n_", num_tips[i])
  these_files = list.files(this_dir, pattern="tree.tre", recursive=TRUE, full.names = TRUE)
  these_trees = lapply(these_files, read.tree)
  these_tree_lengths = sapply(these_trees, function(tree) sum(tree$edge.length))
  tree_lengths[[i]] = these_tree_lengths
}

mean_tree_lengths = sapply(tree_lengths, mean)

# specify rates so that the expected number of changes is 5
rates = 5 / mean_tree_lengths
names(rates) = num_tips

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
  
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_num_traits = this_row[[2]]
  this_tree       = this_row[[3]]
  
  # read the tree
  this_dir = paste0("n_",this_num_tips,"/k_",this_num_traits,"/t_", this_tree)
  tree = read.tree(paste0(this_dir, "/tree.tre"))
  
  # get the rate
  this_rate = rates[as.character(this_num_tips)]
  
  # simulate the character
  history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
  
  # make sure at least 20% of the tips are in either state
  while (! (mean(history$states == "0") > 0.2 & (mean(history$states == "1") > 0.2) ) ) {
    history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
    num_rejections[as.character(this_num_tips)] = num_rejections[as.character(this_num_tips)] + 1
  }
  num_simulations[as.character(this_num_tips)] = num_simulations[as.character(this_num_tips)] + 1
  
  # make state-0 and state-1 trees
  # the branch length of state-i trees is the proportional
  # amount of time each branch spends in state i
  
  maps = history$mapped.edge[,c("0","1")]
  
  state_0_tree = tree
  state_0_tree$edge.length = maps[,1] / tree$edge.length
  
  state_1_tree = tree
  state_1_tree$edge.length = maps[,2] / tree$edge.length
  
  # save these trees
  write.tree(state_0_tree, file=paste0(this_dir,"/state_0.tre"))
  write.tree(state_1_tree, file=paste0(this_dir,"/state_1.tre"))
  
  # save the discrete trait as a nexus file
  writeCharacterData(t(t(history$states)), file=paste0(this_dir,"/disc_trait.nex"), type="Standard")
  
  # save the character history
  save(history, file=paste0(this_dir,"/history.Rda"))
  
  # write a pdf
  pdf(paste0(this_dir,"/history.pdf"))
  plot(history, col=colors)
  dev.off()
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
  
}

cat("\n")
