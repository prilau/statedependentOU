library(TESS)
library(phytools)
source("scripts/readWriteCharacterData.R")

# 3 STATES
cat("simulating discrete characters.\n")

# simulation parameters
num_sim       = 1800

#grid = expand.grid(models=models, tree=1:reps, stringsAsFactors=FALSE)

# first, we need to compute the tree length of each replicate
#tree_lengths = vector("list", length(num_tips))
#for(i in 1:length(num_tips)) {
#  this_dir = paste0("../sdOU_local/IRT3/data/simulation/")
#  these_files = list.files(this_dir, pattern="tree.tre", recursive=TRUE, full.names = TRUE)
#  these_trees = lapply(these_files, read.tree)
#  these_tree_lengths = sapply(these_trees, function(tree) sum(tree$edge.length))
#  tree_lengths[[i]] = these_tree_lengths
#}

tree <- read.tree("data/2_simulation/mammal_diet_height1_n500.tre")
tree_length = sum(tree$edge.length)
#mean_tree_lengths = sapply(tree_lengths, mean)

# specify rates so that the expected number of changes is 200
this_rate = 100 / tree_length
#names(rates) = models

# specify the Mk2 rate matrix
Q = matrix(0.5, 3, 3)
diag(Q) = -1
rownames(Q) = colnames(Q) = 1:3 

# simulate the discrete characters
# track the number of rejected simulations based on proportional
# representation
colors = c("1"="#44aa99", "2"="#ddcc77", "3"="#882255")
#num_rejections = numeric(length(num_tips))
#num_simulations = numeric(length(num_tips))
#names(num_rejections) = names(num_simulations) = num_tips

bar = txtProgressBar(style=3, width=40)
for(i in 201:(200+num_sim)) {
  
  # simulate the character
  history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
  
  # make sure at least 10% of the tips are in either state
  while (!(mean(history$states == "1") > 0.10 & mean(history$states == "2") > 0.10 & mean(history$states == "3") > 0.10) ) {
    history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
  }
  
  maps = history$mapped.edge[,c("1","2","3")]

  this_sub_dir = paste0("data/2_simulation/triState/sim_", i)
  if ( !dir.exists(this_sub_dir) ) {
    dir.create(this_sub_dir, recursive=TRUE, showWarnings=FALSE)
  }

  # save the discrete trait as a nexus file
  writeCharacterData(t(t(history$states)), file=paste0(this_sub_dir, "/discrete.nex"), type="Standard")
  
  # save the character history
  save(history, file=paste0(this_sub_dir, "/history.Rda"))
  
  # write a pdf
  pdf(paste0(this_sub_dir, "/history.pdf"))
  plot(history, col=colors)
  dev.off()
  
  # increment the progress bar
  setTxtProgressBar(bar, i / num_sim)
  
}

cat("\n")











# BINARY STATE
cat("simulating discrete characters.\n")

# simulation parameters
num_sim = 1000

#grid = expand.grid(models=models, tree=1:reps, stringsAsFactors=FALSE)


# first, we need to compute the tree length of each replicate
#for(i in 1:length(num_tips)) {
tree <- read.tree("data/2_simulation/mammal_diet_height1_n500.tre")
tree_length = sum(tree$edge.length)
#}

# specify rates so that the expected number of changes is 50
rate = 100 / tree_length
#names(rates) = models

# specify the Mk2 rate matrix
#Q = matrix(0.5, 3, 3)
#diag(Q) = -1
#rownames(Q) = colnames(Q) = 1:3 - 1
Q = matrix(1, 2, 2)
diag(Q) = -1
rownames(Q) = colnames(Q) = 1:2


# simulate the discrete characters
# track the number of rejected simulations based on proportional
# representation
#colors = c("1"="#44aa99", "2"="#ddcc77", "2"="#882255")
colors = c("1"="#44aa99", "2"="#882255")
#num_rejections = numeric(length(num_tips))
#num_simulations = numeric(length(num_tips))
#names(num_rejections) = names(num_simulations) = num_tips

bar = txtProgressBar(style=3, width=40)
for(i in 1:num_sim) {
  
  # read the tree
  #this_dir = paste0("../sdOU_local/IRT3/data/simulation/", this_model, "/t", this_tree)
  #tree = read.tree(paste0(this_dir, "/tree.tre"))
  
  # get the rate
  #this_rate = rates[as.character(this_model)]
  #this_rate = rate
  
  # simulate the character
  history = sim.history(tree, rate * Q, nsim=1, message=FALSE)
  
  # make sure at least 20% of the tips are in either state
  #while (! (mean(history$states == "1") > 0.125 & mean(history$states == "2") > 0.125 & mean(history$states == "2") > 0.125) ) {
  #  history = sim.history(tree, rate * Q, nsim=1, message=FALSE)
  #  #num_rejections[as.character(num_tips)] = num_rejections[as.character(num_tips)] + 1
  #}
  while (! (mean(history$states == "1") > 0.2 & mean(history$states == "2") > 0.2 ) ) {
    history = sim.history(tree, rate * Q, nsim=1, message=FALSE)
    #num_rejections[as.character(num_tips)] = num_rejections[as.character(num_tips)] + 1
  }
  #num_simulations[as.character(num_tips)] = num_simulations[as.character(num_tips)] + 1
  
  # the branch length of state-i trees is the proportional
  # amount of time each branch spends in state i
  
  #maps = history$mapped.edge[,c("1","2","2")]
  maps = history$mapped.edge[,c("1","2")]
  
  # save these trees
  this_sub_dir = paste0("data/2_simulation/false_positive/sim_", i)
  if ( !dir.exists(this_sub_dir) ) {
    dir.create(this_sub_dir, recursive=TRUE, showWarnings=FALSE)
  }
  
  # save the discrete trait as a nexus file
  writeCharacterData(t(t(history$states)), file=paste0(this_sub_dir, "/discrete.nex"), type="Standard")
  
  # save the character history
  save(history, file=paste0(this_sub_dir, "/history.Rda"))
  
  # write a pdf
  #pdf(paste0(this_sub_dir, "/history.pdf"))
  #plot(history, col=colors)
  #dev.off()
  
  # increment the progress bar
  setTxtProgressBar(bar, i / num_sim)
  
} 
cat("\n")
