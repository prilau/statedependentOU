library(TESS)

cat("simulating trees.\n")

# simulation parameters
num_tips   = c(100, 250, 500)
reps       = 5

grid = expand.grid(num_tips=num_tips, tree=1:reps, stringsAsFactors=FALSE)

# simulate the trees
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_tree       = this_row[[2]]
  
  # create the directories if necessary
  this_dir = paste0("data/n",this_num_tips, "/t", this_tree)
  if ( !dir.exists(this_dir) ) {
    dir.create(this_dir, recursive=TRUE, showWarnings=FALSE)
  }
  
  # simulate the tree
  tree = ladderize(tess.sim.taxa(1, this_num_tips, 10, 1, 0.5)[[1]])
  
  # rescale the tree
  tree$edge.length = tree$edge.length / max(branching.times(tree))
  
  # write the tree
  write.tree(tree, file=paste0(this_dir,"/tree.tre"))
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
  
}

cat("\n")
