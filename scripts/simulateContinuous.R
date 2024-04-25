library(ape)
library(phytools)
library(geiger)
library(TESS)

simulateContinuous = function(tree, halflife, theta, sigma2) {

  ## Re-parameterization
  alpha <- log(2) / halflife
  #sigma2 <- stationaryvar * 2 * alpha

  preorder <- rev(postorder(tree))
  edges <- tree$edge
  root_node <- length(tree$tip.label) + 1
  state = tree$node.states[root_node]
  expected_mu <- rep(0, length(tree$node.states))
  expected_mu[root_node] <- theta[[state]]

  for (edge_index in preorder){
    sub_edges <- tree$maps[[edge_index]]
    parent_node <- edges[edge_index, 1]
    y <- expected_mu[parent_node]
    for (j in 1:length(sub_edges)) {
      state <- names(sub_edges[j])
      mu <- y * exp(-alpha[[state]] * sub_edges[[j]]) + theta[[state]] * (1 - exp(-alpha[[state]] * sub_edges[[j]]))
      v <- sigma2[[state]] / (2 * alpha[[state]]) * (1 - exp(-2 * alpha[[state]] * sub_edges[[j]]))
      y <- rnorm(n=1, mu, sqrt(v))    
    }
    desc_node <- edges[edge_index, 2]
    expected_mu[desc_node] <- y
  }
  
  cont_list <- list()
  for (i in 1:length(tree$tip.label)){
    tip <- tree$tip.label[i]
    cont_list[[tip]] <- expected_mu[tip]
  }
  
  return(cont_list)
}

stateless_halflife <- c(1, 1, 1)
sd_halflife <- c(0.5, 1, 5)
stateless_theta <- c(2, 2, 2)
sd_theta <- c(2, 4, 6)
stateless_sigma2 <- c(1, 1, 1)
sd_sigma2 <- c(0.05, 1, 2)

names(stateless_halflife) = 
  names(sd_halflife) = 
  names(stateless_theta) =
  names(sd_theta) =  
  names(stateless_sigma2) =  
  names(sd_sigma2) = 
  c("0", "1", "2")

cat("simulating continuous characters.\n")

sd_combination   = c("xxx", "xxt", "xsx", "axx", "asx", "axt", "xst", "ast")
reps       = 80

grid = expand.grid(combination=sd_combination, tree=1:reps,
                   stringsAsFactors=FALSE)
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_combo    = this_row[[1]]
  this_tree        = this_row[[2]]

  # read the history
  this_dir = paste0("data/2_simulation/2a_state_dependency/",this_combo, "/t", this_tree)
  load(paste0(this_dir, "/history.Rda"))
  
  if (isTRUE(grepl("a", this_combo))){
    halflife <- sd_halflife
  } else {
    halflife <- stateless_halflife
  }
  
  if (isTRUE(grepl("s", this_combo))){
    sigma2 <- sd_sigma2
  } else {
    sigma2 <- stateless_sigma2
  }
  
  if (isTRUE(grepl("t", this_combo))){
    theta <- sd_theta
  } else {
    theta <- stateless_theta
  }
  
  cont_states <- simulateContinuous(tree = history,
                                    halflife, theta, sigma2)
  
  write.nexus.data(cont_states, 
                   file = paste0(this_dir, "/continuous.nex"),
                   format="Continuous")
  setTxtProgressBar(bar, i / nrow(grid))
}
cat("\n")



