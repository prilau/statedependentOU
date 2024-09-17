library(ape)
library(phytools)
library(geiger)
library(TESS)
library(tidyverse)

drawAlpha <- function(state_dependent){
  if(state_dependent == T){
    alpha <- sort(rlnorm(n=3, meanlog=-2.5, sdlog=1))
  } else {
    alpha <- rep(rlnorm(n=1, meanlog=-3, sdlog=1), 3)
  }
  names(alpha) = c("1", "0", "2")
  return(alpha)
} 

drawSigma2 <- function(state_dependent){
  if(state_dependent == T){
    sigma2 <- rlnorm(n=3, meanlog=-1, sdlog=1)
  } else {
    sigma2 <- rep(rlnorm(n=1, meanlog=-0.5, sdlog=1), 3)
  }
  names(sigma2) = c("0", "1", "2")
  return(sigma2)
}

drawTheta <- function(state_dependent){
  if(state_dependent == T){
    theta <- c(rnorm(n=1, mean=-5, sd=1),
               rnorm(n=1, mean=0, sd=1),
               rnorm(n=1, mean=5, sd=1))
  } else {
    theta <- rep(rnorm(n=1, mean=0, sd=1), 3)
  }
  names(theta) = c("0", "1", "2")
  return(theta)
}

#drawRootState <- function(alphaRoot, sigma2Root, thetaRoot){
#  rootState <- rnorm(n=1, mean=thetaRoot, sd = sigma2Root/(2* alphaRoot))
#  return(rootState)
#}

simulateContinuous = function(tree, stateDependencies) {

  ## Re-parameterization
  #alpha <- log(2) / halflife
  #sigma2 <- stationaryvar * 2 * alpha

  
  alpha <- drawAlpha(state_dependent = T)
  sigma2 <- drawSigma2(state_dependent = T)
  theta <- drawTheta(state_dependent = T)
  
  preorder <- rev(postorder(tree))
  edges <- tree$edge
  root_node <- length(tree$tip.label) + 1
  state = tree$node.states[root_node]
  mu_at_nodes <- rep(0, length(tree$node.states))
  mu_at_nodes[root_node] <- theta[[state]]

  
  
  #draw root state
  #alpha_root <- drawAlpha(state_dependent = stateDependencies[1])[[state]]
  #sigma2_root <- drawAlpha(state_dependent = stateDependencies[2])[[state]]
  #theta_root <- drawTheta(state_dependent = stateDependencies[3])[[state]]
  #mu_at_nodes[root_node] <- drawRootState(alpha_root, sigma2_root, theta_root)
  
  #if (rootState == "optimum"){
  #  mu_at_nodes[root_node] <- theta[[state]]
  #} else if (rootState == "equilibrium") {
  #  mu_at_nodes[root_node] <- rnorm(1, mean = theta[[state]],
  #                                  sd = sqrt(sigma2[[state]]/(2*alpha[[state]])))
  #} else if (rootState == "parameter"){
  #  mu_at_nodes[root_node] <- runif(1, min = min(theta)/2, max = max(theta)*2)
  #}
  
  for (edge_index in preorder){
    sub_edges <- tree$maps[[edge_index]]
    parent_node <- edges[edge_index, 1]
    y <- mu_at_nodes[parent_node]
    for (j in 1:length(sub_edges)) {
      #alpha <- drawAlpha(state_dependent = stateDependencies[1])
      #sigma2 <- drawAlpha(state_dependent = stateDependencies[2])
      #theta <- drawAlpha(state_dependent = stateDependencies[3])
      
      state <- names(sub_edges[j])
      mu <- y * exp(-alpha[[state]] * sub_edges[[j]]) + theta[[state]] * (1 - exp(-alpha[[state]] * sub_edges[[j]]))
      v <- sigma2[[state]] / (2 * alpha[[state]]) * (1 - exp(-2 * alpha[[state]] * sub_edges[[j]]))
      y <- rnorm(n=1, mu, sqrt(v))    
    }
    desc_node <- edges[edge_index, 2]
    mu_at_nodes[desc_node] <- y
  }
  
  cont_list <- list()
  for (i in 1:length(tree$tip.label)){
    tip <- tree$tip.label[i]
    cont_list[[tip]] <- mu_at_nodes[i]
  }
  
  return(list(cont_list, alpha, sigma2, theta))
}

par <- tibble(alpha_0 = 0, alpha_1 = 0, alpha_2 = 0,
              sigma2_0 = 0, sigma2_1 = 0, sigma2_2 = 0,
              theta_0 = 0, theta_1 = 0, theta_2 = 0)
#sim = vector("list", length = 50)
for (i in 1:50){
  filename <- paste0("data/2_simulation/false_positive_sd/artiodactyla/sim_",
                     i, "/history.Rda")
  load(filename)
  sim <- simulateContinuous(history, c(T,T,T))
  par[i,1] <- unname(sim[[2]][which(names(sim[[2]]) == "0")])
  par[i,2] <- unname(sim[[2]][which(names(sim[[2]]) == "1")])
  par[i,3] <- unname(sim[[2]][which(names(sim[[2]]) == "2")])
  par[i,4] <- unname(sim[[3]][which(names(sim[[3]]) == "0")])
  par[i,5] <- unname(sim[[3]][which(names(sim[[3]]) == "1")])
  par[i,6] <- unname(sim[[3]][which(names(sim[[3]]) == "2")])
  par[i,7] <- unname(sim[[4]][which(names(sim[[4]]) == "0")])
  par[i,8] <- unname(sim[[4]][which(names(sim[[4]]) == "1")])
  par[i,9] <- unname(sim[[4]][which(names(sim[[4]]) == "2")])
  
  this_dir <- paste0("data/2_simulation/false_positive_sd/artiodactyla/sim_", i)
  write.nexus.data(sim[[1]], file = paste0(this_dir, "/continuous.nex"),
                   format="Continuous")
}

save(par, file="data/2_simulation/false_positive_sd/artiodactyla/par.Rda")


#simulateContinuous_bm = function(tree, rootState=0) {
#  preorder <- rev(postorder(tree))
#  edges <- tree$edge
#  root_node <- length(tree$tip.label) + 1
#  state = tree$node.states[root_node]
#  mu_at_nodes <- rep(0, length(tree$node.states))
#  
#  sigma2 <- drawSigma2(state_dependent = T)
#  mu_at_nodes[root_node] <- rnorm(1, mean = rootState,
#                                    sd = sqrt(sigma2[[state]]))
#  
#  
#  for (edge_index in preorder){
#    sub_edges <- tree$maps[[edge_index]]
#    parent_node <- edges[edge_index, 1]
#    mu <- mu_at_nodes[parent_node]
#    for (j in 1:length(sub_edges)) {
#      state <- names(sub_edges[j])
#      sigma2 <- drawSigma2(state_dependent = T)
#      v <- sigma2[[state]] * sub_edges[[j]]
#      y <- rnorm(n=1, mu, sqrt(v))  
#      mu <- y
#    }
#    desc_node <- edges[edge_index, 2]
#    mu_at_nodes[desc_node] <- y
#  }
#  
#  cont_list <- list()
#  for (i in 1:length(tree$tip.label)){
    tip <- tree$tip.label[i]
    cont_list[[tip]] <- mu_at_nodes[i]
  }
#  
#  return(cont_list)
#}
#
#
#
#
#
#models   = c("bm", "xxx", "xxt", "xsx", "axx", "asx", "axt", "xst", "ast")
#reps       = 50
#
#grid = expand.grid(models=models, tree=1:reps, stringsAsFactors=FALSE)
#bar = txtProgressBar(style=3, width=40)
#
#
###have to mark down parameters used to simulate somehow?
#for(i in 1:nrow(grid)) {
#  
#  this_row = grid[i,]
#  this_model    = this_row[[1]]
#  this_tree        = this_row[[2]]
#
#  # read the history
#  this_dir = paste0("data/2_simulation/2a_state_dependency/", this_model, "/t", this_tree)
#  load(paste0(this_dir, "/history.Rda"))
#  
#  if (this_model == "bm"){
#    cont_states <- simulateContinuous_bm(tree=history, rootState=0)
#  } else {
#    state_dependencies <- rep(F, 3)
#    if (isTRUE(grepl("a", this_model))){
#      state_dependencies[1] <- T
#    }
#    
#    if (isTRUE(grepl("s", this_model))){
#      state_dependencies[2] <- T
#    }
#    
#    if (isTRUE(grepl("t", this_model))){
#      state_dependencies[3] <- T
#    }
#    
#    cont_states <- simulateContinuous(tree = history, stateDependencies = state_dependencies)
#  }
#  
#  write.nexus.data(cont_states, 
#                   file = paste0(this_dir, "/continuous.nex"),
#                   format="Continuous")
#  setTxtProgressBar(bar, i / nrow(grid))
#}
#cat("\n")


