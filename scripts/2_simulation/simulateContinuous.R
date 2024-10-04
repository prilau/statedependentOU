library(ape)
library(phytools)
library(geiger)
library(TESS)
library(tidyverse)

drawHalflife <- function(state_dependent=T){
  if(state_dependent == T){
    halflife <- rlnorm(n=3, meanlog=4.349757, sdlog=1.044495)
    names(halflife) = sample(c("0", "1", "2"))
  } else {
    halflife <- rep(rlnorm(n=1, meanlog=4.349757, sdlog=1.044495), 3)
    names(halflife) = c("0", "1", "2")
  }
  return(halflife)
}

drawStv <- function(state_dependent=T){
  if(state_dependent == T){
    stv <- rlnorm(n=3, meanlog=log(12.39783716), sdlog=0.587405)
    names(stv) = sample(c("0", "1", "2"))
  } else {
    stv <- rep(rlnorm(n=1, meanlog=log(12.39783716), sdlog=0.587405), 3)
    names(stv) = c("0", "1", "2")
  }
  return(stv)
}


drawTheta <- function(state_dependent=T){
  if(state_dependent == T){
    theta <- c(runif(n=3, -10, 10))
  } else {
    theta <- rep(runif(n=1, -10, 10), 3)
  }
  names(theta) = c("0", "1", "2")
  return(theta)
}








#drawRootState <- function(alphaRoot, sigma2Root, thetaRoot){
#  rootState <- rnorm(n=1, mean=thetaRoot, sd = sigma2Root/(2* alphaRoot))
#  return(rootState)
#}


simulateContinuous = function(tree, stateDependencies=c(halflife=T, stv=T, theta=T)) {

  ## Re-parameterization
  #alpha <- log(2) / halflife
  #sigma2 <- stationaryvar * 2 * alpha
  if (isTRUE(stateDependencies["halflife"])){
    halflife <-  drawHalflife(state_dependent = T)
  } else {
    halflife <-  drawHalflife(state_dependent = F)
  }
  
  if (isTRUE(stateDependencies["stv"])){
    stv <-  drawStv(state_dependent = T)
  } else {
    stv <-  drawStv(state_dependent = F)
  }
  
  if (isTRUE(stateDependencies["theta"])){
    theta <-  drawTheta(state_dependent = T)
  } else {
    theta <-  drawTheta(state_dependent = F)
  }
  
  alpha <- log(2) / halflife
  sigma2 <- 2 * alpha * stv
  
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

pars_sd <- tibble(alpha_0 = 0,  alpha_1 = 0,
                  sigma2_0 = 0, sigma2_1 = 0,
                  theta_0 = 0,  theta_1 = 0)
pars_stateless <- tibble(alpha_0 = 0,  alpha_1 = 0,
                         sigma2_0 = 0, sigma2_1 = 0,
                         theta_0 = 0,  theta_1 = 0)
#sim = vector("list", length = 50)
for (i in 1:3){
  filename <- paste0("data/2_simulation/convergence/sim_",
                     i, "/history.Rda")
  load(filename)
  sim_sd <- simulateContinuous(history, c(halflife=T, stv=T, theta=T))
  sim_stateless <- simulateContinuous(history, c(halflife=F, stv=F, theta=F))
  pars_sd[i,1] <- unname(sim_sd[[2]][which(names(sim_sd[[2]]) == "0")])
  pars_sd[i,2] <- unname(sim_sd[[2]][which(names(sim_sd[[2]]) == "1")])
  pars_sd[i,3] <- unname(sim_sd[[3]][which(names(sim_sd[[3]]) == "0")])
  pars_sd[i,4] <- unname(sim_sd[[3]][which(names(sim_sd[[3]]) == "1")])
  pars_sd[i,5] <- unname(sim_sd[[4]][which(names(sim_sd[[4]]) == "0")])
  pars_sd[i,6] <- unname(sim_sd[[4]][which(names(sim_sd[[4]]) == "1")])

  pars_stateless[i,1] <- unname(sim_stateless[[2]][which(names(sim_stateless[[2]]) == "0")])
  pars_stateless[i,2] <- unname(sim_stateless[[2]][which(names(sim_stateless[[2]]) == "1")])
  pars_stateless[i,3] <- unname(sim_stateless[[3]][which(names(sim_stateless[[3]]) == "0")])
  pars_stateless[i,4] <- unname(sim_stateless[[3]][which(names(sim_stateless[[3]]) == "1")])
  pars_stateless[i,5] <- unname(sim_stateless[[4]][which(names(sim_stateless[[4]]) == "0")])
  pars_stateless[i,6] <- unname(sim_stateless[[4]][which(names(sim_stateless[[4]]) == "1")])
  
  this_dir <- paste0("data/2_simulation/convergence/sim_", i)
  write.nexus.data(sim_sd[[1]], file = paste0(this_dir, "/continuous_sd.nex"),
                   format="Continuous")
  write.nexus.data(sim_stateless[[1]], file = paste0(this_dir, "/continuous_stateless.nex"),
                   format="Continuous")
} 

root_age <- max(node.depth.edgelength(history))
pars_sd <- pars_sd %>%
  mutate(halflife_0 = log(2) / alpha_0,
         halflife_1 = log(2) / alpha_1,
         stv_0 = sigma2_0 / (2 * alpha_0),
         stv_1 = sigma2_1 / (2 * alpha_1),
         rho_0 = 1 - ( 1 - exp( -2 * alpha_0 * root_age ) ) / ( 2 * alpha_0 * root_age ),
         rho_1 = 1 - ( 1 - exp( -2 * alpha_1 * root_age ) ) / ( 2 * alpha_1 * root_age ))
pars_stateless <- pars_stateless %>%
  mutate(halflife_0 = log(2) / alpha_0,
         halflife_1 = log(2) / alpha_1,
         stv_0 = sigma2_0 / (2 * alpha_0),
         stv_1 = sigma2_1 / (2 * alpha_1),
         rho_0 = 1 - ( 1 - exp( -2 * alpha_0 * root_age ) ) / ( 2 * alpha_0 * root_age ),
         rho_1 = 1 - ( 1 - exp( -2 * alpha_1 * root_age ) ) / ( 2 * alpha_1 * root_age ))
save(pars_sd, file="data/2_simulation/convergence/pars_sd.Rda")
save(pars_stateless, file="data/2_simulation/convergence/pars_stateless.Rda")



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


