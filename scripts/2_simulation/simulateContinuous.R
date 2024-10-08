library(ape)
library(phytools)
library(geiger)
library(TESS)
library(tidyverse)

drawHalflife <- function(state_dependent, num_state, root_age){
  if(state_dependent == T){
    halflife <- runif(n=num_state, 0.1*root_age, root_age)
  } else {
    #halflife <- rep(rlnorm(n=1, meanlog=4.349757, sdlog=1.044495), 3)
    halflife <- rep(runif(n=1, 0.1*root_age, root_age), 3)
  }
  names(halflife) = c(1:(num_state))
  return(halflife)
}

drawStv <- function(state_dependent, num_state){
  if(state_dependent == T){
    stv <- rlnorm(n=num_state, meanlog=log(12.39783716), sdlog=0.587405)
  } else {
    stv <- rep(rlnorm(n=1, meanlog=log(12.39783716), sdlog=0.587405), 3)
  }
  names(stv) = c(1:(num_state))
  return(stv)
}


drawTheta <- function(state_dependent, num_state){
  if(state_dependent == T){
    theta <- c(rnorm(n=num_state, 0, 4))
  } else {
    theta <- rep(runif(n=1, -10, 10), 3)
  }
  names(theta) = c(1:(num_state))
  return(theta)
}

simulateContinuous = function(tree, stateDependencies=c(halflife=T, stv=T, theta=T), num_state) {
  root_age = max(node.depth.edgelength(tree))
  
  ## Re-parameterization
  #alpha <- log(2) / halflife
  #sigma2 <- stationaryvar * 2 * alpha
  if (isTRUE(stateDependencies["halflife"])){
    halflife <-  drawHalflife(state_dependent = T, num_state, root_age)
  } else {
    halflife <-  drawHalflife(state_dependent = F, num_state, root_age)
  }
  
  if (isTRUE(stateDependencies["stv"])){
    stv <-  drawStv(state_dependent = T, num_state)
  } else {
    stv <-  drawStv(state_dependent = F, num_state)
  }
  
  if (isTRUE(stateDependencies["theta"])){
    theta <-  drawTheta(state_dependent = T, num_state)
  } else {
    theta <-  drawTheta(state_dependent = F, num_state)
  }
  
  alpha <- log(2) / halflife
  sigma2 <- 2 * alpha * stv
  rho <- 1 - ( 1 - exp( -2 * alpha * root_age ) ) / ( 2 * alpha * root_age )
  
  preorder <- rev(postorder(tree))
  edges <- tree$edge
  root_node <- length(tree$tip.label) + 1
  state = tree$node.states[root_node]
  mu_at_nodes <- rep(0, length(tree$node.states))
  mu_at_nodes[root_node] <- theta[[state]]

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
  
  return(list(alpha, halflife, sigma2, stv, rho, theta, cont_list))
}


createParTable <- function(num_state){
  pars = c("alpha", "halflife", "sigma2", "stv", "rho", "theta")
  headers <- c()
  for (par in pars){
    for (i in 1:num_state){
      headers <- append(headers, paste0(c(par, i), collapse = "_"))
    }
  }
  headers <- append(headers, c("sim", "state"))
  
  par_values <- tibble()
  for (header in headers){
    par_values[[header]] = 0
  }
  return(par_values)
}

enterParTable <- function(par_table, sim_sd, sim_sx, num_state, sim){
  j=sim*2-1
  k=sim*2
  
  for (m in 1:num_state){
    par_values[j,0*num_state+m] <- unname(sim_sd[[1]][m])
    par_values[j,1*num_state+m] <- unname(sim_sd[[2]][m])
    par_values[j,2*num_state+m] <- unname(sim_sd[[3]][m])
    par_values[j,3*num_state+m] <- unname(sim_sd[[4]][m])
    par_values[j,4*num_state+m] <- unname(sim_sd[[5]][m])
    par_values[j,5*num_state+m] <- unname(sim_sd[[6]][m])
    
    par_values[k,0*num_state+m] <- unname(sim_sx[[1]][m])
    par_values[k,1*num_state+m] <- unname(sim_sx[[2]][m])
    par_values[k,2*num_state+m] <- unname(sim_sx[[3]][m])
    par_values[k,3*num_state+m] <- unname(sim_sx[[4]][m])
    par_values[k,4*num_state+m] <- unname(sim_sx[[5]][m])
    par_values[k,5*num_state+m] <- unname(sim_sx[[6]][m])
  }
  
  par_values$sim[j] = par_values$sim[k] = i
  par_values$state[j] <- "sd"
  par_values$state[k] <- "sx"
  return(par_values)
}


num_sim = 200
num_state = 3
dir_in = "data/2_simulation/triState/"
dir_out = "data/2_simulation/triState/"

for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_",
                     i, "/history.Rda")
  load(file_in)
  
  sim_sd <- simulateContinuous(history, c(halflife=T, stv=T, theta=T), num_state)
  sim_sx <- simulateContinuous(history, c(halflife=F, stv=F, theta=F), num_state)
  
  par_values <- enterParTable(par_table, sim_sd, sim_sx, num_state, sim=i)
  
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim_sd[[7]], file = paste0(this_dir, "/continuous_sd.nex"),
                   format="Continuous")
  write.nexus.data(sim_sx[[7]], file = paste0(this_dir, "/continuous_sx.nex"),
                   format="Continuous")
} 

write.csv(par_values, file="data/2_simulation/triState/pars.csv")