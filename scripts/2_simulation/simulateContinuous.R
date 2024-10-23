library(ape)
library(phytools)
library(geiger)
library(TESS)
library(tidyverse)

drawHalflife <- function(linked, state_dependent, num_state, root_age){
  if (isTRUE(linked)){
    # can it be extended to more than 2 states?
    halflife <- c(runif(1, 0.1*root_age, root_age))
    if(num_state == 2){
      x = rnorm(1, 0, 0.2)
      while (abs(x) > halflife[1] | halflife[1]+x > 1 ){
        x = rnorm(1, 0, 0.2)
      }
      halflife[2] <- halflife[1] + x
    } else if (num_state == 3){
      x <- rnorm(1, 0, 0.2)
      while (abs(x) > halflife[1] | halflife[1]+x > 1 | halflife[1]-x < 0.05){
        x = rnorm(1, 0, 0.2)
      }
      halflife[2] <- halflife[1] + x
      halflife[3] <- halflife[1] - x
    } else {
      print(paste("Linked stv only supports 2 and 3-state currently."))
    }
  } else {
    if(state_dependent == T){
      halflife <- runif(n=num_state, 0.1*root_age, root_age)
    } else {
      #halflife <- rep(rlnorm(n=1, meanlog=4.349757, sdlog=1.044495), 3)
      halflife <- rep(runif(n=1, 0.1*root_age, root_age), 3)
    }
  }
  names(halflife) = c(1:(num_state))
  return(halflife)
}

drawStv <- function(linked, state_dependent, num_state){
  #not yet finished
  if (isTRUE(linked)){
    # can it be extended to more than 2 states?
    stv <- c(runif(1, -8, 8))
    if(num_state == 2){
      stv[2] <- stv[1] + rnorm(1, 0, 4)
    } else if (num_state == 3){
      x <- rnorm(1, 0, 4)
      stv[2] <- stv[1] + x
      stv[3] <- stv[1] - x
    } else {
      print(paste("Linked stv only supports 2 and 3-state currently."))
    }
  } else {
    if(isTRUE(state_dependent)){
      stv <- rlnorm(n=num_state, meanlog=log(12.39783716), sdlog=0.587405)
    } else {
      stv <- rep(rlnorm(n=1, meanlog=log(12.39783716), sdlog=0.587405), num_state)
    }
  }
  names(stv) = c(1:(num_state))
  return(stv)
}


drawTheta <- function(linked, state_dependent, num_state){
  if (isTRUE(linked)){
    # can it be extended to more than 2 states?
    theta <- c(runif(1, -10, 10))
    if(num_state == 2){
      x <- rnorm(1, 0, 4)
      #while (abs(x+theta[1]) > 8){
      #  x <- rnorm(1, 0, 4)
      #}
      theta[2] <- theta[1] + x
    } else if (num_state == 3){
      x <- rnorm(1, 0, 4)
      #while (abs(x+theta[1]) > 8 | abs(x-theta[1]) > 8){
      #  x <- rnorm(1, 0, 4)
      #}
      theta[2] <- theta[1] + x
      theta[3] <- theta[1] - x
    } else {
      print(paste("Linked theta only supports 2 and 3-state currently."))
    }
  } else {
    if(isTRUE(state_dependent)){
      theta <- c(runif(n=num_state, -8, 8))
    } else {
      theta <- rep(runif(n=1, -8, 8), num_state)
    }
  }
  names(theta) = c(1:(num_state))
  return(theta)
}

# spread of cont trait also depends on alpha --> different priors for theta in inference?



simulateContinuous = function(tree, alpha, sigma2, theta) {
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
  
  return(cont_list)
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

# sim for false positive
num_sim = 1000
num_state = 3
dir_in = "data/2_simulation/triState/"
dir_out = "data/2_simulation/triState/"
par_values <- createParTable(num_state)

bar = txtProgressBar(style=3, width=40)
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
  setTxtProgressBar(bar, i / num_sim)
} 

write.csv(par_values, file="data/2_simulation/triState/pars_1000.csv")



# simulation for power of theta

num_sim = 200
num_state = 2
emp = 12.39783716
stv = c(0.5 * emp, emp, 2 * emp)
halflife = c(0.1, 0.3, 0.6)
root_age = max(node.depth.edgelength(tree))
grid = expand.grid(sim=1:num_sim, stv=stv, halflife=halflife)
grid$alpha = log(2) / grid$halflife
grid$sigma2 = 2 * grid$alpha * grid$stv
grid$rho = 1 - ( 1 - exp( -2 * grid$alpha * root_age ) ) / ( 2 * grid$alpha * root_age )

dir_in = "data/2_simulation/power_theta/"
dir_out = "data/2_simulation/power_theta/"

bar = txtProgressBar(style=3, width=40)
for (i in 1:nrow(grid)){
  file_in <- paste0(dir_in, "sim_",
                    i, "/history.Rda")
  load(file_in)
  
  this_theta = drawTheta(linked=T, num_state=num_state)
  for (j in 1:num_state){
    theta_state = paste0("theta_", j)
    grid[[theta_state]][i] <- this_theta[j]
  }
  this_alpha <- rep(grid$alpha[i], num_state)
  names(this_alpha) = c(1:num_state)
  this_sigma2 <- rep(grid$sigma2[i], num_state)
  names(this_sigma2) = c(1:num_state)

  sim <- simulateContinuous(tree=history, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  grid$delta_cont[i] = max(unlist(sim)) - min(unlist(sim))
  
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/continuous.nex"),
                   format="Continuous")
  
  setTxtProgressBar(bar, i / nrow(grid))
  
} 

write.csv(grid, file="data/2_simulation/power_theta/sim_pars.csv")

