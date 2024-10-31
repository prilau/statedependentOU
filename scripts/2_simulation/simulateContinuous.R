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

drawStv <- function(linked, state_dependent, num_state, emp){
  #not yet finished
  if (isTRUE(linked)){
    # can it be extended to more than 2 states?
    stv <- c(runif(1, 0.5*emp, 2*emp))
    if(num_state == 2){
      x <- rnorm(1, 0, 0.5*emp)
      while (stv[1]-abs(x) < 0){
        x = rnorm(1, 0, 0.5*emp)
      }
      stv[2] <- stv[1] + x
    } else {
      print(paste("Linked stv only supports binary state currently."))
    }
  } else {
    if(isTRUE(state_dependent)){
      stv <- rlnorm(n=num_state, meanlog=log(emp), sdlog=0.587405)
    } else {
      stv <- rep(rlnorm(n=1, meanlog=log(emp), sdlog=0.587405), num_state)
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


#simulation for exploring power of each OU parameter
tree <- read.tree("data/2_simulation/mammal_diet_height1_n500.tre")
root_age = max(node.depth.edgelength(tree))

num_sim = 200 # #simulations for each combination of background OU parameters
num_state = 2 # #discrete states
emp = 12.39783716 # empirical variance of ln(body size in kg) of 1190 mammals

# simulation for power of theta
stv = c(0.5 * emp, emp, 2 * emp)
halflife = c(0.1, 0.3, 0.6)
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
write.csv(grid, file="data/2_simulation/power_theta/pars.csv")




# simulate continuous traits for power_alpha
stv = c(0.5 * emp, emp, 2 * emp)
theta_1 = c(1, 3, 5)
grid = expand.grid(sim=1:num_sim, stv=stv, theta_1=theta_1,
                   alpha_1=NA, alpha_2=NA, sigma2_1=NA, sigma2_2=NA)
grid$theta_2 = -grid$theta_1

dir_in = "data/2_simulation/power_alpha/"
dir_out = "data/2_simulation/power_alpha/"

bar = txtProgressBar(style=3, width=40)
for (i in 1:nrow(grid)){
  file_in <- paste0(dir_in, "sim_",
                    i, "/history.Rda")
  load(file_in)
  
  this_alpha = log(2) / drawHalflife(linked=T, num_state=num_state, root_age=root_age)
  this_sigma2 <- c()
  for (j in 1:num_state){
    alpha_state = paste0("alpha_", j)
    sigma2_state = paste0("sigma2_", j)
    grid[[alpha_state]][i] <- this_alpha[j]
    grid[[sigma2_state]][i] <- 2 * this_alpha[j] * grid$stv[i]
    this_sigma2 <- append(this_sigma2, grid[[sigma2_state]][i])
  }
  this_theta <- c(grid$theta_1[i], grid$theta_2[i])
  names(this_theta) = c(1:num_state)
  names(this_sigma2) = c(1:num_state)
  
  sim <- simulateContinuous(tree=history, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  grid$delta_cont[i] = max(unlist(sim)) - min(unlist(sim))
  
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/continuous.nex"),
                   format="Continuous")
  
  setTxtProgressBar(bar, i / nrow(grid))
  
} 

grid$rho_1 = 1 - ( 1 - exp( -2 * grid$alpha_1 * root_age ) ) / ( 2 * grid$alpha_1 * root_age )
grid$rho_2 = 1 - ( 1 - exp( -2 * grid$alpha_2 * root_age ) ) / ( 2 * grid$alpha_2 * root_age )

write.csv(grid, file="data/2_simulation/power_alpha/pars.csv")


# simulate continuous traits for power_sigma2
halflife=c(0.1, 0.3, 0.6)
theta_1 = c(1, 3, 5)
grid = expand.grid(sim=1:num_sim, halflife=halflife, theta_1=theta_1,
                   stv_1=NA, stv_2=NA, sigma2_1=NA, sigma2_2=NA)
grid$theta_2 = -grid$theta_1
grid$alpha = log(2) / grid$halflife
grid$rho = 1 - ( 1 - exp( -2 * grid$alpha * root_age ) ) / ( 2 * grid$alpha * root_age )

dir_in = "data/2_simulation/power_sigma2/"
dir_out = "data/2_simulation/power_sigma2/"

bar = txtProgressBar(style=3, width=40)
for (i in 1:nrow(grid)){
  file_in <- paste0(dir_in, "sim_",
                    i, "/history.Rda")
  load(file_in)
  
  this_alpha = c(rep(grid$alpha[i], num_state))
  names(this_alpha) = c(1:num_state)
  this_theta <- c(grid$theta_1[i], grid$theta_2[i])
  names(this_theta) = c(1:num_state)
  
  this_stv = drawStv(linked=T, num_state=num_state, emp=emp)
  grid$stv_1[i] = this_stv[[1]]
  grid$stv_2[i] = this_stv[[2]]
  
  this_sigma2 <- 2 * this_alpha * this_stv
  grid$sigma2_1[i] = this_sigma2[[1]]
  grid$sigma2_2[i] = this_sigma2[[2]]

  
  sim <- simulateContinuous(tree=history, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  grid$delta_cont[i] = max(unlist(sim)) - min(unlist(sim))
  grid$var_cont[i] = var(unlist(sim))
  
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/continuous.nex"),
                   format="Continuous")
  
  setTxtProgressBar(bar, i / nrow(grid))
  
} 

write.csv(grid, file="data/2_simulation/power_sigma2/pars.csv")
