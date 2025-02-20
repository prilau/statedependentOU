library(ape)
library(phytools)
library(geiger)
library(TESS)
library(tidyverse)
source("scripts/5_miscellaneous/functions.R")

# sim for false positive
num_sim = 1000
num_state = 2
tree <- read.tree("data/2_simulation/mammal_diet_height1_n500.tre")
root_age = max(node.depth.edgelength(tree))
emp = 12.39783716
grid = expand.grid(sim=1:num_sim, stv=NA, halflife=NA,
                   theta=NA, delta_cont=NA, var_cont=NA)

dir_in = "data/2_simulation/false_positive/"
dir_out = "data/2_simulation/false_positive/"

bar = txtProgressBar(style=3, width=40)
for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_",
                     i, "/history.Rda")
  load(file_in)
  
  this_halflife <- drawHalflife(state_dependent=F, num_state, root_age)
  grid$halflife[i] <- this_halflife[[1]]
  this_alpha <- log(2) / this_halflife
  this_stv <- drawStv(state_dependent=F, num_state, emp)
  grid$stv[i] <- this_stv[[1]]
  this_sigma2 <- 2 * this_alpha * this_stv
  this_theta <- drawTheta(linked=F, state_dependent=F, num_state)
  grid$theta[i] <- this_theta[[1]]
  
  sim <- simulateContinuous(tree=history, alpha=this_alpha, sigma2=this_sigma2,
                               theta=this_theta)
  grid$delta_cont[i] <- max(unlist(sim)) - min(unlist(sim))
  grid$var_cont[i] = var(unlist(sim))
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/continuous.nex"),
                   format="Continuous")
  setTxtProgressBar(bar, i / num_sim)
} 
grid$alpha = log(2) / grid$halflife
grid$sigma2 = 2 * grid$alpha * grid$stv
grid$rho = 1 - ( 1 - exp( -2 * grid$alpha * root_age ) ) / ( 2 * grid$alpha * root_age )
write.csv(grid, file=paste0(dir_out, "pars.csv"))


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
grid = expand.grid(sim=1:num_sim, stv=stv, theta_1=theta_1, halflife_1=NA,
                   halflife_2=NA, alpha_1=NA, alpha_2=NA, sigma2_1=NA, sigma2_2=NA)
grid$theta_2 = -grid$theta_1

dir_in = "data/2_simulation/power_alpha/"
dir_out = "data/2_simulation/power_alpha/"

bar = txtProgressBar(style=3, width=40)
for (i in 1:nrow(grid)){
  file_in <- paste0(dir_in, "sim_",
                    i, "/history.Rda")
  load(file_in)
  
  # column 4:5 are halflifes
  grid[i,4:5] = drawHalflife(state_dependent = T, num_state=num_state, root_age=root_age)
  # column 6:7 are alphas
  grid[i,6:7] = log(2) / c(grid$halflife_1[i], grid$halflife_2[i])
  # column 8:9 are sigma2s
  grid[i,8:9] <- 2 * this_alpha * grid$stv[i]
  
  this_alpha = c(grid$alpha_1[i], grid$alpha_2[i])
  this_sigma2 <- c(grid$sigma2_1[i], grid$sigma2_2[i])
  this_theta <- c(grid$theta_1[i], grid$theta_2[i])
  names(this_alpha) = names(this_theta) = names(this_sigma2) = c(1:num_state)
  
  sim <- simulateContinuous(tree=history, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  grid$delta_cont[i] = max(unlist(sim)) - min(unlist(sim))
  grid$var_cont[i] = var(unlist(sim))
  
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
  
  this_stv = drawStv(state_dependent=T, num_state=num_state, emp=emp)
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
