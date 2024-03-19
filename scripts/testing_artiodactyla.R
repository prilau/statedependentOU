library(ape)
library(pracma)
library(phytools)
library(tibble) 
library(tidyverse)


# write data
#parameter_1000 <- tibble("alpha_Br" = rep(0,1000),
#                         "alpha_Gr" = rep(0,1000),
#                         "alpha_MF" = rep(0,1000),
#                         "sigma2_Br" =rep(0,1000),
#                         "sigma2_Gr" =rep(0,1000),
#                         "sigma2_MF" =rep(0,1000),
#                         "theta_Br" = rep(0,1000),
#                         "theta_Gr" = rep(0,1000),
#                         "theta_MF" = rep(0,1000))
#
#for (i in 1:1000){
#  parameter_1000[["alpha_Br"]][i] = rgamma(n=1, shape=1, rate=10)
#  parameter_1000[["alpha_Gr"]][i] = rgamma(n=1, shape=1, rate=10)
#  parameter_1000[["alpha_MF"]][i] = rgamma(n=1, shape=1, rate=10)
#  parameter_1000[["sigma2_Br"]][i] = rgamma(n=1, shape=2, rate=10)
#  parameter_1000[["sigma2_Gr"]][i] = rgamma(n=1, shape=2, rate=10)
#  parameter_1000[["sigma2_MF"]][i] = rgamma(n=1, shape=2, rate=10)
#  parameter_1000[["theta_Br"]][i] = rnorm(1, mean = 0, sd = 3)
#  parameter_1000[["theta_Br"]][i] = rnorm(1, mean = 0, sd = 3)
#  parameter_1000[["theta_MF"]][i] = rnorm(1, mean = 0, sd = 3)
#}
#
#write.csv(parameter_1000, "data/1_validation/testing_artiodactyla/ou_parameters_all.csv")
#
#tree <- artiodactyla
#brain <- neocortex$brain_mass_g_log_mean
#names(brain) <- tree$tip.label
#continuous_char <- list()
#for (i in 1:length(brain)){
#  sp = names(brain)[i]
#  continuous_char[[sp]] = brain[[i]]
#}
#
#diet <- as.character(neocortex$diet)
#names(diet) <- neocortex$species
#disc_char <- list()
#for (i in 1:length(diet)){
#  sp = names(diet)[i]
#  disc_char[[sp]] = diet[[i]]
#}
#
#write.nexus.data(continuous_char, "data/1_validation/testing_artiodactyla/continuous_log_brain.nex", format = "continuous")
#write.nexus.data(disc_char, "data/1_validation/testing_artiodactyla/discrete_diet.nex", format = "standard")
#
#for (i in 1:1000){
#  filename = paste0("data/1_validation/artiodactyla_trees/artiodactyla_", i, ".tre")
#  tree <- make.simmap(artiodactyla, diet)
#  write.simmap(tree, file = filename)
#}


# copy of state-dependent pruning function
sd_postorder <- function(node_index, edge, tree, continuousChar,
                         μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta){
  ntip = length(tree$tip.label)
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- sd_postorder(left, edge, tree, continuousChar,
                                μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta)
    μ <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- sd_postorder(right, edge, tree, continuousChar,
                                 μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta)
    μ <- output_right[[1]]
    V <- output_right[[2]]
    log_norm_factor <- output_right[[3]]
    
    
    sub_bl_left = subedges_lengths[[left_edge]] # all subedges of left child edge
    sub_bl_right = subedges_lengths[[right_edge]] # all subedges of right child edge
    
    # for the sake of readability, computation of variance, mean, and log_nf are done in separate loops
    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
    ## Is 'delta_left* exp(2.0 * alpha * bl_left)' added in each sub-edge?
    
    delta_left = V[left]
    v_left = 0 # initialise v_left
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left[i])
      v_left = sigma2[[state]]/(2*alpha[[state]]) *expm1(2.0*alpha[[state]]
                                                         *sub_bl_left[[i]])
      delta_left = v_left + delta_left * exp(2.0 * alpha[[state]] * sub_bl_left[[i]])
    }
    
    delta_right = V[right]
    v_right = 0 # initialise v_right
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right[i])
      v_right = sigma2[[state]]/(2*alpha[[state]]) *expm1(2.0*alpha[[state]]*sub_bl_right[[i]])
      delta_right = v_right + delta_right * exp(2.0 * alpha[[state]] * sub_bl_right[[i]])
    }
    
    var_left = delta_left
    var_right = delta_right
    
    # 2) mean of the normal variable
    mean_left = μ[left]
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left[i])
      mean_left = exp(alpha[[state]]*sub_bl_left[[i]])*(mean_left - theta[[state]]) + theta[[state]]
    }
    
    mean_right = μ[right]
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right[i])
      mean_right = exp(alpha[[state]]*sub_bl_right[[i]])*(mean_right - theta[[state]]) + theta[[state]]
    }
    
    
    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    μ[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node
    
    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    ## this is the problem. I think in RevBayes we compute log_nf with the oldest sub-edge only
    log_nf_left = 0
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left[i])
      log_nf_left = log_nf_left + sub_bl_left[[i]] * alpha[[state]]
    }
    
    log_nf_right = 0
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right[i])
      log_nf_right = log_nf_right + sub_bl_right[[i]] * alpha[[state]]
    }
    
    contrast = mean_left - mean_right
    a = -(contrast*contrast / (2*(var_left+var_right)))
    b = log(2*pi*(var_left+var_right))/2.0
    #b = log(2*pi)/2.0 + log(var_left+var_right)/2.0
    log_nf = log_nf_left + log_nf_right + a - b
    log_norm_factor[node_index] = log_nf
    
    return(list(μ, V, log_norm_factor))
  }
  
  
  # if is tip
  else{
    species = tree$tip.label[node_index]
    
    μ[node_index] = as.numeric(continuousChar[[which(names(continuousChar) == species)]])
    V[node_index] = 0.0 ## if there is no observation error
    
    return(list(μ, V, log_norm_factor))
  }
}

sd_logL_pruning <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes
  
  V = numeric(max_node_index)
  μ = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)
  
  subedges_lengths = tree$maps
  
  root_index = ntip + 1
  
  output <- sd_postorder(root_index, edge, tree, continuousChar,
                         μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta)
  μ <- output[[1]]
  V <- output[[2]]
  log_norm_factor <- output[[3]]
  
  ## assume root value equal to theta
  μ_root = μ[root_index]
  v_root = V[root_index]
  left_edge_from_root <- which(edge[,1] == root_index)[1] # obtain left child edge index of root node
  left_subedges_from_root <- subedges_lengths[[left_edge_from_root]] # obtain sub-edge lengths
  ### note here
  root_state = names(head(left_subedges_from_root, 1)) # obtain root state, assuming it equals last state at left child edge
  lnl = dnorm(theta[[root_state]], mean = μ_root, sd = sqrt(v_root), log = TRUE)
  
  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}





# read data
data("artiodactyla")
data("neocortex")
neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
tree <- artiodactyla
brain <- neocortex$brain_mass_g_log_mean
names(brain) <- tree$tip.label



all_trees <- read.simmap("data/1_validation/testing_artiodactyla/artiodactyla_all.tre", format="phylip",version=1)
parameter_csv <- read.csv("data/1_validation/testing_artiodactyla/ou_parameters_all.csv")

pruning_likelihoods = c()
for (i in 1:10){
  tree <- all_trees[[i]]
  alpha = c(parameter_csv$alpha_Br[i], parameter_csv$alpha_Gr[i], parameter_csv$alpha_MF[i])
  names(alpha) = c("Br", "Gr", "MF")
  sigma2 = c(parameter_csv$sigma2_Br[i], parameter_csv$sigma2_Gr[i], parameter_csv$sigma2_MF[i])
  names(sigma2) = c("Br", "Gr", "MF")
  theta = c(parameter_csv$theta_Br[i], parameter_csv$theta_Gr[i], parameter_csv$theta_MF[i])
  names(theta) = c("Br", "Gr", "MF")
  pruning_likelihoods[i] <- sd_logL_pruning(tree, brain, alpha, sigma2, theta)
}

pruning_likelihoods
