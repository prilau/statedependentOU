library(ape)
library(pracma)
library(phytools)
library(slouch)
library(tibble) 
library(tidyverse)

###################################################
#                                                 #
#                 Stateless pruning               #
#                                                 #
###################################################

## Postorder function (stateless)
postorder <- function(node_index, edge, tree, continuousChar,
                      μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta){
  ntip = length(tree$tip.label)
  #root_node = ntip + 1
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- postorder(left, edge, tree, continuousChar,
                         μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
    μ <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- postorder(right, edge, tree, continuousChar,
                             μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
    μ <- output_right[[1]]
    V <- output_right[[2]]
    log_norm_factor <- output_right[[3]]
    
    
    bl_left = branch_lengths[left_edge] # all branch of left child edge
    bl_right = branch_lengths[right_edge] # all branch of right child edge
    
    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
    
    v_left = sigma2/(2*alpha) *expm1(2.0*alpha*bl_left)
    var_left = v_left + V[left] * exp(2.0 * alpha * bl_left)
    
    v_right = sigma2/(2*alpha) *expm1(2.0*alpha*bl_right)
    var_right = v_right + V[right] * exp(2.0 * alpha * bl_right)
    
    # 2) mean of the normal variable
    mean_left = exp(alpha*bl_left)*(μ[left] - theta) + theta
    mean_right = exp(alpha*bl_right)*(μ[right] - theta) + theta
    
    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    μ[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node
    
    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    log_nf_left = bl_left * alpha
    log_nf_right = bl_right * alpha

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
    #edge_index = which(edge[,2] == node_index) # find edge index by tip node index
    #subedge = tree$maps[[edge_index]]
    species = tree$tip.label[node_index]
    
    μ[node_index] = as.numeric(continuousChar[[which(names(continuousChar) == species)]])
    V[node_index] = 0.0 ## if there is no observation error
    
    return(list(μ, V, log_norm_factor))
  }
}


## Pruning method (stateless)
logL_pruning <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes
  
  V = numeric(max_node_index)
  μ = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)
  
  branch_lengths = tree$edge.length
  
  root_index = ntip + 1
  
  
  
  
  output <- postorder(root_index, edge, tree, continuousChar,
                      μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
  μ <- output[[1]]
  V <- output[[2]]
  log_norm_factor <- output[[3]]
  
  ## assume root value equal to theta
  μ_root = μ[root_index]
  v_root = V[root_index]
  lnl = dnorm(theta, mean = μ_root, sd = sqrt(v_root), log = TRUE) # are \theta and \mu in correct positions?
  
  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}



###################################################
#                                                 #
#            State-dependent pruning              #
#                                                 #
###################################################

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



###################################################
#                                                 #
#                  Stateless vcv                  #
#                                                 #
###################################################

logL_vcv <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip <- length(tree$tip.label)
  mrca1 <- ape::mrca(tree) # get the ancestral node label for each pair of tips
  times <- ape::node.depth.edgelength(tree) # get time at each node from root
  ta <- matrix(times[mrca1], nrow=ntip, dimnames = list(tree$tip.label, tree$tip.label)) # get time of divergence for each pair of tips
  T.term <- times[1:ntip] # get time at tips
  tia <- times[1:ntip] - ta
  tja <- t(tia)
  tij <- tja + tia # distance in time unit between two tips
  
  vy = sigma2 / (2*alpha)
  
  #V = vy * (1 - exp(-2 * alpha * ta)) * exp(-alpha * tij)
  V = vy * -1 * expm1(-2 * alpha * ta) * exp(-alpha * tij) ### ta = time tgt; tij = time not tgt (sum of two branches)
  
  X = matrix(1, ntip)
  
  C = chol(V) # upper triangular matrix
  L = t(C) # lower triangular matrix
  log_det_V = 0
  for (i in 1:ntip){
    log_det_V = log_det_V + log(L[i,i])
  }
  log_det_V = log_det_V *2.0 # equals to julia implementation to 12 sig. fig.
  
  y = NULL
  for (species in tree$tip.label){
    y[species] = as.numeric(continuousChar[species])
  }
  
  r = solve(L) %*% y - solve(L) %*% X * theta # what does de-correlated residuals mean?
  
  # res = - (n/2) * log(2*pi) - 0.5 * log_det_V - 0.5 * dot(r, r)
  #     = exp(-n/2)^(2*pi) * exp(-0.5)^det_V * exp(-0.5)^dot(r, r) ?
  res = 0.0
  res = res - 0.5 * ntip * log(2*pi)
  res = res - 0.5 * log_det_V
  res = res - 0.5 * dot(r, r) # what is r and what is  dot product of r?
  
  return(res)
}


###################################################
#                                                 #
#               State-dependent vcv               #
#                                                 #
###################################################

# find parent node by providing child node
parentNode <- function(tree, x){
  m <- which(tree$edge[, 2] == x)
  return(tree$edge[m, 1])
}

# find nodes along a lineage towards root node by providing initial child (presumably tip) node
nodesAlongLineage <- function(tree, old_node, young_node){
  k <- young_node
  while(young_node != old_node){
    k <- c(k, parentNode(tree, young_node))
    young_node <- tail(k, n = 1)
  }
  return(k)
}

# find subedges of a lineage
lineage.constructor <- function(tree, root_node, e){
  nodes <- nodesAlongLineage(tree, root_node, e)
  edges <- which(tree$edge[,2] %in% nodes) # from root to tip
  subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # tip to root

  state_changes <- names(subedge_lengths) # from tip to root
  #state_changes <- c(state_changes[1], state_changes) # add root state, assuming root state equals the state of the closest subedge
  
  #lineage$state_indicator <- lapply(all_states, function(x) {res <- match(lineage$state_changes, x); res[is.na(res)] <- 0; return(res)})
  #names(lineage$state_indicator) <- all_states
  
  # recording time-related numbers of each subedge (root is a subedge with length = 0)
  #times <- cumsum(unname(subedge_lengths))
  #time_tip <- tail(times, n = 1)
  
  #time_begin <- time_tip - c(0, head(times, n = -1))
  #time_end <- time_tip - times
  #time_span <- time_begin - time_end
  
  return(tibble(state = state_changes,
                #time_begin = time_begin,
                #time_end = time_end,
                time_span = subedge_lengths))
}

# not yet finished - weight matrix function
## need updates
weights.lineage <- function(tree, alpha, e){
  root_node = length(tree$tip.label) + 1
  lineage <- lineage.constructor(tree, root_node, e)
  lineage[["alpha"]] = alpha[lineage[["state"]]]
  
  W = matrix(0, ncol = length(alpha), nrow = 1)
  colnames(W) = sort(names(alpha))
  
  if (length(lineage[[1]]) > 1){
    lineage <- lineage %>%
      mutate(
        exp1 = -1 * expm1(-1 * alpha * time_span),
        sum2_temp = -1 * alpha * time_span)
    lineage$exp1[length(lineage$exp1)] = 1
    lineage$sum2 = 0

    for (i in 2:length(lineage[[1]])){
      lineage$sum2[i] = lineage$sum2_temp[i-1]
      lineage$sum2_temp[i] = lineage$sum2[i] + lineage$sum2_temp[i]
    }

    weights = lineage %>% mutate(exp_final = exp1 * exp(sum2)) %>% 
      group_by(state) %>% 
      summarise(weight = sum(exp_final))
    
    for (i in 1:nrow(weights)){
      W[, weights$state[i]] = weights$weight[i]
    }
  } else {
    W[, weights$state[1]] = 1
  }
  
  return(W)
}

# combine to form weight matrix
weight.matrix <- function(tree, alpha){
  ntip = length(tree$tip.label)
  weight_matrix = matrix(0, nrow = ntip, ncol = length(alpha))
  rownames(weight_matrix) <- tree$tip.label
  colnames(weight_matrix) <- c(sort(names(alpha)))
  for (i in 1:ntip){
    weight_matrix[i, ] <- weights.lineage(tree, alpha, i)
  }
  return(weight_matrix)
}


cov.accum <- function(tree, mrca_node, alpha, sigma2){
  root_node = length(tree$tip.label) + 1
  if (mrca_node == root_node){
    cov_accum = 0.0
  } else {
    nodes <- nodesAlongLineage(tree, root_node, mrca_node)
    edges <- which(tree$edge[,2] %in% nodes) # from root to mcra_node
    subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from mcra_node to root
    
    subedge_lengths <- tibble(state = names(subedge_lengths),
                              time_span = subedge_lengths,
                              alpha = alpha[names(subedge_lengths)],
                              sigma2 = sigma2[names(subedge_lengths)]) %>% 
      mutate(exp1 = -1 * expm1(-2 * alpha * time_span),
             sum2_temp = -2 * alpha * time_span)
    subedge_lengths$sum2= 0
    
    if (length(subedge_lengths[[1]]) == 1){
      subedge_lengths = subedge_lengths %>% 
        mutate(cov = sigma2 / (2 * alpha) * exp1)
      cov_accum = subedge_lengths$cov[[1]]
    } else {
      for (i in 2:length(subedge_lengths[[1]])){
        subedge_lengths$sum2[i] = subedge_lengths$sum2_temp[i-1]
        subedge_lengths$sum2_temp[i] = subedge_lengths$sum2[i] + subedge_lengths$sum2_temp[i]
      }
      cov_accum = subedge_lengths %>% mutate(exp3 = exp1 * exp(sum2)) %>% 
        group_by(state) %>% 
        summarise(sum4 = sum(sigma2 / (2 * alpha) * exp3)) %>% 
        reframe(sum_final = sum(sum4)) %>% 
        unlist() %>% 
        unname()
    }
  }
  return(cov_accum)
}

cov.loss <- function(tree, mrca_node, alpha, tip){
  if (mrca_node == tip){
    cov_loss_rate = 0
  } else {
    nodes <- nodesAlongLineage(tree, mrca_node, tip)
    nodes <- head(nodes, n = -1)
    edges <- which(tree$edge[,2] %in% nodes) # from root to mcra_node
    subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from mcra_node to root
    subedge_lengths <- tibble(time_span = subedge_lengths,
                              alpha = alpha[names(subedge_lengths)])
    cov_loss_rate = subedge_lengths %>% 
      mutate(sum1 = -1 * alpha * time_span) %>% 
      reframe(sum_final = sum(sum1))
  }
  return(cov_loss_rate)
}

vcv.pairwise <- function(tree, alpha, sigma2, tip1, tip2){
  mrca_node <- ape::mrca(tree)[tip1, tip2]
  cov_accum = cov.accum(tree, mrca_node, alpha, sigma2)
  cov_loss1 = cov.loss(tree, mrca_node, alpha, tip1)
  cov_loss2 = cov.loss(tree, mrca_node, alpha, tip2)
  cov = cov_accum * exp(cov_loss1 + cov_loss2)
  return(unlist(unname(cov)))
}


vcv.matrix <- function(tree, alpha, sigma2){
  ntip <- length(tree$tip.label)
  V <- matrix(nrow = ntip, ncol = ntip)
  j = ntip
  while (j != 0){
    for (i in 1:ntip){
      V[i,j] <- vcv.pairwise(tree, alpha, sigma2, i, j)
      V[j,i] <- V[i,j]
    }
    j = j-1
  }
  colnames(V) <- tree$tip.label
  rownames(V) <- tree$tip.label
  return(V)
}


sd_logL_vcv <- function(tree, continuousChar, alpha, sigma2, theta){
  alpha = alpha[sort(names(alpha))]
  sigma2 = sigma2[sort(names(sigma2))]
  theta = theta[sort(names(theta))]
  theta = as.matrix(theta, nrow = 3)
  
  ntip <- length(tree$tip.label)
  V = vcv.matrix(tree, alpha, sigma2)
  
  W = weight.matrix(tree, alpha)
  
  C = chol(V) # upper triangular matrix
  L = t(C) # lower triangular matrix
  log_det_V = 0
  for (i in 1:ntip){
    log_det_V = log_det_V + log(L[i,i])
  }
  log_det_V = log_det_V * 2.0 # equals to julia implementation to 12 sig. fig.
  
  y = NULL
  for (species in tree$tip.label){
    y[species] = as.numeric(continuousChar[species])
  }
  
  # inverse of L
  r = solve(L) %*% y - solve(L) %*% W %*% theta # what does de-correlated residuals mean?
  
  # res = - (n/2) * log(2*pi) - 0.5 * log_det_V - 0.5 * dot(r, r)
  #     = exp(-n/2)^(2*pi) * exp(-0.5)^det_V * exp(-0.5)^dot(r, r) ?
  res = 0.0
  res = res - (ntip/2) * log(2*pi)
  res = res - 0.5 * log_det_V
  res = res - 0.5 * dot(r, r) # is it dot product? what is dot product of r?
  
  return(res)
}



###################################################
#                                                 #
#                 Load test data                  #   
#                                                 #
###################################################

# From local
#dummy_tree <- read.simmap("data/1_validation/dummy_threetaxon_simmap2.tre",
#                          format="phylip",version=1)

# read.nexus.data() gives lists of species with character values, need to change to input format for our functions
#char <- read.nexus.data("data/1_validation/dummy_threetaxon_Continuous.nex")
#continuousChar <- c()
#for (species in char){
#  continuousChar <- append(continuousChar, as.numeric(species))
#}
#names(continuousChar) <- names(char)


# From slouch
data("artiodactyla")
data("neocortex")
neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]

## Guideline for stateless methods
# tree: a tree
# continuousChar: a vector with names of elements == species name
# alpha, sigma2, theta: scalars

# Example
tree <- artiodactyla
brain <- neocortex$brain_mass_g_log_mean
names(brain) <- tree$tip.label


## Guideline for state-dependent methods
# tree: a mapped tree
# continuousChar: a vector with names of elements == species name
# alpha, sigma2, theta: vectors with names of elements == names of discrete character states

# Example
diet <- as.character(neocortex$diet)
names(diet) <- neocortex$species
discrete_states <- unique(diet)
set.seed(123)
tree <- make.simmap(artiodactyla, diet)
plot(tree)

#alpha = c(rep(rgamma(n=1, shape=1, rate=10), 3))
alpha = rgamma(n=3, shape=1, rate=10)
names(alpha) = discrete_states

#sigma2 = c(rep(rgamma(n=1, shape=2, rate=10), 3))
sigma2 = rgamma(n=3, shape=2, rate=10)
names(sigma2) = discrete_states

#theta = c(rep(rnorm(1, mean = 0, sd = 3), 3))
theta = rnorm(3, mean = 0, sd = 3)
names(theta) = discrete_states




# RUN
logL_vcv(tree, brain, alpha[[1]], sigma2[[1]], theta[[1]])
logL_pruning(tree, brain, alpha[[1]], sigma2[[1]], theta[[1]])

sd_logL_vcv(tree, brain, alpha, sigma2, theta)
sd_logL_pruning(tree, brain, alpha, sigma2, theta)


###################################################
#                                                 #
#                     Plots                       #   
#                                                 #
###################################################

likelihood_difference = c()
for (i in 1:20){
  alpha = rgamma(n=3, shape=1, rate=10)
  sigma2 = rgamma(n=3, shape=2, rate=10)
  theta = rnorm(3, mean = 0, sd = 3)
  names(alpha) = discrete_states
  names(sigma2) = discrete_states
  names(theta) = discrete_states
  l1 = sd_logL_vcv(tree, brain, alpha, sigma2, theta)
  l2 = sd_logL_pruning(tree, brain, alpha, sigma2, theta)
  likelihood_difference[i] = l1 - l2
}

ggplot(as.data.frame(likelihood_difference)) +
  geom_point(aes(x = 1:20, y = likelihood_difference))
