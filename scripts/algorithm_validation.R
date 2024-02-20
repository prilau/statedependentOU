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
                      μ, V, log_norm_factor, branch_lengths, sigma2, alpha, theta){
  ntip = length(tree$tip.label)
  #root_node = ntip + 1
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- postorder(left, edge, tree, continuousChar,
                         μ, V, log_norm_factor, branch_lengths, sigma2, alpha, theta)
    μ <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- postorder(right, edge, tree, continuousChar,
                             μ, V, log_norm_factor, branch_lengths, sigma2, alpha, theta)
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
logL_pruning <- function(tree, continuousChar, sigma2, alpha, theta){
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
                      μ, V, log_norm_factor, branch_lengths, sigma2, alpha, theta)
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
                      μ, V, log_norm_factor, subedges_lengths, sigma2, alpha, theta){
  ntip = length(tree$tip.label)

  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- sd_postorder(left, edge, tree, continuousChar,
                             μ, V, log_norm_factor, subedges_lengths, sigma2, alpha, theta)
    μ <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- sd_postorder(right, edge, tree, continuousChar,
                              μ, V, log_norm_factor, subedges_lengths, sigma2, alpha, theta)
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

sd_logL_pruning <- function(tree, continuousChar, sigma2, alpha, theta){
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
                         μ, V, log_norm_factor, subedges_lengths, sigma2, alpha, theta)
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

logL_vcv <- function(tree, continuousChar, sigma2, alpha, theta){
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
nodesAlongLineage <- function(tree, x){
  k <- x
  N <- length(tree$tip.label)
  while(x != N + 1){
    k <- c(k, parentNode(tree, x))
    x <- tail(k, n = 1)
  }
  return(k)
}

# find subedges of a lineage
lineage.constructor <- function(tree, e){
  nodes <- nodesAlongLineage(tree, e)
  
  ## Simmap splits up each edge into sub-edges, depending on the split. So, we use edges instead of nodes, and introduce sub-edges
  edges <- which(tree$edge[,2] %in% nodes) # from tip to root
  subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from tip to root

  #all_states <- colnames(tree$mapped.edge)

  state_changes <- names(subedge_lengths) # from tip to root
  state_changes <- c(state_changes, state_changes[length(state_changes)]) # add root state, assuming root state equals the state of the closest subedge
  
  #lineage$state_indicator <- lapply(all_states, function(x) {res <- match(lineage$state_changes, x); res[is.na(res)] <- 0; return(res)})
  #names(lineage$state_indicator) <- all_states
  
  # recording time-related numbers of each subedge (root is a subedge with length = 0)
  time_point <- cumsum(c(0, unname(subedge_lengths))) # tip
  time_begin <- c(tail(time_point, n = -1), tail(time_point, n = 1)) ## older end of subedge
  time_end <- time_point ## younger end of subedge
  time_span <- c(unname(subedge_lengths), 0) # time at root and each subedge
  
  
  return(tibble(state_changes = state_changes,
                time_begin = time_begin,
                time_end = time_end,
                time_span = time_span))
}

# not yet finished - weight matrix function
weights.lineage <- function(tree, named_alpha, e){
  lineage <- lineage.constructor(tree, e)
  lineage[["alpha"]] = named_alpha[lineage[["state_changes"]]]
  lineage <- lineage %>% mutate(exp_1 = exp(-alpha * time_end) - exp(-alpha * time_begin),
                     exp_2 = alpha * lineage$time_span)
  ancestral_state <- lineage$state_changes[length(lineage$state_changes)]
  weight_ancestral = lineage %>%
    summarise(ancestral = exp(-1 * sum(exp_2))) %>% 
    unlist()
  names(weight_ancestral) <- ancestral_state
  
  weights <- lineage %>%
    group_by(state_changes) %>% 
    summarise(weight = sum(exp_1) * exp(-1 * sum(exp_2)))
  weights$weight[which(weights$state_changes == ancestral_state)] = weights$weight[which(weights$state_changes == ancestral_state)] + weight_ancestral
  
  
  #weight_states <- weights$weight
  #names(weight_states) <- weights$state_changes
  #weight_states[which(names(weight_states) == ancestral_state)] = weight_states[which(names(weight_states) == ancestral_state)] + weight_ancestral
  
  #weight_matrix = weight_matrix / sum(weight_matrix)
  weight_matrix <- matrix(nrow = length(named_alpha))
  rownames(weight_matrix) <- c(names(named_alpha))
  for (rowname in rownames(weight_matrix)){
    if (rowname %in% weights$state_changes){
      weight_matrix[rowname,] = weights$weight[which(weights$state_changes == rowname)]
    }
    else {
      weight_matrix[rowname,] = 0
    }
  }
  
  # normalise the weights so that sum(weights) == 1
  weight_matrix = weight_matrix/sum(weight_matrix)
  
  return(weight_matrix)
}

# combine to form weight matrix
weight.matrix <- function(tree, named_alpha){
  ntip = length(tree$tip.label)
  weight_matrix = matrix(nrow = ntip, ncol = length(named_alpha))
  rownames(weight_matrix) <- tree$tip.label
  colnames(weight_matrix) <- c(names(named_alpha))
  for (i in 1:ntip){
    weight_matrix[i,] <- weights.lineage(tree, named_alpha, i)
  }
  return(weight_matrix)
}



nodesBeforeDiverge <- function(tree, tip1, tip2){
  k <- ape::mrca(tree)[tip1, tip2]
  x <- k

  N <- length(tree$tip.label)
  while(x != N + 1){
    k <- c(k, parentNode(tree, x))
    x <- tail(k, n = 1)
  }
  return(k)
}

# after diverge
v.sum2 <- function(tree, tip, named_alpha, mcra_node){
  nodes <- nodesAlongLineage(tree, tip)
  index <- which(nodes == mcra_node)
  nodes <- head(nodes, n = index-1) # retain the nodes from mcra to tip
  if (length(nodes) == 0){
    return(sum2 = 0)
  }
  else{
    edges <- which(tree$edge[,2] %in% nodes) # from tip to mcra_node
    subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from tip to root
    
    subedge_lengths <- tibble(time_span = subedge_lengths,
                              alpha = named_alpha[names(subedge_lengths)])
    #edges_common <- which(tree$edge[,2] %in% nodes_common) # from tip to root
    #subedge_lengths_common <- rev(unlist(lapply(edges_common, function(i) tree$maps[[i]])))
    
    
    sum2 <- subedge_lengths %>% 
      mutate(sum2 = time_span * alpha) %>% 
      reframe(sum = sum(sum2)) %>% 
      unlist() %>% 
      unname()
    
    return(sum2)
  }
}

#before diverge
v.sum1 <- function(tree, tip1, tip2, named_alpha, named_sigma2){
  nodes <- nodesBeforeDiverge(tree, tip1, tip2)
  if (length(nodes) == 1){
    return(sum1 = 0)
  }
  else{
    edges <- which(tree$edge[,2] %in% nodes) # from tip to root
    subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from mcra_node to root
    
    mcra_time <- sum(subedge_lengths)
    time_point <- cumsum(c(mcra_time, unname(-subedge_lengths))) # from mcra_node to root
    tb <- c(tail(time_point, n = -1)) ## older end of subedge (smaller positive value)
    te <- c(head(time_point, n = -1)) ## younger end of subedge (larger positive value)
    times <- tibble(time_begin = tb,
                    time_end = te,
                    alpha = named_alpha[names(subedge_lengths)],
                    sigma2 = named_sigma2[names(subedge_lengths)])
    sum1 <- times %>% 
      mutate(exp = sigma2 / (2 * alpha) * (exp(-2 * alpha * time_begin) - exp(-2 * alpha * time_end))) %>% 
      reframe(sum1 = sum(exp)) %>% 
      unlist() %>% 
      unname()
    return(sum1)
  }
}

vcv.pairwise <- function(tree, named_alpha, named_sigma2, tip1, tip2){
  mcra_node <- ape::mrca(tree)[tip1, tip2]
  sum2_tip1 <- v.sum2(tree, tip1, named_alpha, mcra_node)
  sum2_tip2 <- v.sum2(tree, tip2, named_alpha, mcra_node)
  
  exp2 <- exp(-1 * (sum2_tip1 + sum2_tip2))
  
  sum1 <- v.sum1(tree, tip1, tip2, named_alpha, named_sigma2)
  
  v_ij <- exp2 * sum1
  
  return(v_ij)
}

vcv.matrix <- function(tree, named_alpha, named_sigma2){
  ntip <- length(tree$tip.label)
  V <- matrix(nrow = ntip, ncol = ntip)
  j = ntip
  while (j != 0){
    for (i in 1:ntip){
      V[i,j] <- vcv.pairwise(tree, named_alpha, named_sigma2, i, j)
      V[j,i] <- V[i,j]
    }
    j = j-1
  }
  colnames(V) <- tree$tip.label
  rownames(V) <- tree$tip.label
  return(V)
}


sd_logL_vcv <- function(tree, continuousChar, named_alpha, named_sigma2, named_theta){
  ntip <- length(tree$tip.label)
  V = vcv.matrix(tree, named_alpha, named_sigma2)
  
  W = weight.matrix(tree, named_alpha)
  
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
  r = solve(L) %*% y - solve(L) %*% W %*% named_theta # what does de-correlated residuals mean?
  
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

# From slouch
data("artiodactyla")
data("neocortex")
neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]

# From local
dummy_tree <- read.simmap("data/1_validation/dummy_threetaxon_simmap.tre",
                          format="phylip",version=1)

# read.nexus.data() gives lists of species with character values, need to change to input format for our functions
char <- read.nexus.data("data/1_validation/dummy_threetaxon_Continuous.nex")
continuousChar <- c()
for (species in char){
  continuousChar <- append(continuousChar, as.numeric(species))
}
names(continuousChar) <- names(char)


###################################################
#                                                 #
#         Guideline for logL_pruning()            #   
#                                                 #
###################################################
logL_pruning(smaptree, brain, sigma2 = 2, alpha = 1, theta = 6)
###################################################
#                                                 #
#       Guideline for sd_logL_pruning()           #   
#                                                 #
###################################################
sd_logL_pruning(smaptree, brain, sigma2 = named_sigma2, alpha = named_alpha, theta = named_theta)
###################################################
#                                                 #
#             Guideline for logL_vcv()            #   
#                                                 #
###################################################

# INPUT
# tree: a tree
# continuousChar: a vector with names of elements == species name
# alpha, sigma2, theta: scalars

# Example
## tree
tree <- artiodactyla

## continuousChar
brain <- neocortex$brain_mass_g_log_mean
names(brain) <- smaptree$tip.label

# RUN
logL_vcv(tree, brain, named_sigma2[[1]], named_alpha[[1]], named_theta[[1]])
###################################################
#                                                 #
#         Guideline for sd_logL_vcv()             #   
#                                                 #
###################################################

# INPUT
# tree: a mapped tree
# continuousChar: a vector with names of elements == species name
# named_alpha, named_sigma2, named_theta: vectors with names of elements == names of discrete character states

# Example
## tree and dummy discrete character history
diet <- as.character(neocortex$diet)
names(diet) <- neocortex$species
set.seed(123)
smaptree <- make.simmap(artiodactyla, diet)
plot(smaptree)

## continuousChar
brain <- neocortex$brain_mass_g_log_mean
names(brain) <- smaptree$tip.label

## parameters
### option 1: set up parameter values and name it after
discrete_states <- unique(diet)
named_sigma2 <- c(2.1, 2.2, 2.3)
names(named_sigma2) <- discrete_states
# option 2: match parameter values and discrete states in one line
named_alpha <- c("MF" = 1.1, "Gr" = 1.2, "Br" = 1.3)
named_theta <- c("MF" = 6.1, "Gr" = 6.2, "Br" = 6.3)

# RUN
sd_logL_vcv(smaptree, brain, named_alpha, named_sigma2, named_theta)
