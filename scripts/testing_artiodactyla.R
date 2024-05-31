library(ape)
library(KScorrect)
library(pracma)
library(phytools)
library(slouch)
library(tibble)
library(tidyverse)

# write data
parameter_1000 <- tibble("alpha_Br" = rep(0,1000),
                         "alpha_Gr" = rep(0,1000),
                         "alpha_MF" = rep(0,1000),
                         "sigma2_Br" =rep(0,1000),
                         "sigma2_Gr" =rep(0,1000),
                         "sigma2_MF" =rep(0,1000),
                         "theta_Br" = rep(0,1000),
                         "theta_Gr" = rep(0,1000),
                         "theta_MF" = rep(0,1000))

for (i in 1:1000){
  parameter_1000[["alpha_Br"]][i] = rexp(n=1, rate=19.61)
  parameter_1000[["alpha_Gr"]][i] = rexp(n=1, rate=19.61)
  parameter_1000[["alpha_MF"]][i] = rexp(n=1, rate=19.61)
  parameter_1000[["sigma2_Br"]][i] = rlunif(n=1, min=1e-5, max=10)
  parameter_1000[["sigma2_Gr"]][i] = rlunif(n=1, min=1e-5, max=10)
  parameter_1000[["sigma2_MF"]][i] = rlunif(n=1, min=1e-5, max=10)
  parameter_1000[["theta_Br"]][i] = runif(n=1, min=0, max=10)
  parameter_1000[["theta_Br"]][i] = runif(n=1, min=0, max=10)
  parameter_1000[["theta_MF"]][i] = runif(n=1, min=0, max=10)
}
#
write.csv(parameter_1000, "data/1_validation/testing_artiodactyla/ou_parameters_all.csv")
#
# data("artiodactyla")
# data("neocortex")
# neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
# brain <- neocortex$brain_mass_g_log_mean
# names(brain) <- artiodactyla$tip.label
# metadata <- read.csv("data/1_validation/testing_artiodactyla/ou_parameters_all.csv")
# alltrees <- read.simmap("data/1_validation/testing_artiodactyla/artiodactyla_all.tre", format="phylip")
# root_state <- c()


# root_state <- as.numeric(root_state)
# metadata$root_state <- root_state
# metadata <- metadata %>% mutate(X = NULL)
# write.csv(metadata, "data/1_validation/testing_artiodactyla/ou_parameters_all.csv")

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
    
    all_weights = lineage %>% mutate(exp_final = exp1 * exp(sum2)) %>% 
      group_by(state) %>% 
      summarise(weight = sum(exp_final))
    
    for (i in 1:nrow(all_weights)){
      W[, all_weights$state[i]] = all_weights$weight[i]
    }
  } else {
    W[, lineage$state[1]] = 1
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




# read data
data("artiodactyla")
data("neocortex")
neocortex <- neocortex[match(artiodactyla$tip.label, neocortex$species), ]
tree <- artiodactyla
brain <- neocortex$brain_mass_g_log_mean
names(brain) <- tree$tip.label



alltrees <- read.simmap("data/1_validation/testing_artiodactyla/artiodactyla.trees", format="phylip",version=1)
metadata <- read.csv("data/1_validation/testing_artiodactyla/ou_parameters_all.csv")
logL <- tibble(rb = rep(0,1000), R_pruning = rep(0,1000), R_vcv = rep(0,1000))

bar = txtProgressBar(style=3, width=40)
for (i in 181:1000){
  tree <- alltrees[[i]]
  #root_state[i] <- names(tree$maps[[1]][1])
  alpha <- c()
  sigma2 <- c()
  theta <- c()
  alpha[1] <- metadata$alpha_Br[i]
  alpha[2] <- metadata$alpha_Gr[i]
  alpha[3] <- metadata$alpha_MF[i]
  sigma2[1] <- metadata$sigma2_Br[i]
  sigma2[2] <- metadata$sigma2_Gr[i]
  sigma2[3] <- metadata$sigma2_MF[i]
  theta[1] <- metadata$theta_Br[i]
  theta[2] <- metadata$theta_Gr[i]
  theta[3] <- metadata$theta_MF[i]
  names(alpha) <- c("0", "1", "2")
  names(sigma2) <- c("0", "1", "2")
  names(theta) <- c("0", "1", "2")
  logL$R_pruning[i] <- sd_logL_pruning(tree, brain, alpha, sigma2, theta)
  if (logL$R_pruning[i] < -1500) next
  if (logL$R_pruning[i] > -800) next
  logL$R_vcv[i] <- sd_logL_vcv(tree, brain, alpha, sigma2, theta)
  setTxtProgressBar(bar, i / 1000)
}



logL$rb <- c(-2468.309, -24785.64, -71.82691, -15329.81, -58661.87, -2123543, -42495.89, -7302.432, -172.1363, -4696.353, -28648.66, -155.0363, -158.352, -9411.262, -194.0261, -43540.71, -1547.616, -7989.36, -43778.86, -82.57928, -34823.08, -113.0799, -213.6456, -402.4085, -190.7875, -1754.08, -6131.339, -9361.487, -19829.72, -952.6882, -76.28786, -109.6679, -98.17324, -20241.6, -2423.318, -148.1917, -254.4509, -252.64, -749.2503, -671.281, -7444.577, -72.2839, -1888.717, -3511.308, -139.1061, -139750.8, -176.0028, -454.3128, -1276.297, -5604.959, -123209.8, -129.4021, -1954.423, -10568.25, -1543.898, -93.74094, -5006.088, -34076.87, -179.703, -5487.865, -122.2083, -97.86571, -4127.089, -4752.476, -101.1673, -785.7062, -700.3579, -2155.765, -280.4455, -200.2608, -111.8447, -8300.504, -322.3717, -70521.62, -1990.863, -20636.85, -835.2134, -6248.922, -274.8449, -136.1804, -1197.375, -68604.43, -1236.02, -244.9317, -86.00942, -86.44113, -1608.097, -11862.46, -6848.522, -147658.1, -108.6507, -108.7744, -18468.82, -480.5396, -994.2138, -964.259, -2701.651, -2712.28, -1124.771, -113.1102, -186.6221, -6219.034, -10201.15, -202340.7, -9197.522, -27653.79, -676.9042, -97.61144, -131380.8, -1703.002, -1204.223, -498.8149, -87.66188, -778.3148, -1729.456, -823.1209, -5758.08, -70.19222, -25073.45, -10091.17, -281.8471, -3576.776, -2430.324, -694.294, -326.6676, -259.7692, -245.2008, -184.8155, -2283.345, -46733.02, -49852.48, -1108.954, -4317.342, -1113.521, -236896.2, -109889.5, -417531.7, -2908.726, -73.60076, -13433.46, -25267.11, -74.16325, -77.93937, -93.87129, -21386.4, -48681.1, -1519.626, -3586.94, -22215.13, -1662.686, -245.1733, -1387.944, -8005.337, -2909.926, -2192.3, -2220.578, -658592.4, -184020.3, -83.3847, -20781.9, -259560.7, -84.89376, -578.7574, -414531.2, -10316.56, -347.6206, -10181.84, -46545.62, -297.9335, -863.1451, -71.66552, -4788.858, -91.04848, -9227.65, -139.6391, -92056.5, -194.641, -357.216, -608.8089, -765.6133, -907.2016, -2624.104, -5095.292, -1537.053, -675.0415, -1126.227, -38173.72, -1045.898, -67.26265, -909.823, -273.2337, -10791.65, -14020.52, -6048940, -2360.61, -1225.192, -1707.42, -575.7298, -139.3867, -6023.979, -593.1978, -117.2578, -2659.998, -5470.951, -231.9196, -3217.241, -8085.772, -669.0042, -8165.008, -215.7044, -20836.1, -82.80844, -828.3606, -219.5671, -101.1537, -556.6609, -186.1744, -59.24563, -111937.9, -1195.213, -1811.133, -254196.7, -47021.27, -88.04631, -140.6299, -1007.931, -15314.61, -20910.57, -123.2814, -1021.145, -264.8904, -65.33322, -8190.397, -1698.571, -67.80352, -89.60548, -82124.73, -1106.526, -93.65357, -139.0982, -123.569, -99.92973, -45546.63, -12872.13, -791.3115, -498969.3, -1222.286, -974.6344, -90.32058, -113.386, -6053.239, -1964.064, -3180.723, -2231.125, -178.5594, -253.9954, -6791.264, -10550.34, -20203.38, -122260.4, -1310.466, -6744.235, -1054.657, -36151.57, -240.8374, -1220.792, -2215.422, -82897.31, -3367.081, -143119, -24985.48, -939.8726, -104.1299, -103.6559, -77212.29, -1504.463, -4181.103, -1709.592, -4571.835, -343.46, -628.4964, -2970.583, -19778.65, -17605.6, -105.984, -564.7477, -471.2721, -22442.44, -1359.008, -94.26098, -999.3083, -6968.794, -3309.319, -392.7962, -266.8367, -7290.38, -168.6918, -6931.837, -165439.6, -298034.4, -4068.849, -3129.433, -76.08835, -1975.972, -23788.72, -169.9826, -193.8486, -110.5678, -2441.189, -59322.03, -7638.382, -131.236, -1341.379, -90.29671, -17589.57, -37376.42, -1303.987, -776.4087, -842.5263, -25430.42, -166151.2, -57339.91, -111.7103, -101.7984, -5187.845, -194885.2, -22492.63, -159.6764, -56207.82, -495.6891, -1323.163, -2315.059, -827.1985, -4561.37, -4383.214, -95.13464, -90410.02, -6412.101, -90.08965, -845.381, -21874.13, -150.0889, -116.4369, -2236.859, -14291.96, -20101.77, -5347.099, -19637.3, -1183.001, -3747.325, -594592.3, -12362.06, -387300.8, -19480.52, -171.064, -271.0721, -139.1004, -16650.7, -377.4222, -4522.904, -43368.85, -75.72753, -229.0154, -11612.07, -12659.23, -558.9058, -51.90345, -52693.17, -1372.567, -2799.788, -76.10589, -79.5465, -99.91038, -231.3285, -4380.447, -220.0571, -9890.169, -99.76831, -1023.894, -2114.68, -79.22842, -99.53013, -351024.9, -82.26634, -367896.7, -877.5043, -269.0318, -7448.836, -385.4031, -1851.051, -1316.278, -43996.97, -217.3361, -5647.709, -544.1557, -193.7115, -2459.585, -7968.056, -997.3561, -5367.899, -26144.77, -1395.61, -707.476, -13706.13, -24665.98, -289379.7, -1900.299, -935.7048, -26241.13, -102.4442, -112.6146, -62.18937, -114.7662, -338159.5, -20384.75, -94.95013, -111.5499, -104.3752, -98.67674, -4838.26, -12961.35, -76.13244, -113190, -392.7882, -27314.2, -90.5363, -13094.03, -13182.4, -1540.571, -1294.894, -4447.888, -115642.4, -8395.352, -508.5636, -122805.3, -92.67694, -591.1157, -135.7391, -87.18301, -5193.167, -1596.367, -466.5727, -151313.2, -102.1246, -1519.999, -2092.368, -13086.48, -14568.78, -81507.61, -5082.133, -5471.129, -1934.782, -951.373, -326.2728, -392.734, -79.31406, -19235.58, -3343.736, -230.0774, -7300.938, -5053.941, -272445, -75993.06, -68.58372, -414.033, -4202.581, -20046.88, -64.60776, -34292.75, -36122.34, -80.76973, -3240.164, -74690.25, -10956.32, -1000.084, -187.9762, -42387.83, -1682.852, -361.935, -9447.001, -95.55955, -7765.663, -7867.183, -473.6806, -602.3167, -263.3384, -3504.329, -339.8025, -30474.05, -98.99632, -1453.583, -1585.568, -93224.67, -327378.1, -1168.371, -77.1681, -16972.04, -74.20005, -111.3711, -343.9008, -12337.94, -5684.732, -12547.16, -565.6885, -9428.258, -71420.29, -88309.85, -1239599, -16560.44, -16731.22, -136.466, -570.6465, -7557.46, -70.77728, -411.7413, -1320.759, -758.0833, -77.49672, -44112.68, -131.8243, -365.4208, -266406.3, -202.4806, -532.8785, -435.3524, -78.3466, -82.28127, -580.6738, -9649.277, -874729.7, -34871.16, -92649.65, -18066.17, -192.6317, -27022.82, -128.6621, -165329.5, -43697.51, -85.39797, -3314.774, -52140.66, -1405.005, -116.7125, -125.3974, -117.2878, -30103.86, -701.5208, -742.656, -79.65396, -4411.855, -183238.3, -2189.142, -187.0501, -5503.695, -56749.21, -12719.08, -121.5813, -76.14091, -3324.813, -7403.094, -499.4629, -83.64997, -251.7048, -248.1994, -76.44276, -207.5242, -141708.6, -4878.274, -555.8122, -107.5794, -126431, -156.0855, -669.0742, -19952.26, -7654.829, -64.18408, -351.4218, -81652, -109.2368, -15442.21, -808.372, -12865.84, -8710.168, -110.5007, -122023.5, -336.6635, -8203.527, -58380.68, -156.2505, -2338.398, -107.5846, -5009.466, -171853.6, -93.54497, -194.8298, -102990.6, -1679.328, -775938.8, -142.5444, -6951.069, -279420.6, -1660.914, -67.99879, -863.3858, -87.42688, -374.3599, -2173.251, -768.1221, -21159.61, -74.84289, -387.8638, -8022.294, -498.477, -2237.988, -1062.301, -3869.645, -535.1727, -97930.75, -5021.036, -1331467, -669.4864, -87.78199, -548.9137, -443.5905, -15936.09, -25494.2, -39852.28, -56.77172, -15019.4, -60.48401, -914.2256, -1184.175, -121835.4, -1699.207, -4538.229, -15191.79, -85302.39, -58.49719, -3338.33, -7081.995, -6922.418, -22703.69, -20125.91, -97.12972, -2474.878, -5968.256, -19158.44, -19289.69, -64957.87, -9172.522, -74209.36, -476.098, -34407.15, -83973.35, -1863315, -95.53453, -20436.53, -66.68353, -1377.57, -261.4095, -434449.2, -1120.396, -717.3585, -2910.08, -251.7864, -1815.391, -23370.38, -319.5003, -61.67435, -3028.021, -867250.2, -48540.43, -17939.74, -3360.287, -5782.058, -7671.41, -49.00358, -2920.31, -1971.249, -114.4329, -28416.2, -192.9177, -85372.3, -136.6364, -19187.28, -92779.37, -1308.804, -34925.67, -68.2023, -621.9975, -2037.096, -607.184, -79650.63, -700.1796, -20286.5, -4475.651, -37418.28, -14618.07, -2265.217, -7157.143, -483.8295, -25612.1, -17652.42, -3803.657, -170.6488, -26032.71, -116.6561, -81.25173, -4226.206, -12557.03, -167.6812, -1380.949, -18031.59, -23651.27, -1078.956, -594.2906, -35429.94, -327.2146, -956.1028, -3505.399, -24319.64, -5529.386, -13895.64, -1065.987, -66017.75, -475.9713, -92.26074, -1089.405, -1050.367, -98.89347, -15749.9, -21839.95, -3956.427, -75.77045, -115770, -384.9312, -135.9743, -2183.129, -226866.8, -34755.31, -190.1122, -1206.107, -115.4162, -86.78504, -784.317, -432.6703, -103.5678, -1427357, -261.8427, -876.5903, -28586.97, -91.14235, -2290.367, -143857.8, -23086.52, -1436.566, -43328.82, -87665.63, -79.01486, -127.0879, -76.46109, -201.3914, -1203.808, -307.7069, -33676.82, -13793.08, -73.01101, -16805.43, -8049.452, -2792.49, -94.27467, -371.5071, -13890.43, -16710.35, -80.55452, -1083.435, -683.6131, -52464.04, -1904.727, -5969.139, -74.43074, -250.5571, -17026.11, -1134.071, -207.0683, -1563.953, -80.55026, -3135.378, -2040.464, -911.2491, -1329.283, -2608.516, -199.9166, -244919.6, -112.7075, -8287.189, -803.2223, -11633.28, -40872.46, -83.86762, -8372.236, -2133.093, -88.3741, -868.3634, -30390.52, -466.5618, -190.8839, -20148.24, -16123.53, -709.2931, -90397.33, -109.1027, -37439.89, -34545.31, -19585.45, -316.6868, -136.2689, -16714.75, -355.5827, -15315.15, -228.5048, -85.29223, -23810.91, -5650.944, -397.3732, -148.4061, -179.772, -2922.255, -8834.61, -89749.59, -5459.09, -217.7853, -1105.074, -107.1693, -183.2307, -109.8995, -678.3385, -1787.374, -30774.33, -236.8801, -2117.742, -54587.08, -11958, -120.1034, -989.6925, -242.039, -497.2269, -403.7786, -365.9512, -284.6871, -763.1387, -88.4524, -1638.476, -11899.25, -81.48788, -5425.348, -25692.94, -61530.44, -13327.23, -336936.7, -100.2651, -3131.402, -1715.068, -1543.072, -92.68563, -491.5025, -119460.1, -1218857, -205.1812, -10202.2, -86.65182, -43117.54, -19014.05, -22526.69, -1316.426, -1614.11, -42418.9, -174485.9, -20508.16, -21855.14, -1260.887, -132.1892, -4421.233, -60888.33, -53218.19, -12389.38, -2567.597, -378.1138, -101.6393, -6718.996, -695.4572, -4329.147, -351.6469, -22589.5, -21204.37, -190720, -194397, -394.4499, -93.76421, -94366.25, -20188.58, -98.036, -244.7626, -45626.44, -1758.46, -41504.26, -4326.759, -76.11329, -6017.42, -758.0301, -91.4672, -2481.053, -5265.768, -197.5942, -431.6644, -110.2835, -1636.495, -130.5628, -96.05346, -81582.01, -2277.747, -571788.5, -11509.96, -14041.21, -5692.134, -2442.867, -1036.794, -681.9481, -54248.4, -18341.1, -71.07288, -45569.58, -103.1367, -179.1236, -1404.083, -162.7678, -208.777, -867.0285, -35937.55, -46371.89, -420533.9, -59896.12, -192479.1, -77084.66, -766.3161, -97.27417, -86996.19, -8936.262, -22352.23, -3279.133, -781.5249, -91.0711, -60394.22, -163.0844, -1822.133, -9881.913, -165575.1, -1079.097, -785.0172, -239.6403, -3247.461, -99.49352, -314.9271, -186.5805, -32492.2, -133.0312, -406924.6, -230.9771, -3955.492, -1036.511, -4509.952, -161309.4, -25721.11, -77.70173, -76952.21, -80839.13, -114.2025, -11654.48, -2391.084, -19614.84, -419.8353, -979.5218, -310.127, -6486.243, -10866.73, -9400.825, -564581.9, -13896.06, -77.77657, -5926.53, -511.5995, -148.9691, -95.13065, -84133.18, -143666.8, -186.098, -634.202, -168.2896, -7535.401, -7646.822)

save(logL, file="output/likelihood_comparison_across_methods.Rda")




