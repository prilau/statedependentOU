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
for (i in 1:1000){
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
  logL$R_vcv[i] <- sd_logL_vcv(tree, brain, alpha, sigma2, theta)
  delta = log10(abs(logL$R_vcv[i] - logL$R_pruning[i]))
  cat(delta)
  setTxtProgressBar(bar, i / nrow(grid))
}


rb <- c(-72.6084, -115.7578, -110.6516, -76.95967, -1429.058, -503.656, -208.5575, -1562.999, -120.6112, -2229.951, -213.3462, -58.41681, -1009.345, -82.87116, -88.70166, -68.77792, -359.4803, -1972.743, -292.9592, -809.7714, -490.1009, -85.06839, -256.1854, -279.5419, -5393.097, -169.0738, -1138.529, -3188.177, -164.2645, -165.5264, -96.95484, -1167.898, -3002.669, -312.2576, -760.9732, -1104.367, -147.514, -249.922, -533.8657, -233.3733, -1082.198, -630.2748, -78.95238, -890.5392, -113.148, -411.866, -195.99, -1038.725, -165.3171, -543.2153, -389.3023, -91.00754, -165.2065, -163.1596, -136.4967, -1053.435, -1458.746, -118.0226, -244.4299, -399.6847, -300.1411, -173.9038, -453.937, -692.6883, -79.55002, -92.33177, -145.3483, -109.0905, -1393.902, -953.8958, -640.2978, -98.96109, -198.5391, -52.1253, -253.8197, -665.859, -2270.801, -307.0706, -404.8608, -57.83114, -196.4837, -754.426, -177.3719, -116.0304, -719.9016, -251.4545, -88.01059, -176.4526, -189.4593, -1110.627, -120.377, -3129.735, -7303.782, -104.5744, -12419.38, -76.35507, -1054.975, -136.0114, -501.0487, -244.0038, -80.08015, -974.2094, -229.021, -376.667, -1924.743, -301.9298, -64.17509, -220.4051, -85.16314, -242.8656, -62.40016, -259.1913, -1213.011, -504.2052, -425.7589, -521.8702, -184.7197, -704.3817, -686.4647, -948.1779, -83.26507, -6852.886, -142.122, -15105.53, -614.201, -121.4475, -319.5397, -88.23793, -1518.323, -1302.267, -103.0897, -1139.53, -735.8119, -1528.292, -183.3851, -84.71453, -69.49504, -95.90522, -481.35, -201.4511, -2673.279, -1314.797, -372.2157, -1100.619, -104.4755, -2146.022, -722.6698, -3149.955, -196.8765, -1467.163, -576.2565, -662.8426, -633.8212, -1295.749, -650.3496, -266.1061, -257.4996, -117.6567, -636.7839, -91.57611, -122.5528, -3987.102, -661.6744, -77.74248, -509.2256, -1551.901, -693.4916, -379.0426, -431.6432, -230.9031, -113.641, -79.83114, -154.3835, -1219.981, -205.0331, -300.2428, -1439.811, -133.3234, -111.0483, -282.8639, -83.51152, -121.6148, -291.267, -163.0739, -789.5075, -487.8332, -2693.179, -2160.412, -381.4345, -478.6083, -1073.537, -147.7478, -147.5896, -78.65341, -823.7121, -349.9383, -274.3976, -1387.8, -1775.937, -467.5556, -226.6769, -38712.36, -996.4071, -24463.26, -142.6405, -499.9702, -153.8974, -211.0895, -380.6434, -225.6757, -505.2023, -81.22607, -614.1131, -520.4537, -226.4177, -1878.644, -575.1909, -406.4115, -138.3111, -2109.98, -2570.728, -219.4211, -149.549, -1515.67, -295.5001, -120.8119, -8982.949, -88.71524, -784.1032, -967.3047, -80.79116, -108.0172, -291.5954, -299.5143, -893.4025, -2633.575, -256.3647, -1563.662, -298.7354, -273.7816, -136.3288, -528.8955, -2208.058, -106.6761, -268.6226, -291.1955, -98.40984, -606.1053, -375.9234, -891.1528, -8431.133, -242.3761, -516.8657, -1001.997, -296.2329, -61.8671, -605.2926, -113.5323, -282.0316, -47.80026, -5468.696, -463.7983, -604.5281, -290.6094, -127.5883, -95.82287, -626.3233, -435.2305, -92.18613, -395.5591, -887.9665, -433.5466, -2626.55, -647.3487, -353.658, -111.1401, -2286.883, -967.5471, -94.87544, -250.5969, -1269.302, -125.4254, -254.614, -467.2063, -325.5747, -1459.333, -793.5226, -90.71461, -226.8787, -230.3609, -319.1255, -170.5312, -901.2164, -3497.891, -91.67427, -348.7308, -1071.482, -554.1172, -1188.142, -585.0328, -85.39332, -138.6284, -64.28428, -134.6149, -250.2669, -278.6646, -170.3969, -79.03994, -781.5455, -1107.617, -372.2599, -85.02427, -402.7739, -322.4952, -365.1317, -762.4887, -254.9658, -12003.44, -2375.114, -80.88069, -1435.383, -74.71699, -651.5182, -189.5931, -167.6364, -211.2516, -314.4362, -1395.313, -281.7056, -213.1261, -559.8139, -585.6671, -199.4711, -461.9483, -975.64, -458.2564, -278.316, -563.6219, -1844.374, -287.8289, -62.75711, -550.9088, -486.0746, -846.0885, -489.9724, -362.8985, -92.58692, -310.2197, -3628.61, -104.4414, -192.5205, -1190.743, -55.39641, -62.67156, -238.9486, -642.6084, -1664.482, -2478.245, -355.0641, -187.1507, -3623.958, -849.7552, -215.2667, -707.6404, -101.1395, -1367.337, -112.0163, -1778.921, -1186.609, -11653.85, -459.9344, -100.6277, -193.9926, -391.5045, -107.568, -167.4799, -1081.353, -358.8642, -136.4074, -134.5665, -2077.592, -859.3716, -252.3003, -858.5004, -396.0648, -1371.412, -590.8836, -103.9659, -509.3442, -799.3551, -122.0305, -3523.586, -136.7719, -411.4929, -70.37486, -70.16089, -3127.276, -120.3226, -2440.124, -148.5661, -2524.479, -1923.271, -547.7499, -303.4869, -120.2148, -91.68679, -2003.861, -128.9044, -527.93, -520.4251, -889.1725, -395.8041, -1913.987, -1606.636, -1515.534, -150.5786, -204.7869, -524.0419, -479.0508, -721.702, -275.358, -125.8917, -308.7593, -340.1085, -321.1537, -129.1094, -285.4776, -125.4655, -406.802, -417.917, -1244.409, -1596.872, -2254.815, -104.274, -205.6163, -3360.199, -866.204, -180.2049, -208.5065, -108.2573, -487.4651, -950.6609, -368.379, -178.6288, -1987.805, -56.73706, -199.0143, -88.879, -1860.097, -1480.69, -101.4169, -65.22493, -2576.902, -86.91894, -216.1693, -382.605, -328.3156, -206.9194, -78.00559, -119.9794, -278.0624, -441.8648, -272.9098, -249.8891, -596.538, -1178.478, -105.893, -90.29012, -69.78515, -80.70354, -127.2879, -6536.948, -1990.223, -173.8989, -485.1115, -635.7858, -90.23209, -55.03521, -290.557, -1974.277, -262.1233, -102.1257, -128.3929, -277.9724, -109.1693, -1864.928, -2577.397, -81.61406, -177.0018, -2289.865, -157.939, -201.3856, -152.2712, -211.6862, -7628.667, -286.0267, -155.6391, -283.5129, -471.097, -608.598, -244.2433, -1925.915, -261.2833, -110.2581, -262.5933, -424.1004, -133.6782, -239.5821, -305.7592, -1206.627, -2019.297, -194.5939, -124.0444, -125.3262, -447.2799, -168.7282, -1044.18, -123.9101, -179.9481, -1651.455, -107.4459, -49.48312, -86.08222, -133.5139, -251.2532, -231.581, -166.2052, -253.0048, -255.4748, -831.6746, -112.233, -315.5293, -395.3134, -1284.592, -204.4282, -66.22585, -541.4795, -629.8227, -555.6435, -113.8276, -220.3263, -100.7546, -686.1153, -1095.177, -176.5733, -275.6774, -668.678, -72.94058, -288.9827, -793.5472, -17557.89, -134.9291, -741.4648, -935.6007, -462.2662, -51.75288, -681.8188, -913.0113, -61.54388, -1207.316, -3728.455, -104.104, -996.7056, -838.2277, -186.2725, -275.3731, -265.4134, -410.9137, -246.4696, -856.1375, -112.1293, -578.2347, -116.6373, -128.3384, -365.5808, -80.95993, -181.149, -1467.096, -332.4481, -302.785, -297.0789, -395.2671, -5596.201, -79.12301, -103.4401, -548.9491, -261.3308, -83.97941, -1437.994, -3324.996, -1652.922, -310.9529, -92.40432, -252.4747, -258.2836, -66.23457, -169.1102, -711.7949, -1025.496, -900.3612, -2458.395, -372.3501, -722.4776, -131.4971, -269.105, -1000.881, -1689.133, -361.7238, -577.354, -959.7902, -140.4487, -95.33098, -553.4327, -267.9743, -405.0268, -794.4293, -72.79931, -67.3711, -717.0661, -462.8371, -2090.054, -227.5958, -1391.793, -282.1498, -112.7481, -232.0747, -385.8978, -1798.234, -298.7625, -256.8379, -1665.06, -622.4999, -79.44645, -120.3957, -332.1108, -70.18955, -1314.662, -70.23976, -3450.377, -238.5076, -865.0258, -249.4036, -1911.877, -3176.051, -204.3569, -1624.178, -335.9285, -326.7263, -1940.985, -143.5088, -122.5147, -100.7017, -68.79675, -783.6338, -64.69339, -72.35534, -176.2895, -1703.96, -1133.058, -147.2694, -206.217, -140.1455, -178.8262, -408.3117, -481.0835, -156.0995, -331.6632, -83.5937, -141.5341, -64.31167, -84.16823, -55.01532, -616.4504, -248.7804, -56.2827, -208.4407, -497.8503, -300.4507, -212.0126, -6638.687, -206.1212, -124.0641, -103.8091, -184.9815, -301.4083, -352.7355, -461.4986, -157.6176, -1768.745, -78.84064, -71.52378, -58.31214, -242.2961, -460.6399, -4269.637, -155.6005, -231.7428, -875.9894, -421.101, -9620.441, -4199.719, -168.3963, -62.14634, -505.4514, -160.7173, -49.69849, -506.7923, -3766.717, -970.6283, -150.784, -298.9809, -88.39019, -170.9529, -118.6055, -169.2212, -530.4092, -2489.157, -2866.941, -225.5562, -74.85914, -1014.882, -69.1297, -260.6774, -394.5474, -298.0274, -205.6889, -2692.37, -657.3967, -276.8851, -85.92205, -581.5474, -809.3542, -275.5792, -105.8173, -1703.438, -184.623, -616.5315, -69.73673, -346.9783, -142.3676, -81.61594, -65.2959, -256.2883, -81.62579, -103.4641, -326.2043, -103.9219, -118.7736, -245.7786, -359.8995, -954.9234, -492.2972, -91.1018, -83.18666, -136.6819, -119.2989, -880.8135, -98.37718, -1232.469, -1070.137, -245.1617, -597.6486, -568.1739, -121.9879, -1980.293, -1623.228, -844.2115, -1367.707, -542.3496, -326.5735, -460.3797, -296.718, -211.8756, -139.9416, -103.114, -8388.961, -385.2305, -257.7392, -115.7436, -397.4533, -110.1391, -168.3126, -141.8974, -118.7038, -563.7371, -173.471, -1258.892, -369.3086, -851.46, -103.979, -269.3621, -770.826, -307.299, -51.51854, -204.1689, -89.75644, -401.46, -726.9818, -328.1161, -287.1132, -97.99124, -183.4501, -167.0729, -1270.095, -885.2495, -580.9281, -646.5767, -983.4913, -259.6585, -583.1818, -381.0457, -568.6367, -566.8843, -1097.945, -253.2765, -266.0308, -156.4282, -179.6311, -1189.163, -1360.104, -159.1665, -94.1644, -355.1784, -283.9889, -125.6403, -138.3646, -531.6108, -184.1899, -366.6027, -95.12953, -88.42796, -655.526, -1868.554, -978.0002, -749.3163, -548.4953, -1665.685, -74.06027, -502.3353, -64.90022, -713.9455, -135.9738, -1209.077, -519.2485, -82.97283, -69.05714, -232.8566, -74.02649, -87.14702, -88.11667, -483.9038, -1394.981, -233.7866, -96.43266, -1265.281, -145.3686, -67.989, -136.4632, -201.5485, -671.4932, -230.2411, -1856.939, -116.2901, -581.5985, -4556.075, -504.3689, -78.876, -132.9511, -227.3375, -484.0491, -90.52034, -153.9265, -63.55313, -344.2174, -222.332, -1429.627, -87.02938, -939.2005, -2810.588, -1530.535, -269.8682, -135.5592, -185.2329, -1373.603, -146.6717, -730.7011, -78.25983, -669.4654, -158.6627, -309.0842, -185.7215, -117.9195, -1496.672, -347.8573, -100.0521, -203.7242, -178.8271, -112.7003, -141.0474, -176.2554, -345.1437, -121.7695, -434.331, -66.66499, -337.6801, -2978.765, -380.378, -320.4482, -186.2968, -215.2111, -211.4467, -304.8275, -733.1509, -9955.575, -272.8835, -245.5059, -824.3364, -425.3833, -347.025, -233.9808, -374.9493, -245.6659, -3899.338, -666.0091, -140.6756, -117.6499, -520.3422, -715.2694, -713.63, -95.17385, -796.9712, -103.814, -9468.662, -143.053, -153.3667, -1081.789, -302.3421, -203.5861, -67.79594, -317.1495, -335.163, -1253.604, -140.6282, -272.1235, -1509.443, -365.4692, -456.3895, -252.3761, -83.5412, -1209.944, -171.073, -184.7504, -1526.525, -2080.294, -268.018, -2019.855, -482.0634, -114.12, -532.4495, -81.33575, -3152.489, -2599.814, -840.335, -572.8467, -1354.9, -1388.357, -282.5812, -385.5017, -85.89238, -1218.933, -80.60332, -152.6688, -61.04869, -2147.474, -157.6311, -5191.152, -77.54848, -292.5079, -582.2005, -199.8316, -2355.938, -667.6958, -80.78741, -78.60421, -203.1053, -96.07314, -130.2534, -642.5582, -138.3447, -208.7831, -410.4116, -157.5356, -3760.343, -122.3402, -3829.73, -642.9106, -59.17004, -86.42986, -1812.698, -554.674)

logL$rb <- rb

save(logL, file="output/likelihood_comparison_across_methods.Rda")


ggplot(logL) +
  geom_point(aes(x=1:1000,y=sort(rb)), alpha=0.5, shape=1, col="red") +
  geom_point(aes(x=1:1000,y=sort(R_pruning)), alpha=0.5, shape=4) +
  geom_line(aes(x=1:1000,y=sort(R_vcv)), alpha=0.2, linewidth=1.8)


