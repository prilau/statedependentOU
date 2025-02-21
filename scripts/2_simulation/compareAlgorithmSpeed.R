library(ape)
library(phytools)
library(TESS)
library(pracma)
library(slouch)
library(tibble) 
#library(tidyverse)
library(ggdist)
#source("scripts/5_miscellaneous/functions.R")

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

# simulate trees of different sizes
num_tips   = c(100, 200, 500, 1000, 5000, 10000, 100000)
#num_tips = 100000
reps       = 10

grid = expand.grid(num_tips=num_tips, tree=1:reps, lik_rb=NA, time_rb=NA,
                   lik_pruning=NA, time_pruning=NA, lik_vcv=NA, time_vcv=NA, 
                   stringsAsFactors=FALSE)

# simulate the trees
#bar = txtProgressBar(style=3, width=40)
#for(i in 1:nrow(grid)) {
#  
#  this_row = grid[i,]
#  this_num_tip      = this_row[[1]]
#  this_tree       = this_row[[2]]
#  
#  # create the directories if necessary
#  this_dir = paste0("data/2_simulation/algorithm_speed/n", this_num_tip)
#  if ( !dir.exists(this_dir) ) {
#    dir.create(this_dir, recursive=TRUE, showWarnings=FALSE)
#  }
#  
#  # simulate the tree
#  tree = ladderize(tess.sim.taxa(1, this_num_tip, 10, 1, 0.5)[[1]])
#  
#  # rescale the tree
#  tree$edge.length = tree$edge.length / max(branching.times(tree))
#  
#  # write the tree
#  write.tree(tree, file=paste0(this_dir,"/t", this_tree,".tre"))
#  
#  # increment the progress bar
#  setTxtProgressBar(bar, i / nrow(grid))
#  
#}

# I use fixed OU parameters since the speed of likelihood calculation is not dependent on parameter values
alpha <- c()
sigma2 <- c()
theta <- c()
alpha[1] <- 1
alpha[2] <- 0.1
sigma2[1] <- 4
sigma2[2] <- 8
theta[1] <- 5
theta[2] <- -5
names(alpha) <- names(sigma2) <- names(theta) <- c("0", "1")


# simulate 1 discrete regime history per tree


#Q = matrix(1, 2, 2)
#diag(Q) = -1
#rownames(Q) = colnames(Q) = 0:1
#
#bar = txtProgressBar(style=3, width=40)
#for(i in 1:nrow(grid)) {
#  #if(grid[i,1]!=1e+05) next
#
#  this_row = grid[i,]
#  this_num_tip      = this_row[[1]]
#  this_tree       = this_row[[2]]
#  
#  this_dir = paste0("data/2_simulation/algorithm_speed/n", this_num_tip)
#  if(paste0("t", this_tree, "_simmap.tre") %in% list.files(this_dir)) next
#
#  tree <- read.tree(paste0(this_dir, "/t", this_tree, ".tre"))
#  cat("successfully read tree.\n")
#  
#  tree_length = sum(tree$edge.length)
#  rate = 100 / tree_length
#  
#  set.seed(1503)
#  cat("simulating regime history.\n")
#  history = sim.history(tree, rate * Q, nsim=1, message=FALSE)
#  write.simmap(history, file = paste0(this_dir, "/t", this_tree, "_simmap.tre"))
#  setTxtProgressBar(bar, i / nrow(grid))
#}

# I use same continuous trait for tree replicates with same num_tips since the speed of likelihood calculation is not dependent on the trait values. 
for(i in num_tips) {
  this_dir = paste0("data/2_simulation/algorithm_speed/n", i)
  set.seed(1503)
  cont <- rnorm(i, mean=0, sd=4)
  cont_list <- list()
  for (j in 1:i){
    tip <- paste0("t", j)
    cont_list[[tip]] <- cont[j]
  }
  write.nexus.data(cont_list, paste0(this_dir, "/Continuous.nex"), format = "continuous")
}


#timetaken <- tibble(rb=NA, rPrune=NA, rVcv=NA)

bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  this_row = grid[i,]
  this_num_tip = this_row[[1]]
  this_tree = this_row[[2]]
  this_dir = paste0("data/2_simulation/algorithm_speed/n", this_num_tip)
  history <- read.simmap(paste0(this_dir, "/t", this_tree, "_simmap.tre"), format="phylip")
  
  set.seed(1503)
  cont <- rnorm(this_num_tip, mean=0, sd=4)
  names(cont) <- paste0("t", 1:this_num_tip)
  
  cat("running pruning algorithm.\n")
  start.time <- Sys.time()
  grid[i,5] <- sd_logL_pruning(history, cont, alpha, sigma2, theta)
  end.time <- Sys.time()
  grid[i,6] <- end.time - start.time
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
}

write.csv(grid, file="output/2_simulation/algorithm_speed/r_pruning.csv")

bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  this_row = grid[i,]
  this_num_tip = this_row[[1]]
  this_tree = this_row[[2]]
  this_dir = paste0("data/2_simulation/algorithm_speed/n", this_num_tip)
  history <- read.simmap(paste0(this_dir, "/t", this_tree, "_simmap.tre"), format="phylip")
  
  set.seed(1503)
  cont <- rnorm(this_num_tip, mean=0, sd=4)
  names(cont) <- paste0("t", 1:this_num_tip)
  
  cat("running vcv algorithm.\n")
  start.time <- Sys.time()
  grid[i,7] <- sd_logL_vcv(history, cont, alpha, sigma2, theta)
  end.time <- Sys.time()
  grid[i,8] <- end.time - start.time
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
}

write.csv(grid, file="output/2_simulation/algorithm_speed/r_vcv.csv")


#############
## plotting #
#############
#
## Calculate intercept and slope of regression line
#grid_tmp <- grid %>% filter(!is.na(time_pruning)) %>% 
#  mutate(tree_size=log(num_tips),
#         ln_sec_pruning = log(time_pruning * 60)) %>%
#  select(tree_size, ln_sec_pruning)
#coef(lm(grid_tmp$ln_sec_pruning ~ grid_tmp$tree_size))
#
#
#
#grid %>% pivot_longer(cols=c(10,12), names_to = "method",
#                      values_to = "time") %>% 
#  ggplot(grid) +
#  geom_boxplot(aes(y=time, group=c(num_tips, method)))
#
#grid %>%
#  filter(!is.na(time_pruning)) %>% 
#  mutate(tree_size=log(num_tips),
#         sec_pruning = time_pruning * 60) %>% 
#  ggplot() +
#  #ggdist::stat_pointinterval(aes(x = tree_size, y = log(sec_pruning),
#  #                               group = tree_size),
#  #                           size=10) +
#  #  ggdist::stat_dots(aes(x = tree_size, y = log(sec_pruning),
#  #                        group = tree_size, color=as.factor(tree_size), fill=as.factor(tree_size)),
#  #    ## orientation to the left
#  #    side = "left", 
#  #    ## move geom to the left
#  #    justification = 1.2, 
#  #    ## adjust grouping (binning) of observations 
#  #    binwidth = 0.0001,
#  #    alpha=0.3,
#  #    dotsize=1000
#  #  ) +
#  geom_point(aes(x=tree_size, y=log(sec_pruning), group=tree_size, color=as.factor(tree_size)),
#             #position=position_jitter(width = 0.15),
#             size=1.2, alpha=0.4) +
#  geom_abline(linetype="dashed", alpha=0.5,
#              intercept=-9.102700,
#              slope=1.548436 
#  ) +
#  scale_color_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
#  scale_fill_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
#  scale_x_continuous(breaks=log(c(100, 200, 500, 1000, 5000, 10000)), label=c("100", "200", "500", "1k", "5k", "10k")) +
#  scale_y_continuous(breaks=log(c(0.25, 0.5, 1, 5, 0, 30, 60, 180)), label=c("0.25s", "0.5s", "1s", "5s", "10s", "30s", "1min", "3min")) +
#  ylab("Time taken for one iteration") +
#  theme_classic() +
#  ggtitle("Likelihood calculation using pruning algorithm") +
#  theme(plot.title = element_text(hjust=0.5)) +
#  xlab("Tree size") +
#  theme(legend.position="none")

#grid %>%
#  filter(!is.na(time_pruning)) %>% 
#  mutate(tree_size=log(num_tips),
#         sec_pruning = time_pruning * 60) %>% 
#  ggplot(aes(x = tree_size, y = log(sec_pruning), group = tree_size,
#             color=as.factor(tree_size))) + 
#  ## add half-violin from {ggdist} package
#  geom_boxplot(
#    width = .2, 
#    linewidth=0.2,
#    ## remove outliers
#    outlier.color = NA ## `outlier.shape = NA` or `outlier.alpha = 0` works as well
#  ) +
#  ## add dot plots from {ggdist} package
#  ggdist::stat_dots(
#    aes(fill=as.factor(tree_size)),
#    ## orientation to the left
#    side = "left", 
#    ## move geom to the left
#    justification = 1.2, 
#    ## adjust grouping (binning) of observations 
#    binwidth = 0.005,
#    alpha=0.3,
#    dotsize=20
#  ) + theme_classic() +
#  scale_x_continuous(breaks=log(c(100, 200, 500, 1000, 5000, 10000)), label=c("100", "200", "500", "1k", "5k", "10k")) +
#  scale_y_continuous(breaks=log(c(0.25, 0.5, 1, 5, 0, 30, 60, 180)), label=c("0.25s", "0.5s", "1s", "5s", "10s", "30s", "1min", "3min")) +
#  ylab("Time taken for one iteration") +
#  ggtitle("Likelihood calculation using pruning algorithm") +
#  theme(plot.title = element_text(hjust=0.5)) +
#  xlab("Tree size") +
#  scale_color_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
#  scale_fill_manual(values=c('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377')) +
#  theme(legend.position="none")
  
  
  
  ## remove white space on the sides
  #coord_cartesian(xlim = c(1.3, 2.9))
