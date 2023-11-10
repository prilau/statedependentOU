library(ape)
library(phytools)
library(geiger)
library(TESS)


treeheight <- function(tree) max(node.depth.edgelength(tree))
obtainContinuousStates_ver7 = function(tree, halflife_0,
                                       theta_0,
                                       stationaryvar_0,
                                       initialValue = theta_0,
                                       dt) {
  if (missing(dt)){
    dt <- 0.002 * treeheight(history)
  }
  
  ## Re-parameterization
  alpha_0 <- log(2) / halflife_0 
  sigma_0 <- sqrt(stationaryvar_0 * 2 * alpha_0)

  cont_states <- list()
  preorder <- rev(postorder(tree))
  edges <- tree$edge
  ntips <- length(tree$tip.label)
  root_node <- ntips + 1
  node_values <- list()
  node_values[[root_node]] <- initialValue
  
  #for (i in 1:length(branch_order)) {
  for (edge_index in preorder){
    parent_node <- edges[edge_index, 1]
    xt0 <- node_values[[parent_node]]
    dt_length <- tree$edge.length[[edge_index]] %/% dt
    dt_remainder = tree$edge.length[[edge_index]] %% dt
    for (k in 1:dt_length) {
      xt1 <- xt0 + alpha_0 * (theta_0 - xt0) * dt + sigma_0 * sqrt(dt) * rnorm(1)
      xt0 <- xt1
    }
    xt1 <- xt0 + alpha_0 * (theta_0 - xt0) * dt_remainder + sigma_0 * sqrt(dt_remainder) * rnorm(1)
    xt0 <- xt1
    desc_node <- edges[edge_index, 2]
    node_values[[desc_node]] <- xt0
  }
  for (i in 1:length(edges[,2])) {
    is_tip <- !(edges[i,2] %in% edges[,1])
    if (is_tip) {
      tip_label <- tips(tree, edges[i,2])
      cont_states[[tip_label]] <- unname(node_values[[edges[i,2]]])
    }
  }
  res <- list(
    cont_states
  )
  return(res)
}

plot(history)




cat("simulating continuous characters.\n")
for(i in 1:4) {
  
  # read the history
  history <- read.tree("data/simulation_1_validation/n1600/t1/tree.tre")
  
  cont_states_ver7 <- obtainContinuousStates_ver7(tree = history,
                                                  halflife_0 = 0.35,
                                                  theta_0 = 0.5,
                                                  stationaryvar_0 = 0.0625,
                                                  initialValue = 0.5, dt = 0.001)
  
  write.nexus.data(cont_states_ver7[[1]], 
                   file = paste0("data/simulation_1_validation/n1600/t1/n1600c", i, "_Continuous.nex"),
                   format="Continuous")
}