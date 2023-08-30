library(ape)
library(phytools)
library(geiger)
library(TESS)
source("scripts/readWriteCharacterData.R")

# simulate the tree
num_tips = 50
tree = ladderize(tess.sim.taxa(1, num_tips, 10, 1, 0.5)[[1]])
# rescale the tree
tree$edge.length = tree$edge.length / max(branching.times(tree))

write.tree(tree, file=paste0("data/n50_simulation2.tre"))


tree_lengths <- sum(tree$edge.length)
# specify rates so that the expected number of changes is 5
rates = 8 / tree_lengths
names(rates) = num_tips

# specify the Mk2 rate matrix
Q = matrix(1, 2, 2)
diag(Q) = -1
rownames(Q) = colnames(Q) = 1:2 - 1

# simulate the discrete characters
# track the number of rejected simulations based on proportional
# representation
colors = c("0"="blue","1"="red")
num_rejections = numeric(length(num_tips))
num_simulations = numeric(length(num_tips))
names(num_rejections) = names(num_simulations) = num_tips
this_rate = rates[as.character(num_tips)]
history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
while (! (mean(history$states == "0") > 0.2 & (mean(history$states == "1") > 0.2) ) ) {
  history = sim.history(tree, this_rate * Q, nsim=1, message=FALSE)
  num_rejections[as.character(num_tips)] = num_rejections[as.character(num_tips)] + 1
}
num_simulations[as.character(num_tips)] = num_simulations[as.character(num_tips)] + 1
maps = history$mapped.edge[,c("0","1")]

pdf("data/n50History2.pdf")
plot(history, col=colors)
dev.off()

writeCharacterData(t(t(history$states)), file=paste0("data/n50_simulationDiscrete2.nex"), type="Standard")






# obtain root state
obtainRootState = function(tree) {
  edge1d <- rev(postorder(tree))[1]
  rootState <- names(history$maps[[edge1d]][1])
  rootState <- as.integer(rootState)
  return(rootState)
}

obtainRootState(tree)
## node0 <- basal_node(tree)
## node1d <- Descendants(node0, tree)[1]
## edge1d <- which.edge(node1d, tree)
## names(history$maps[[edge1d]][1]) == "1" or "0"?

#obtainBaseToTipBranches = function(tree) {
#  branches <- data.frame()
#  for (i in 1:length(tree$tip.label)) {
#    currentNode <- tree$tip.label[i]
#    j = 1
#    while (is.na(which.edge(tree, currentNode)) == FALSE) {
#      branches[i,j] <- which.edge(tree, currentNode)
#      currentNode <- Parent(currentNode, tree)
#      j = j+1
#    }
#    rownames(branches)[i] <- tree$tip.label[i]
#  }
#  branches <- rev(branches)
#  branches <- as.data.frame(t(apply(branches, 1, function(x) x[order(is.na(x))])))
#  return(branches)
#}

#branch_names <- obtainBaseToTipBranches(tree)
#
#branch_names
#
#
#obtainBranchLengths = function(branches) {
#  branchLengths <- list()
#  for (i in 1:length(branches)) {
#    branch <- unlist(branches[i])
#    branchLength <- c()
#    for (j in 1:length(branch)) {
#      fragment <- branch[j]
#      # 1/0 or 0/1 is base-to-tip - checked
#      branchLength <- append(branchLength, unlist(history$maps[fragment]))
#    }
#    branchLengths[[names(branches[i])]] <- branchLength
#  }
#  return(branchLengths)
#}
#
#branch_frag_lengths <- obtainBranchLengths(branch_names)
# obtain basal node - basal_node(tree)
# obtain descendents of basal node (2-deg nodes) - Descendants(basal node #, tree)
# obtain edge # of 2-deg nodes - which.edge(tree, 2-deg node #)
# obtain edge lengths of 2-deg nodes' edges - testTree$edge.length[#]
# how to end?





#obtainContinuousStates = function(branchLengths, alphaRoot, alphaAlt, thetaRoot, 
#                                  thetaAlt, sigmaRoot, sigmaAlt, initialState) {
#  contStates <- list()
#  for (i in 1:length(branchLengths)) {
#    branchLength <- unlist(branchLengths[i])
#    xt0 <- initialState
#    j = 1
#    ## the state of fragments can be 1-1 consecutively!!
#    while (j < length(branchLength)) {
#      t01 <- branchLength[j]
#      t12 <- branchLength[j+1]
#      xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * t01 + sigmaRoot * sqrt(t01) * rnorm(1)
#      xt2 <- xt0 + alphaAlt * (thetaAlt - xt0) * t12 + sigmaAlt * sqrt(t12) * rnorm(1)
#      xt0 <- xt2
#      j = j + 2
#    }
#    if (j == length(branchLength)) {
#      t01 <- branchLength[j]
#      xPresent <- xt0 + alphaRoot * (thetaRoot - xt0) * t01 + sigmaRoot * sqrt(t01) * rnorm(1)
#    } 
#    
#    else if (j == length(branchLength) + 1) {
#      t01 <- branchLength[j-2]
#      t12 <- branchLength[j-1]
#      xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * t01 + sigmaRoot * sqrt(t01) * rnorm(1)
#      xPresent <- xt0 + alphaAlt * (thetaAlt - xt0) * t12 + sigmaAlt * sqrt(t12) * rnorm(1)
#    }
#    else {
#      print("sth went wrong!")
#    }
#    contStates[[names(branchLengths[i])]]$trait1 <- xPresent
#  }
#  for (i in 1:length(contStates)) {
#    contStates[i] <- unname(contStates[[i]])
#  }
#  return(contStates)
#}

#cont_states <- obtainContinuousStates(branchLengths = branch_frag_lengths, alphaRoot = 0.5, alphaAlt = 0.5, thetaRoot = 10, 
#                                           thetaAlt = 30, sigmaRoot = 5, sigmaAlt = 5, initialState = 10)
#cont_states



#obtainContinuousStates_ver2 = function(tree, branchLengths, alphaRoot, alphaAlt,
#                                       thetaRoot, thetaAlt, sigmaRoot, sigmaAlt,
#                                       initialState) {
#  contStates <- list()
#  ## obtain root state
#  rootState <- obtainRootState(tree)
#  ##
#  for (i in 1:length(branchLengths)) {
#    branchLength <- branchLengths[[i]]
#    xt0 <- initialState
#    for (j in 1:length(branchLength)) {
#      if (rootState == names(branchLength[j])) {
#        dt <- branchLength[j]
#        xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt + sigmaRoot * sqrt(dt) * rnorm(1)
#        xt0 <- xt1
#      } else {
#        dt <- branchLength[j]
#        xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt + sigmaAlt * sqrt(dt) * rnorm(1)
#        xt0 <- xt1        
#      }
#    }
#    contStates[[names(branchLengths[i])]] <- unname(xt1)
#  }
#  return(contStates)
#}
#
#
#
#
#cont_states_ver2 <- obtainContinuousStates_ver2(tree = tree, branchLengths = branch_frag_lengths, alphaRoot = 2, alphaAlt = 2, thetaRoot = 30, 
#                                                thetaAlt = 10, sigmaRoot = 2, sigmaAlt = 2)
#cont_states_ver2

#obtainContinuousStates_ver3 = function(tree, branchLengths, alphaRoot, alphaAlt,
#                                       thetaRoot, thetaAlt, sigmaRoot, sigmaAlt,
#                                       initialState = thetaRoot, dt = 0.002) {
#  contStates <- list()
#  ## obtain root state
#  rootState <- obtainRootState(tree)
#  ##
#  for (i in 1:length(branchLengths)) {
#    branchLength <- branchLengths[[i]]
#    xt0 <- initialState
#    for (j in 1:length(branchLength)) {
#      dt_length = branchLength[j] %/% dt
#      dt_remainder = branchLength[j] %% dt
#      if (rootState == names(branchLength[j])) {
#        for (k in 1:dt_length) {
#          xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt + sigmaRoot * sqrt(dt) * rnorm(1)
#          xt0 <- xt1
#        }
#        xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt_remainder + sigmaRoot * sqrt(dt_remainder) * rnorm(1)
#        xt0 <- xt1
#      } else {
#        for (k in 1:dt_length) {
#          xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt + sigmaAlt * sqrt(dt) * rnorm(1)
#          xt0 <- xt1
#        }
#        xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt_remainder + sigmaAlt * sqrt(dt_remainder) * rnorm(1)
#        xt0 <- xt1
#      }
#    }
#    contStates[[names(branchLengths[i])]] <- unname(xt1)
#  }
#  return(contStates)
#}
#
#
#
#
#cont_states_ver3 <- obtainContinuousStates_ver3(
#  tree = tree,branchLengths = branch_frag_lengths, alphaRoot = 1, alphaAlt = 1, 
#  thetaRoot = 50, thetaAlt = 20, sigmaRoot = 2, sigmaAlt = 2, initialState = 20)
#
#cont_states_ver3

#rec_fac <- function(x){
#  if(x==0 || x==1)
#  {
#    return(1)
#  }
#  else
#  {
#    return(x*rec_fac(x-1))
#  }
#}


#obtainContinuousStates_ver4 = function(tree, alphaRoot, alphaAlt, thetaRoot,
#                                       thetaAlt, sigmaRoot, sigmaAlt,
#                                       initialState = thetaRoot, dt = 0.002,
#                                       currentNode = NA) {
#  contStates <- list()
#  ## obtain root state
#  rootState <- obtainRootState(tree)
#  
#  ##if it is at root node right now
#  if (is.na(currentNode) == TRUE) {
#    currentNode <- Descendants(basal_node(tree), tree)
#    xt0 <- initialState
#  }
#  for (i in 1: length(currentNode)) {
#    ## if it is at tip node right now
#    if (length(tips(tree, currentNode[i])) == 1) {
#      contStates[[tips(tree, currentNode[i])]] <- initial_state[[currentNode[i]]]
#    }
#    ## if it is any other internal nodes
#    else {
#      edge <- which.edge(tree, currentNode[i])
#      edgeLength <- history$maps[[edge]]
#      parentNode <- Parent(currentNode[i], tree)
#      xt0 <- initial_state[[parentNode]]
#      for (j in 1:length(edgeLength)) {
#        subEdgeLength <- edgeLength[j]
#        dt_length = branchLength[j] %/% dt
#        dt_remainder = branchLength[j] %% dt
#        
#        if (rootState == names(edgeLength[j])) {
#          for (k in 1:dt_length) {
#            xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt + sigmaRoot * sqrt(dt) * rnorm(1)
#            xt0 <- xt1
#          }
#          xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt_remainder + sigmaRoot * sqrt(dt_remainder) * rnorm(1)
#          xt0 <- xt1
#        }
#        else {
#          for (k in 1:dt_length) {
#            xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt + sigmaAlt * sqrt(dt) * rnorm(1)
#            xt0 <- xt1
#          }
#          xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt_remainder + sigmaAlt * sqrt(dt_remainder) * rnorm(1)
#          xt0 <- xt1
#        }
#        initial_state <- xt1
#        next_nodes <- Descendants(currentNode[i], tree)
#        obtainContinuousStates_ver4(tree = tree, alphaRoot = alphaRoot,
#                                    alphaAlt = alphaAlt, thetaRoot = thetaRoot, 
#                                    thetaAlt = thetaAlt, sigmaRoot = sigmaRoot, 
#                                    sigmaAlt = sigmaAlt,
#                                    initialState = initial_state, dt = 0.002,
#                                    currentNode = next_nodes)
#    }
#    }
#  }
#  return(contStates)
#}
#
#
#
#
#cont_states_ver4 <- obtainContinuousStates_ver4(tree = tree, alphaRoot = 1,
#                                                alphaAlt = 1, thetaRoot = 50,
#                                                thetaAlt = 20, sigmaRoot = 2,
#                                                sigmaAlt = 2, initialState = 50,
#                                                currentNode = NA)
#
#cont_states_ver4


#obtainContinuousStates_ver5 = function(tree, alphaRoot, alphaAlt, thetaRoot,
#                                       thetaAlt, sigmaRoot, sigmaAlt,
#                                       initialState = thetaRoot, dt = 0.002,
#                                       currentNode = NULL) {
#  contStates <- list()
#  ## obtain root state
#  rootState <- obtainRootState(tree)
#  
#  ##if it is at root node right now
#  if (is.null(currentNode[i]) == TRUE) {
#    currentNode <- Descendants(basal_node(tree), tree)
#  }
#  xt0 <- initialState
#  for (i in 1: length(currentNode)) {
#    ## if it is at tip node right now
#    if (isTRUE(length(tips(tree, currentNode[i])) == 1)) {
#      contStates[[tips(tree, currentNode[i])]] <- initial_state[[currentNode[i]]]
#    }
#    ## if it is any other internal nodes
#    else {
#      edge <- which.edge(tree, currentNode[i])
#      edgeLength <- history$maps[[edge]]
#      for (j in 1:length(edgeLength)) {
#        subEdgeLength <- edgeLength[[j]]
#        dt_length = subEdgeLength %/% dt
#        dt_remainder = subEdgeLength %% dt
#        
#        if (rootState == as.integer(names(edgeLength[j]))) {
#          for (k in 1:dt_length) {
#            xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt + sigmaRoot * sqrt(dt) * rnorm(1)
#            xt0 <- xt1
#          }
#          xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt_remainder + sigmaRoot * sqrt(dt_remainder) * rnorm(1)
#          xt0 <- xt1
#        }
#        else {
#          for (k in 1:dt_length) {
#            xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt + sigmaAlt * sqrt(dt) * rnorm(1)
#            xt0 <- xt1
#          }
#          xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt_remainder + sigmaAlt * sqrt(dt_remainder) * rnorm(1)
#          xt0 <- xt1
#        }
#        initial_state <- xt1
#        if 
#        next_nodes <- Descendants(currentNode[i], tree)
#        
#        obtainContinuousStates_ver5(tree = tree, alphaRoot = alphaRoot,
#                                    alphaAlt = alphaAlt, thetaRoot = thetaRoot, 
#                                    thetaAlt = thetaAlt, sigmaRoot = sigmaRoot, 
#                                    sigmaAlt = sigmaAlt,
#                                    initialState = initial_state, dt = 0.002,
#                                    currentNode = next_nodes)
#      }
#    }
#  }
#  return(contStates)
#}
#
#cont_states_ver5 <- obtainContinuousStates_ver5(tree = tree, alphaRoot = 1,
#                                                alphaAlt = 1, thetaRoot = 50,
#                                                thetaAlt = 20, sigmaRoot = 2,
#                                                sigmaAlt = 2, initialState = 50,
#                                                currentNode = NULL)

#obtainContinuousStates_ver6 = function(tree, alphaRoot, alphaAlt, thetaRoot,
#                                       thetaAlt, sigmaRoot, sigmaAlt,
#                                       initialState = thetaRoot, dt = 0.002,
#                                       currentNode = NULL, contStates = list()) {
#  cont_states <- contStates
#  ## obtain root state
#  rootState <- obtainRootState(tree)
#  
#  ##if it is at root node right now
#  if (is.null(currentNode[i]) == TRUE) {
#    currentNode <- Descendants(basal_node(tree), tree)
#  }
#  xt0 <- initialState
#  for (i in 1: length(currentNode)) {
#    edge <- which.edge(tree, currentNode[i])
#    edgeLength <- history$maps[[edge]]
#    for (j in 1:length(edgeLength)) {
#      subEdgeLength <- edgeLength[[j]]
#      dt_length = subEdgeLength %/% dt
#      dt_remainder = subEdgeLength %% dt
#      
#      if (rootState == as.integer(names(edgeLength[j]))) {
#        for (k in 1:dt_length) {
#          xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt + sigmaRoot * sqrt(dt) * rnorm(1)
#          xt0 <- xt1
#        }
#        xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt_remainder + sigmaRoot * sqrt(dt_remainder) * rnorm(1)
#        xt0 <- xt1
#      }
#      else {
#        for (k in 1:dt_length) {
#          xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt + sigmaAlt * sqrt(dt) * rnorm(1)
#          xt0 <- xt1
#        }
#        xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt_remainder + sigmaAlt * sqrt(dt_remainder) * rnorm(1)
#        xt0 <- xt1
#      }
#      initial_state <- xt1
#      
#      if (isTRUE(length(tips(tree, currentNode[i])) == 1)) {
#        cont_states[[tips(tree, currentNode[i])]] <- initial_state
#      }
#      else {
#        next_nodes <- Descendants(currentNode[i], tree)
#        return(obtainContinuousStates_ver6(tree = tree, alphaRoot = alphaRoot,
#                                           alphaAlt = alphaAlt, thetaRoot = thetaRoot,
#                                           thetaAlt = thetaAlt, sigmaRoot = sigmaRoot,
#                                           sigmaAlt = sigmaAlt,
#                                           initialState = initial_state, dt = 0.002,
#                                           currentNode = next_nodes,
#                                           contStates = cont_states))
#      }
#    }
#  }
#  return(cont_states)
#}
#
#cont_states_ver6 <- obtainContinuousStates_ver6(tree = tree, alphaRoot = 1,
#                                                alphaAlt = 1, thetaRoot = 50,
#                                                thetaAlt = 20, sigmaRoot = 2,
#                                                sigmaAlt = 2, initialState = 50)


obtainContinuousStates_ver7 = function(tree, history, alphaRoot, alphaAlt,
                                       thetaRoot, thetaAlt, sigmaRoot, sigmaAlt,
                                       initialValue = thetaRoot, dt = 0.002) {
  cont_states <- list()
  ## obtain root state
  root_state <- obtainRootState(tree)
  branch_order <- rev(postorder(tree))
  corr_nodes <- tree$edge
  node0 <- corr_nodes[branch_order[1], 1]
  node_values <- list()
  node_values[[as.character(node0)]] <- initialValue

  for (i in 1:length(branch_order)) {
    sub_edges <- history$maps[[branch_order[i]]]
    parent_node <- corr_nodes[branch_order[i], 1]
    xt0 <- node_values[[as.character(parent_node)]]
    for (j in 1:length(sub_edges)) {
      dt_length = sub_edges[j] %/% dt
      dt_remainder = sub_edges[j] %% dt
      
      if (root_state == as.integer(names(sub_edges[j]))) {
        for (k in 1:dt_length) {
          xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt + sigmaRoot * sqrt(dt) * rnorm(1)
          xt0 <- xt1
        }
        xt1 <- xt0 + alphaRoot * (thetaRoot - xt0) * dt_remainder + sigmaRoot * sqrt(dt_remainder) * rnorm(1)
        xt0 <- xt1
      }
      
      else {
        for (k in 1:dt_length) {
          xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt + sigmaAlt * sqrt(dt) * rnorm(1)
          xt0 <- xt1
        }
        xt1 <- xt0 + alphaAlt * (thetaAlt - xt0) * dt_remainder + sigmaAlt * sqrt(dt_remainder) * rnorm(1)
        xt0 <- xt1
      }
    }
    desc_node <- corr_nodes[branch_order[i], 2]
    node_values[[as.character(desc_node)]] <- xt0
  }
  for (i in 1:length(corr_nodes[,2])) {
    if (!(corr_nodes[i,2] %in% corr_nodes[,1])) {
      tip_label <- tips(tree, corr_nodes[i,2])
      cont_states[[tip_label]] <- unname(node_values[[as.character(corr_nodes[i,2])]])
    }
  }
  return(cont_states)
}
  
  
  
  
  
  
  
  
cont_states_ver7 <- obtainContinuousStates_ver7(tree = tree, history = history,
                                                alphaRoot = 0.5, alphaAlt = 1,
                                                thetaRoot = 50, thetaAlt = 80,
                                                sigmaRoot = 2, sigmaAlt = 2,
                                                initialValue = 60, dt = 0.002)



write.nexus.data(cont_states_ver7, file = "data/n50_simulationContinuous2.nex", format = "continuous")
var(unlist(cont_states_ver7))
