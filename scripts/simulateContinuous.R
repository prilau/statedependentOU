library(ape)
library(phytools)
library(geiger)
library(TESS)
source("scripts/readWriteCharacterData.R")

# simulate the tree
num_tips = 250
tree = ladderize(tess.sim.taxa(1, num_tips, 10, 1, 0.5)[[1]])
# rescale the tree
tree$edge.length = tree$edge.length / max(branching.times(tree))

write.tree(tree, file=paste0("data/n100_simulation.tre"))


tree_lengths <- sum(tree$edge.length)
# specify rates so that the expected number of changes is 5
rates = 10 / tree_lengths
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

pdf("data/n100History.pdf")
plot(history, col=colors)
dev.off()

writeCharacterData(t(t(history$states)), file=paste0("data/n100_simulationDiscrete.nex"), type="Standard")






# obtain root state
obtainRootState = function(tree) {
  edge1d <- rev(postorder(tree))[1]
  rootState <- names(history$maps[[edge1d]][1])
  rootState <- as.integer(rootState)
  return(rootState)
}

obtainRootState(tree)


treeheight <- function(tree) max(node.depth.edgelength(tree))
obtainContinuousStates_ver7 = function(tree, halflifeRoot, halflifeAlt,
                                       thetaRoot, thetaAlt, sigmaRoot, sigmaAlt,
                                       initialValue = thetaRoot, dt = 0.002) {
  if (missing(dt)){
    dt <- 0.002 * treeheight(history)
  }
  ## Re-parameterization
  alphaRoot <- log(2) / halflifeRoot 
  alphaAlt <- log(2) / halflifeAlt
  
  cont_states <- list()
  ## obtain root state
  root_state <- obtainRootState(tree)
  preorder <- rev(postorder(tree))
  edges <- tree$edge
  ntips <- length(tree$tip.label)
  root_node <- ntips + 1
  node_values <- list()
  node_values[[root_node]] <- initialValue
  
  #for (i in 1:length(branch_order)) {
  for (edge_index in preorder){
    sub_edges <- tree$maps[[edge_index]]
    parent_node <- edges[edge_index, 1]
    xt0 <- node_values[[parent_node]]
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
    desc_node <- edges[edge_index, 2]
    node_values[[desc_node]] <- xt0
  }
  disc_states <- list()
  for (i in 1:length(edges[,2])) {
    is_tip <- !(edges[i,2] %in% edges[,1])
    if (is_tip) {
      tip_label <- tips(tree, edges[i,2])
      cont_states[[tip_label]] <- unname(node_values[[edges[i,2]]])
      disc_states[[tip_label]] <- tail(names(history$maps[[i]]), n = 1)
    }
  }
  res <- list(
    cont_states,
    disc_states
  )
  return(res)
}

plot(history)

cont_states_ver7 <- obtainContinuousStates_ver7(tree = history,
                                                halflifeRoot = 0.35, halflifeAlt = 0.35,
                                                thetaRoot = 50, thetaAlt = 20,
                                                sigmaRoot = 5, sigmaAlt = 5,
                                                initialValue = 50, dt = 0.001)

par(mar = c(6,7,3,3))
hist(as.numeric(cont_states_ver7[[1]]), xlab = "tip character values", ylab = "frequency")

treeheight(history)


library(slouch)
df <- data.frame(
  "species" = names(cont_states_ver7[[1]]),
  "y" = as.numeric(cont_states_ver7[[1]]),
  "regime" = unname(unlist(cont_states_ver7[[2]]))
)
df <- df[match(history$tip.label, df$species), ]

m0 <- slouch.fit(
  history,
  hl_values = seq(0.05, 0.25, length.out = 35),
  sigma2_y_values = seq(35, 65, length.out = 35),
  response = df$y,
  species = df$species,
  fixed.fact = df$regime,
  anc_maps = "simmap",
  hillclimb = F
)
m0
summary(m0)
plot(m0)
regimeplot(m0)
summary(m0)


write.nexus.data(cont_states_ver7[[1]]), file = "data/n100_simulationDiscrete.nex", format="Continuous")

var(unlist(cont_states_ver7[[1]]))
