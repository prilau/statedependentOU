library(ape)
library(phytools)
library(geiger)
library(TESS)
#source("scripts/readWriteCharacterData.R")


treeheight <- function(tree) max(node.depth.edgelength(tree))
obtainContinuousStates_ver7 = function(tree, halflife_0, halflife_1,
                                       theta_0, theta_1, stationaryvar_0,
                                       stationaryvar_1, initialValue = theta_0,
                                       dt = 0.002) {
  if (missing(dt)){
    dt <- 0.002 * treeheight(history)
  }
  
  ## Re-parameterization
  alpha_0 <- log(2) / halflife_0 
  alpha_1 <- log(2) / halflife_1
  sigma_0 <- sqrt(stationaryvar_0 * 2 * alpha_0)
  sigma_1 <- sqrt(stationaryvar_1 * 2 * alpha_1)
  
  cont_states <- list()
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
      
      if (as.integer(names(sub_edges[j])) == 0) {
        for (k in 1:dt_length) {
          xt1 <- xt0 + alpha_0 * (theta_0 - xt0) * dt + sigma_0 * sqrt(dt) * rnorm(1)
          xt0 <- xt1
        }
        xt1 <- xt0 + alpha_0 * (theta_0 - xt0) * dt_remainder + sigma_0 * sqrt(dt_remainder) * rnorm(1)
        xt0 <- xt1
      }
      
      else {
        for (k in 1:dt_length) {
          xt1 <- xt0 + alpha_1 * (theta_1 - xt0) * dt + sigma_1 * sqrt(dt) * rnorm(1)
          xt0 <- xt1
        }
        xt1 <- xt0 + alpha_1 * (theta_1 - xt0) * dt_remainder + sigma_1 * sqrt(dt_remainder) * rnorm(1)
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




cat("simulating continuous characters.\n")

num_tips   = c(100, 250, 500)
reps       = 10
num_dtraits = 5
num_ctraits = 5

grid = expand.grid(num_tips=num_tips, tree=1:reps,
                   dtraits=1:num_dtraits, ctraits=1:num_ctraits,
                   stringsAsFactors=FALSE)
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_num_tips    = this_row[[1]]
  this_tree        = this_row[[2]]
  this_num_dtraits = this_row[[3]]
  this_num_ctraits = this_row[[4]]
  
  # read the history
  this_dir = paste0("data/n",this_num_tips, "/t", this_tree, "/d", this_num_dtraits)
  load(paste0(this_dir, "/n", this_num_tips,
              "t", this_tree, "d", this_num_dtraits, "_History.Rda"))
  
  cont_states_ver7 <- obtainContinuousStates_ver7(tree = history,
                                                  halflife_0 = 0.35, halflife_1 = 0.35,
                                                  theta_0 = 0.5, theta_1 = 2,
                                                  stationaryvar_0 = 0.0625, stationaryvar_1 = 0.0625,
                                                  initialValue = 0.5, dt = 0.001)
  
  #log_cont_states <- list()
  #for (i in 1:length(cont_states_ver7[[1]])) {
  #  tiplabel <- names(cont_states_ver7[[1]][i])
  #  log_value <- log(as.double(cont_states_ver7[[1]][i]))
  #  log_cont_states[[tiplabel]] <- log_value
  #}
  
  write.nexus.data(cont_states_ver7[[1]], 
                   file = paste0(this_dir, "/n", this_num_tips,
                                 "t", this_tree, "d", this_num_dtraits,
                                 "c", this_num_ctraits, "_Continuous.nex"),
                   format="Continuous")

  #write.nexus.data(log_cont_states, 
  #                 file = paste0(this_dir, "/n", this_num_tips,
  #                               "t", this_tree, "_logContinuous.nex"),
  #                 format="Continuous")
  setTxtProgressBar(bar, i / nrow(grid))
}
cat("\n")




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
  hl_values = seq(0.05, 0.5, length.out = 25),
  vy_values = seq(0.001, 0.2, length.out = 25),
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

