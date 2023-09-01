library(ape)
library(phytools)
library(geiger)
library(TESS)
#source("scripts/readWriteCharacterData.R")


# obtain root state
obtainRootState = function(tree) {
  edge1d <- rev(postorder(tree))[1]
  rootState <- names(history$maps[[edge1d]][1])
  rootState <- as.integer(rootState)
  return(rootState)
}


treeheight <- function(tree) max(node.depth.edgelength(tree))
obtainContinuousStates_ver7 = function(tree, halflifeRoot, halflifeAlt,
                                       thetaRoot, thetaAlt, stationaryvarRoot,
                                       stationaryvarAlt, initialValue = thetaRoot,
                                       dt = 0.002) {
  if (missing(dt)){
    dt <- 0.002 * treeheight(history)
  }
  
  ## Re-parameterization
  alphaRoot <- log(2) / halflifeRoot 
  alphaAlt <- log(2) / halflifeAlt
  sigmaRoot <- sqrt(stationaryvarRoot * 2 * alphaRoot)
  sigmaAlt <- sqrt(stationaryvarAlt * 2 * alphaAlt)
  
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
                                                  halflifeRoot = 0.35, halflifeAlt = 0.35,
                                                  thetaRoot = 0.5, thetaAlt = 2,
                                                  stationaryvarRoot = 0.0625, stationaryvarAlt = 0.0625,
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


data <- c(
  0.961836956464556,
  1.10837213261247,
  1.0926195360605,
  0.795822732390246,
  0.985575448959426,
  1.056369274489,
  0.987839911351482,
  1.02218819324461,
  1.10071874767439,
  0.603395233048221,
  0.661391519895569,
  0.37657223519301,
  0.662461273470934,
  0.668882152361772,
  0.669237004416171,
  0.472232991425192,
  0.410034503028078,
  0.434754645768031,
  0.420006944134845,
  0.410643193456357,
  0.841620923492179,
  0.750371042573203,
  0.798478984490198,
  0.727478313842834,
  0.873711002012752,
  0.670507891150641,
  0.761558968836325,
  0.85597680360232,
  1.04395211752431,
  0.278215670647991,
  0.275485555652166,
  0.236342206436231,
  0.217685246267636,
  0.338632366667114,
  0.119557982197699,
  0.192439486642763,
  0.261360873939005,
  0.191497774851475,
  0.509536352587912,
  0.431189700959807,
  0.20460492027561,
  0.442863085605562,
  0.755962033137867,
  0.792988985136037,
  0.969995910032833,
  0.531843453340787,
  0.679893311617496,
  0.982600975582826,
  0.913240895438416,
  0.915677146910292,
  0.920370483910389,
  0.678761758755848,
  1.03970615229784,
  1.10597754740447,
  0.98495050131774,
  1.14621719939735,
  0.780011663357825,
  0.81595114952783,
  0.390736774320853,
  0.134670575743913,
  0.131740307603808,
  0.106820661277778,
  0.167635437537387,
  0.594927500762657,
  0.26929719435948,
  0.331744483382552,
  0.803795975541894,
  0.620092054518424,
  0.653166440107539,
  0.740906917658979,
  0.95685708873794,
  0.270384118893878,
  0.271884286149036,
  0.453650166349657,
  0.587292024760133,
  0.601014417343575,
  0.19733596103397,
  0.0871656595398915,
  0.111083621381402,
  1.06234692978812,
  1.12856621987065,
  1.09889154247018,
  1.16475235900096,
  0.868691873674262,
  1.18207992265349,
  1.14223224068461,
  1.20518946483718,
  1.45129709684034,
  1.31951675810736,
  1.25048563111712,
  1.10463954577296,
  1.27534865419882,
  1.23131496844669,
  1.28853437231597,
  1.08398510081215,
  1.03030328207124,
  0.861580913098945,
  0.338733670467138,
  0.323339949174712,
  0.324980035355324
)
var(data)
