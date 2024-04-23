library(phytools)
library(pracma)
library(slouch)
library(TESS)
library(tidyverse)

parent <- function(phy, x){
  m <- which(phy$edge[, 2] == x)
  return(phy$edge[m, 1])
}

lineage.nodes <- function(phy, x){
  k <- x
  N <- length(phy$tip.label)
  while(x != N + 1){
    k <- c(k, parent(phy, x))
    x <- tail(k, n = 1)
  }
  return(k)
}

lineage.constructor <- function(phy, e, anc_maps, regimes, ace){
  nodes <- lineage.nodes(phy, e)
  node_ages <- ape::node.depth.edgelength(phy)[nodes]
  min_age <- min(node_ages)
  max_age <- max(node_ages)
  
  if(anc_maps == "simmap"){
    ## Simmap splits up each edge into sub-edges, depending on the split. So, we use edges instead of nodes, and introduce sub-edges
    edge_is <- which(phy$edge[,2] %in% nodes)
    subedges <- unlist(lapply(edge_is, function(i) phy$maps[[i]]))
    simmap_regimes <- rev(names(subedges))
    
    which.regimes <- lapply(levels(regimes), function(x) {res <- match(simmap_regimes, x); res[is.na(res)] <- 0; return(res)})
    # Problem. simmap does not provide root estimate. Assuming root estimate is equal to the oldest branch estimate
    root <- lapply(which.regimes, function(e) tail(e, n = 1))
    which.regimes <- lapply(seq_along(levels(regimes)), function(x) c(which.regimes[[x]], root[[x]]))
    
    times <- max_age - c(min_age, cumsum(unname(subedges))) ## time from the tip
    timeflip <- rev(times)
    # save the regimes in this lineage
    lineage_regimes <- names(subedges)
  }
  names(which.regimes) <- levels(regimes)
  
  ## time is t=0 at the tip, and increasing toward the past
  t_end <- head(timeflip, n = -1) ## Time from tip to the youngest point of segment(s)
  t_beginning <- tail(timeflip, n = -1) ## Time from tip to oldest point of segment(s)
  ## the "time spent in the root regime" is 0 myr
  regime_time <- c(t_beginning - t_end, 0.0) 
  
  return(list(nodes = nodes, 
              times = times,
              t_end = t_end,
              t_beginning = t_beginning,
              regime_time = regime_time,
              which.regimes = which.regimes,
              lineage_regimes = lineage_regimes))
}

weights_segments <- function(a, lineage){
  res <- c(exp(-a * lineage$t_end) - exp(-a * lineage$t_beginning), 
           exp(-a * lineage$times[1]))
  return(res)
}

weights_regimes <- function(a, lineage) {
  #nt <- lineage$nodes_time
  res <- weights_segments(a, lineage) ## Rcpp wrapper, equivalent to above commented out code
  w <- vapply(lineage$which.regimes, function(e) sum(unlist(e)*res), FUN.VALUE = 0) ## Sum coefficients for which the respective regime is equal
  return(w)
}

weight.matrix <- function(phy, a, lineages){
  res <- t(vapply(lineages, function(x) weights_regimes(a, x), 
                  FUN.VALUE = numeric(length(lineages[[1]]$which.regimes))) ## Specify type of output
  )
  
  rownames(res) <- phy$tip.label
  return(res)
}

theta_dep_logL_vcv <- function(tree, continuousChar, alpha, sigma2, theta){
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
  
  #X = matrix(1, ntip)
  
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
  
  l <- lapply(1:ntip, function(e) lineage.constructor(tree, e, "simmap", regimes = regimes))
  W = weight.matrix(tree, alpha, l)
  
  #r = solve(L) %*% y - solve(L) %*% X * theta # what does de-correlated residuals mean?
  r = solve(L) %*% y - solve(L) %*% W %*% theta
  # res = - (n/2) * log(2*pi) - 0.5 * log_det_V - 0.5 * dot(r, r)
  #     = exp(-n/2)^(2*pi) * exp(-0.5)^det_V * exp(-0.5)^dot(r, r) ?
  res = 0.0
  res = res - 0.5 * ntip * log(2*pi)
  res = res - 0.5 * log_det_V
  res = res - 0.5 * dot(r, r) # what is r and what is  dot product of r?
  
  return(res)
}

num_tips   = c(50, 100, 200, 500, 1000)
reps = 10
grid = expand.grid(num_tips=num_tips, tree=1:reps, stringsAsFactors=FALSE)

# 3 states
discChar <- list()
for (i in 1:length(num_tips)){
  ntip <- num_tips[i]
  discChar[[i]] <- c(rep(0,ntip*0.3), rep(1,ntip*0.6), rep(2,ntip*0.1))
  names(discChar[[i]]) <- paste0("t", as.character(c(1:50)))
}

contChar <- list()
for (i in 1:length(num_tips)){
  contChar[[i]] <- rnorm(num_tips[i], mean=0, sd=1)
  names(contChar[[i]]) <- paste0("t", as.character(c(1:50)))
}

alltrees <- list()
simmap_trees <- list()
# simulate the trees
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_tree       = this_row[[2]]

  alltrees[[i]] = ladderize(tess.sim.taxa(1, this_num_tips, 10, 1, 0.5)[[1]])
  
  # rescale the tree
  alltrees[[i]]$edge.length = alltrees[[i]]$edge.length / max(branching.times(alltrees[[i]]))
  
  j = i %% (length(num_tips)+1)
  simmap_trees[[i]] <- make.simmap(alltrees[[i]], discChar[[j]], model="SYM",nsim=1)
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
}




alpha = list()
sigma2 = list()
theta = list()
for (i in 1:50){
  alpha[[i]] = rgamma(n=3, shape=1, rate=10)
  #alpha[[i]] = rep(rgamma(n=1, shape=1, rate=10),3)
  names(alpha[[i]]) <- c("0", "1", "2")
  sigma2[[i]] = rgamma(n=3, shape=2, rate=10)
  #sigma2[[i]] = rep(rgamma(n=1, shape=3, rate=10),3)
  names(sigma2[[i]]) <- c("0", "1", "2")
  theta[[i]] = rnorm(3, mean = 0, sd = 3)
  names(theta[[i]]) <- c("0", "1", "2")
  }


time.taken <- tibble(ntips = num_tips, pruning_theta = 0, pruning_all = 0, vcv_theta=0)
for (i in length(num_tips)){
  j = i %% (length(num_tips)+1)
  #start.time <- Sys.time()
  #theta_dep_logL_vcv(simmap_trees[[i]], contChar[[j]], alpha[[i]][1], sigma2[[i]][1], theta[[i]])
  #theta_dep_logL_vcv(simmap_trees[[i+5]], contChar[[j]], alpha[[i+5]][1], sigma2[[i+5]][1], theta[[i+5]])
  #theta_dep_logL_vcv(simmap_trees[[i+10]], contChar[[j]], alpha[[i+10]][1], sigma2[[i+10]][1], theta[[i+10]])
  #theta_dep_logL_vcv(simmap_trees[[i+15]], contChar[[j]], alpha[[i+15]][1], sigma2[[i+15]][1], theta[[i+15]])
  #theta_dep_logL_vcv(simmap_trees[[i+20]], contChar[[j]], alpha[[i+20]][1], sigma2[[i+20]][1], theta[[i+20]])
  #theta_dep_logL_vcv(simmap_trees[[i+25]], contChar[[j]], alpha[[i+25]][1], sigma2[[i+25]][1], theta[[i+25]])
  #theta_dep_logL_vcv(simmap_trees[[i+30]], contChar[[j]], alpha[[i+30]][1], sigma2[[i+30]][1], theta[[i+30]])
  #theta_dep_logL_vcv(simmap_trees[[i+35]], contChar[[j]], alpha[[i+35]][1], sigma2[[i+35]][1], theta[[i+35]])
  #theta_dep_logL_vcv(simmap_trees[[i+40]], contChar[[j]], alpha[[i+40]][1], sigma2[[i+40]][1], theta[[i+40]])
  #theta_dep_logL_vcv(simmap_trees[[i+45]], contChar[[j]], alpha[[i+45]][1], sigma2[[i+45]][1], theta[[i+45]])
  #end.time <- Sys.time()
  #time.taken$vcv_theta[i] <- (end.time - start.time)/10
  
  start.time <- Sys.time()
  sd_logL_pruning(simmap_trees[[i]], contChar[[j]],
                  alpha=c("0"=alpha[[i]][1], "1"=alpha[[i]][1], "2"=alpha[[i]][1]),
                  sigma2=rep(sigma2[[i]][1], 3),
                  theta[[i]])
  sd_logL_pruning(simmap_trees[[i+5]], contChar[[j]], alpha[[i+5]][1], sigma2[[i+5]][1], theta[[i+5]])
  sd_logL_pruning(simmap_trees[[i+10]], contChar[[j]], alpha[[i+10]][1], sigma2[[i+10]][1], theta[[i+10]])
  sd_logL_pruning(simmap_trees[[i+15]], contChar[[j]], alpha[[i+15]][1], sigma2[[i+15]][1], theta[[i+15]])
  sd_logL_pruning(simmap_trees[[i+20]], contChar[[j]], alpha[[i+20]][1], sigma2[[i+20]][1], theta[[i+20]])
  sd_logL_pruning(simmap_trees[[i+25]], contChar[[j]], alpha[[i+25]][1], sigma2[[i+25]][1], theta[[i+25]])
  sd_logL_pruning(simmap_trees[[i+30]], contChar[[j]], alpha[[i+30]][1], sigma2[[i+30]][1], theta[[i+30]])
  sd_logL_pruning(simmap_trees[[i+35]], contChar[[j]], alpha[[i+35]][1], sigma2[[i+35]][1], theta[[i+35]])
  sd_logL_pruning(simmap_trees[[i+40]], contChar[[j]], alpha[[i+40]][1], sigma2[[i+40]][1], theta[[i+40]])
  sd_logL_pruning(simmap_trees[[i+45]], contChar[[j]], alpha[[i+45]][1], sigma2[[i+45]][1], theta[[i+45]])
  end.time <- Sys.time()
  time.taken$pruning_theta[i] <- (end.time - start.time)/10
  
  start.time <- Sys.time()
  sd_logL_pruning(simmap_trees[[i]], contChar[[j]], alpha[[i]], sigma2[[i]], theta[[i]])
  sd_logL_pruning(simmap_trees[[i+5]], contChar[[j]], alpha[[i+5]], sigma2[[i+5]], theta[[i+5]])
  sd_logL_pruning(simmap_trees[[i+10]], contChar[[j]], alpha[[i+10]], sigma2[[i+10]], theta[[i+10]])
  sd_logL_pruning(simmap_trees[[i+15]], contChar[[j]], alpha[[i+15]], sigma2[[i+15]], theta[[i+15]])
  sd_logL_pruning(simmap_trees[[i+20]], contChar[[j]], alpha[[i+20]], sigma2[[i+20]], theta[[i+20]])
  sd_logL_pruning(simmap_trees[[i+25]], contChar[[j]], alpha[[i+25]], sigma2[[i+25]], theta[[i+25]])
  sd_logL_pruning(simmap_trees[[i+30]], contChar[[j]], alpha[[i+30]], sigma2[[i+30]], theta[[i+30]])
  sd_logL_pruning(simmap_trees[[i+35]], contChar[[j]], alpha[[i+35]], sigma2[[i+35]], theta[[i+35]])
  sd_logL_pruning(simmap_trees[[i+40]], contChar[[j]], alpha[[i+40]], sigma2[[i+40]], theta[[i+40]])
  sd_logL_pruning(simmap_trees[[i+45]], contChar[[j]], alpha[[i+45]], sigma2[[i+45]], theta[[i+45]])
  end.time <- Sys.time()
  time.taken$pruning_all[i] <- (end.time - start.time)/10
}


# pruning all sd, pruning theta sd, vcv theta sd

save(time.taken, file = "data/1_validation/testing_speed/time.taken.Rda")







speed_plot <- ggplot(data=time.taken) +
  geom_line(aes(x=ntips, y=vcv_theta), color="#44AA99") +
  geom_line(aes(x=ntips, y=pruning_theta), color="#882255") +
  geom_point(aes(x=ntips, y=vcv_theta), color="#44AA99") +
  geom_point(aes(x=ntips, y=pruning_theta), color="#882255") +
  
  #geom_line(aes(x=ntips, y=pruning_all), color="#882255") +
  ggtitle("") +
  theme(panel.grid = element_blank(),
        #panel.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10))

pdf("speed_plot.pdf")
speed_plot
dev.off()
