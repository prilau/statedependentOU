library(phytools)
library(slouch)
library(TESS)
library(tidyverse)

num_tips   = c(50, 100, 200, 500, 1000)
reps = 10
grid = expand.grid(num_tips=num_tips, tree=1:reps, stringsAsFactors=FALSE)
this_dir = paste0("data/1_validation/testing_speed/")
alltrees <- list()

# simulate the trees
bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_tree       = this_row[[2]]

  alltrees[[i]] = ladderize(tess.sim.taxa(1, this_num_tips, 10, 1, 0.5)[[1]])
  
  # rescale the tree
  alltrees[[i]]$edge.length = tree$edge.length / max(branching.times(tree))
  
  # increment the progress bar
  setTxtProgressBar(bar, i / nrow(grid))
}

# 3 states
disc_50 <- c(rep(0,15), rep(1,30), rep(2,5))
names(disc_50) <- paste0("t", as.character(c(1:50)))
disc_100 <- c(rep(0,30), rep(1,60), rep(2,10))
names(disc_100) <- paste0("t", as.character(c(1:100)))
disc_200 <- c(rep(0,60), rep(1,125), rep(2,15))
names(disc_200) <- paste0("t", as.character(c(1:200)))
disc_500 <- c(rep(0,150), rep(1,325), rep(2,25))
names(disc_500) <- paste0("t", as.character(c(1:500)))
disc_1000 <- c(rep(0,300), rep(1,650), rep(2,50))
names(disc_1000) <- paste0("t", as.character(c(1:1000)))

simmap_trees <- list()
trees_names <- tibble(it = 1:50, treename = "0")
for(i in 1:nrow(grid)) {
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_tree       = this_row[[2]]
  
  if (this_num_tips==50) {
    simmap_trees[[i]] <- make.simmap(alltrees[[i]], disc_50, model="SYM",nsim=1)
  }
  else if (this_num_tips==100) {
    simmap_trees[[i]] <- make.simmap(alltrees[[i]], disc_100, model="SYM",nsim=1)
  }
  else if (this_num_tips==200) {
    simmap_trees[[i]] <- make.simmap(alltrees[[i]], disc_200, model="SYM",nsim=1)
  }
  else if (this_num_tips==500) {
    simmap_trees[[i]] <- make.simmap(alltrees[[i]], disc_500, model="SYM",nsim=1)
  }
  else if (this_num_tips==1000) {
    simmap_trees[[i]] <- make.simmap(alltrees[[i]], disc_1000, model="SYM",nsim=1)
  }
  trees_names$treename <- paste0(this_num_tips, "_", this_tree)
}


cont_50 <- rnorm(50, mean=0, sd=1)
names(cont_50) <- paste0("t", as.character(c(1:50)))
cont_100 <- rnorm(100, mean=0, sd=1)
names(cont_100) <- paste0("t", as.character(c(1:100)))
cont_200 <- rnorm(200, mean=0, sd=1)
names(cont_200) <- paste0("t", as.character(c(1:200)))
cont_500 <- rnorm(500, mean=0, sd=1)
names(cont_500) <- paste0("t", as.character(c(1:500)))
cont_1000 <- rnorm(1000, mean=0, sd=1)
names(cont_1000) <- paste0("t", as.character(c(1:1000)))

alpha = list()
sigma2 = list()
theta = list()
for (i in 1:50){
  alpha[[i]] = rgamma(n=3, shape=1, rate=10)
  names(alpha[[i]]) <- c("0", "1", "2")
  sigma2[[i]] = rgamma(n=3, shape=2, rate=10)
  names(sigma2[[i]]) <- c("0", "1", "2")
  theta[[i]] = rnorm(3, mean = 0, sd = 3)
  names(theta[[i]]) <- c("0", "1", "2")
  }


time.taken.all <- tibble(ntips = c(50,100,200,500,1000), mean_time = 0)

start.time <- Sys.time()
sd_logL_pruning(simmap_trees[[1]], cont_50, alpha[[1]], sigma2[[1]], theta[[1]])
sd_logL_pruning(simmap_trees[[6]], cont_50, alpha[[6]], sigma2[[6]], theta[[6]])
sd_logL_pruning(simmap_trees[[11]], cont_50, alpha[[11]], sigma2[[11]], theta[[11]])
sd_logL_pruning(simmap_trees[[21]], cont_50, alpha[[21]], sigma2[[21]], theta[[21]])
sd_logL_pruning(simmap_trees[[31]], cont_50, alpha[[31]], sigma2[[31]], theta[[31]])
sd_logL_pruning(simmap_trees[[41]], cont_50, alpha[[41]], sigma2[[41]], theta[[41]])
sd_logL_pruning(simmap_trees[[16]], cont_50, alpha[[16]], sigma2[[16]], theta[[16]])
sd_logL_pruning(simmap_trees[[26]], cont_50, alpha[[26]], sigma2[[26]], theta[[26]])
sd_logL_pruning(simmap_trees[[36]], cont_50, alpha[[36]], sigma2[[36]], theta[[36]])
sd_logL_pruning(simmap_trees[[46]], cont_50, alpha[[46]], sigma2[[46]], theta[[46]])
end.time <- Sys.time()
time.taken.all$mean_time[1] <- (end.time - start.time)/10

start.time <- Sys.time()
sd_logL_pruning(simmap_trees[[2]], cont_100, alpha[[2]], sigma2[[2]], theta[[2]])
sd_logL_pruning(simmap_trees[[7]], cont_100, alpha[[7]], sigma2[[7]], theta[[7]])
sd_logL_pruning(simmap_trees[[12]], cont_100, alpha[[12]], sigma2[[12]], theta[[12]])
sd_logL_pruning(simmap_trees[[22]], cont_100, alpha[[22]], sigma2[[22]], theta[[22]])
sd_logL_pruning(simmap_trees[[32]], cont_100, alpha[[32]], sigma2[[32]], theta[[32]])
sd_logL_pruning(simmap_trees[[42]], cont_100, alpha[[42]], sigma2[[42]], theta[[42]])
sd_logL_pruning(simmap_trees[[17]], cont_100, alpha[[17]], sigma2[[17]], theta[[17]])
sd_logL_pruning(simmap_trees[[27]], cont_100, alpha[[27]], sigma2[[27]], theta[[27]])
sd_logL_pruning(simmap_trees[[37]], cont_100, alpha[[37]], sigma2[[37]], theta[[37]])
sd_logL_pruning(simmap_trees[[47]], cont_100, alpha[[47]], sigma2[[47]], theta[[47]])
end.time <- Sys.time()
time.taken.all$mean_time[2] <- (end.time - start.time)/10

start.time <- Sys.time()
sd_logL_pruning(simmap_trees[[3]], cont_200, alpha[[3]], sigma2[[3]], theta[[3]])
sd_logL_pruning(simmap_trees[[8]], cont_200, alpha[[8]], sigma2[[8]], theta[[8]])
sd_logL_pruning(simmap_trees[[13]], cont_200, alpha[[13]], sigma2[[13]], theta[[13]])
sd_logL_pruning(simmap_trees[[23]], cont_200, alpha[[23]], sigma2[[23]], theta[[23]])
sd_logL_pruning(simmap_trees[[33]], cont_200, alpha[[33]], sigma2[[33]], theta[[33]])
sd_logL_pruning(simmap_trees[[43]], cont_200, alpha[[43]], sigma2[[43]], theta[[43]])
sd_logL_pruning(simmap_trees[[18]], cont_200, alpha[[18]], sigma2[[18]], theta[[18]])
sd_logL_pruning(simmap_trees[[28]], cont_200, alpha[[28]], sigma2[[28]], theta[[28]])
sd_logL_pruning(simmap_trees[[38]], cont_200, alpha[[38]], sigma2[[38]], theta[[38]])
sd_logL_pruning(simmap_trees[[48]], cont_200, alpha[[48]], sigma2[[48]], theta[[48]])
end.time <- Sys.time()
time.taken.all$mean_time[3] <- (end.time - start.time)/10

start.time <- Sys.time()
sd_logL_pruning(simmap_trees[[4]], cont_500, alpha[[4]], sigma2[[4]], theta[[4]])
sd_logL_pruning(simmap_trees[[9]], cont_500, alpha[[9]], sigma2[[9]], theta[[9]])
sd_logL_pruning(simmap_trees[[14]], cont_500, alpha[[14]], sigma2[[14]], theta[[14]])
sd_logL_pruning(simmap_trees[[24]], cont_500, alpha[[24]], sigma2[[24]], theta[[24]])
sd_logL_pruning(simmap_trees[[34]], cont_500, alpha[[34]], sigma2[[34]], theta[[34]])
sd_logL_pruning(simmap_trees[[44]], cont_500, alpha[[44]], sigma2[[44]], theta[[44]])
sd_logL_pruning(simmap_trees[[19]], cont_500, alpha[[19]], sigma2[[19]], theta[[19]])
sd_logL_pruning(simmap_trees[[29]], cont_500, alpha[[29]], sigma2[[29]], theta[[29]])
sd_logL_pruning(simmap_trees[[39]], cont_500, alpha[[39]], sigma2[[39]], theta[[39]])
sd_logL_pruning(simmap_trees[[49]], cont_500, alpha[[49]], sigma2[[49]], theta[[49]])
end.time <- Sys.time()
time.taken.all$mean_time[4] <- (end.time - start.time)/10

start.time <- Sys.time()
sd_logL_pruning(simmap_trees[[5]], cont_1000, alpha[[5]], sigma2[[5]], theta[[5]])
sd_logL_pruning(simmap_trees[[10]], cont_1000, alpha[[10]], sigma2[[10]], theta[[10]])
sd_logL_pruning(simmap_trees[[15]], cont_1000, alpha[[15]], sigma2[[15]], theta[[15]])
sd_logL_pruning(simmap_trees[[25]], cont_1000, alpha[[25]], sigma2[[25]], theta[[25]])
sd_logL_pruning(simmap_trees[[35]], cont_1000, alpha[[35]], sigma2[[35]], theta[[35]])
sd_logL_pruning(simmap_trees[[45]], cont_1000, alpha[[45]], sigma2[[45]], theta[[45]])
sd_logL_pruning(simmap_trees[[50]], cont_1000, alpha[[50]], sigma2[[50]], theta[[50]])
sd_logL_pruning(simmap_trees[[20]], cont_1000, alpha[[20]], sigma2[[20]], theta[[20]])
sd_logL_pruning(simmap_trees[[30]], cont_1000, alpha[[30]], sigma2[[30]], theta[[30]])
sd_logL_pruning(simmap_trees[[40]], cont_1000, alpha[[40]], sigma2[[40]], theta[[40]])
end.time <- Sys.time()
time.taken.all$mean_time[5] <- (end.time - start.time)/10

# pruning all sd, pruning theta sd, vcv theta sd
time.taken <- tibble(ntips = c(50,100,200,500,1000), pruning_theta = 0, pruning_all = 0)
time.taken$pruning_theta <- time.taken.theta$mean_time
time.taken$pruning_all <- time.taken.all$mean_time
save(time.taken, file = "data/1_validation/testing_speed/time.taken.Rda")







ggplot(data=time.taken) +
  geom_point(aes(x=ntips, y=pruning_theta), color="#0077bb") +
  geom_point(aes(x=ntips, y=pruning_all), color="#ee7733") +
  geom_point(aes(x=ntips, y=slouch_theta), color="#050aaa")
