library(ape)
library(phytools)
library(TESS)
library(pracma)
library(slouch)
library(tibble) 
library(tidyverse)
library(ggdist)
source("scripts/5_miscellaneous/functions.R")


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
grid_tmp <- grid %>%
  mutate(tree_size=log(num_tips),
         log_time_rb = log(time_rb)) %>%
  select(tree_size, log_time_rb)
coef(lm(grid_tmp$log_time_rb ~ grid_tmp$tree_size))
#grid_tmp <- grid %>% filter(!is.na(time_pruning)) %>% 
#  mutate(tree_size=log(num_tips),
#         ln_sec_pruning = log(time_pruning * 60)) %>%
#  select(tree_size, ln_sec_pruning)
#coef(lm(grid_tmp$ln_sec_pruning ~ grid_tmp$tree_size))



#grid %>% pivot_longer(cols=c(10,12), names_to = "method",
#                      values_to = "time") %>% 
#  ggplot(grid) +
#  geom_boxplot(aes(y=time, group=c(num_tips, method)))

#grid %>%
#  filter(!is.na(time_pruning)) %>% 
#  mutate(tree_size=log(num_tips),
#         sec_pruning = time_pruning * 60) %>% 
#  ggplot() +
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

p <- ggplot(grid) + 
  geom_point(aes(x=tree_size, y=log(time_rb), group=tree_size, color=as.factor(tree_size)),
           #position=position_jitter(width = 0.15),
           size=2, alpha=0.8) +
  geom_abline(linetype="dashed", alpha=0.5,
              intercept=-6.1151586,
              slope=0.8807052) + 
  scale_color_manual(values=c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99')) +
  scale_fill_manual(values=c('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE', '#882255', '#44AA99')) +
  scale_x_continuous(breaks=log(c(100, 200, 500, 1000, 5000, 10000, 100000)), label=c("100", "200", "500", "1k", "5k", "10k", "100k")) +
  scale_y_continuous(breaks=log(c(0.5, 1, 5, 10, 50)), label=c("0.5", "1", "5", "10", "50")) +
  ylab("Time taken per calculation (millisecond)") +
  theme_classic() +
  ggtitle("Likelihood algorithm in RevBayes") +
  theme(plot.title = element_text(hjust=0.5)) +
  xlab("Tree size") +
  theme(legend.position="none")

ggsave("figures/1_validation/rb_speed.pdf", p, height=3.4, width=3.4, unit="in")


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
