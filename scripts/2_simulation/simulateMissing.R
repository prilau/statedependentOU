library(ape)
library(phytools)
library(tidyverse)
source("scripts/5_miscellaneous/readWriteCharacterData.R")

dir_in = "data/2_simulation/missing_state/"
dir_out = "data/2_simulation/missing_state/"

num_sim=10

tree <- read.tree("data/1_validation/artiodactyla/artiodactyla.tree")
tree_length = sum(tree$edge.length)
rate = 20 / tree_length # expected number of changes is 50

Q2 = matrix(1, 2, 2, dimnames=list(1:2, 1:2))
diag(Q2) = -1

Q3 = matrix(c(-0.75, 0.75, 0,
              0.75, -1.50, 0.75,
              0, 0.75, -0.75),
            nrow=3,
            dimnames=list(1:3, 1:3))



Q4 = matrix(c(-2/3, 2/3, 0, 0,
              2/3, -4/3, 2/3, 0,
              0, 2/3, -4/3, 2/3,
              0, 0, 2/3, -2/3),
            nrow=4,
            dimnames=list(1:4, 1:4))

bar = txtProgressBar(style=3, width=40)
for(i in 1:num_sim) {

  history2 = sim.history(tree, rate * Q2, nsim=1, message=FALSE)
  while (! (mean(history2$states == "1") > 0.2 & mean(history2$states == "2") > 0.2 ) ) {
    history2 = sim.history(tree, rate * Q2, nsim=1, message=FALSE)
  }
  maps = history2$mapped.edge[,c("1","2")]
  
  history3 = sim.history(tree, rate * Q3, nsim=1, message=FALSE)
  while (! (mean(history3$states == "1") > 0.2 & mean(history3$states == "2") > 0.2 & mean(history3$states == "3") > 0.2) ) {
    history3 = sim.history(tree, rate * Q3, nsim=1, message=FALSE)
  }
  maps = history3$mapped.edge[,c("1","2","3")]
  
  history4 = sim.history(tree, rate * Q4, nsim=1, message=FALSE)
  while (! (mean(history4$states == "1") > 0.2 & mean(history4$states == "2") > 0.2 & mean(history4$states == "3") > 0.2 & mean(history4$states == "4") > 0.2) ) {
    history4 = sim.history(tree, rate * Q4, nsim=1, message=FALSE)
  }
  maps = history4$mapped.edge[,c("1","2","3","4")]
  
  
  # save these trees
  this_sub_dir = paste0("data/2_simulation/missing_state/sim_", i)
  if ( !dir.exists(this_sub_dir) ) {
    dir.create(this_sub_dir, recursive=TRUE, showWarnings=FALSE)
  }
  
  # save the discrete trait as a nexus file
  writeCharacterData(t(t(history2$states)), file=paste0(this_sub_dir, "/nstate_2_discrete.nex"), type="Standard")
  writeCharacterData(t(t(history3$states)), file=paste0(this_sub_dir, "/nstate_3_discrete.nex"), type="Standard")
  writeCharacterData(t(t(history4$states)), file=paste0(this_sub_dir, "/nstate_4_discrete.nex"), type="Standard")
  
  # save the character history
  save(history2, file=paste0(this_sub_dir, "/nstate_2_history.Rda"))
  save(history3, file=paste0(this_sub_dir, "/nstate_3_history.Rda"))
  save(history4, file=paste0(this_sub_dir, "/nstate_4_history.Rda"))
  
  
  # increment the progress bar
  setTxtProgressBar(bar, i / num_sim)
  
} 
cat("\n")

# Root age
root_age <- max(node.depth.edgelength(tree))

# 2 missing states
halflives <- c("1"=0.3*root_age, "2"=0.6*root_age)
this_alpha <- log(2) / halflives
this_stv <- c("1"=1, "2"=2)
this_sigma2 <- 2 * this_alpha * this_stv
this_theta <- c("1"=-6, "2"=6)

v2 <- c()

bar = txtProgressBar(style=3, width=40)
for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_",
                    i, "/nstate_2_history.Rda")
  load(file_in)

  sim <- simulateContinuous(tree=history2, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  v2[i] <- var(sim %>% unlist())
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/nstate_2_continuous.nex"),
                   format="Continuous")
  setTxtProgressBar(bar, i / num_sim)
}
write_tsv(as_tibble(v2), file=paste0(dir_in, "nstate_2_var_cont.txt"), col_names = FALSE)

# 3 missing states
halflives <- c("1"=0.3*root_age, "2"=0.6*root_age, "3"=0.6*root_age)
this_alpha <- log(2) / halflives
this_stv <- c("1"=1, "2"=2, "3"=2)
this_sigma2 <- 2 * this_alpha * this_stv
this_theta <- c("1"=-6, "2"=0, "3"=6)

v3 <- c()

bar = txtProgressBar(style=3, width=40)
for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_",
                    i, "/nstate_3_history.Rda")
  load(file_in)

  sim <- simulateContinuous(tree=history3, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  v3[i] <- var(sim %>% unlist())
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/nstate_3_continuous.nex"),
                   format="Continuous")
  setTxtProgressBar(bar, i / num_sim)
}
write_tsv(as_tibble(v3), file=paste0(dir_in, "nstate_3_var_cont.txt"), col_names = FALSE)

# 4 missing states
halflives <- c("1"=0.3*root_age, "2"=0.6*root_age, "3"=0.6*root_age, "4"=0.6*root_age)
this_alpha <- log(2) / halflives
this_stv <- c("1"=1, "2"=2, "3"=2, "4"=2)
this_sigma2 <- 2 * this_alpha * this_stv
this_theta <- c("1"=-6, "2"=-2, "3"=2, "4"=6)

v4 <- c()

bar = txtProgressBar(style=3, width=40)
for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_",
                    i, "/nstate_4_history.Rda")
  load(file_in)
  
  sim <- simulateContinuous(tree=history4, alpha=this_alpha, sigma2=this_sigma2,
                            theta=this_theta)
  v4[i] <- var(sim %>% unlist())
  this_dir <- paste0(dir_out, "sim_", i)
  write.nexus.data(sim, file = paste0(this_dir, "/nstate_4_continuous.nex"),
                   format="Continuous")
  setTxtProgressBar(bar, i / num_sim)
}
write_tsv(as_tibble(v4), file=paste0(dir_in, "nstate_4_var_cont.txt"), col_names = FALSE)