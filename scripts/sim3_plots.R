library(ggridges)
library(ggplot2)
library(tidyverse)

cleanDf <- function(path, simID, expChange) {
  messy_df <- read.csv(path, sep = "\t")
  clean_df <- messy_df %>% 
    filter(Iteration >= 240) %>% 
    select(halfLife, sigma2Stateless, theta.1., theta.2.) %>% 
    rename(theta_state0 = theta.1.,
           theta_state1 = theta.2.) %>% 
    mutate(simulationID = simID) %>%
    pivot_longer(cols = c("halfLife", "sigma2Stateless",
                          "theta_state0", "theta_state1"),
                 names_to = "parameter",
                 values_to = "value") %>% 
    mutate(exp_change = expChange)
}



num_tips   = 500
tree_reps   = 1
rate = c(5, 10, 20, 50, 500)
num_dtraits = 4
num_ctraits = 4

grid = expand.grid(num_tips=num_tips, tree=1:tree_reps,
                   ctraits=1:num_ctraits, dtraits = 1:num_dtraits,
                   rates = rate,
                   stringsAsFactors=FALSE)

this_dir = "output/"
looongPosteriorDf <- data.frame()


for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_tree       = this_row[[2]]
  this_num_ctraits = this_row[[3]]
  this_num_dtraits = this_row[[4]]
  this_rate = this_row[[5]]
  
  
  this_simID = paste0("n", this_num_tips,
                      "t", this_tree,
                      "r", this_rate,
                      "d", this_num_dtraits,
                      "c", this_num_ctraits)
  this_path = paste0(this_dir, "sdOU_simulation3_", this_simID, ".log")
  clean_df <- cleanDf(this_path, this_simID, this_rate)
  looongPosteriorDf <- looongPosteriorDf %>% 
    bind_rows(clean_df)
}






# basic example
mean_parameters <- looongPosteriorDf %>% 
  group_by(parameter, simulationID) %>% 
  summarise(mean = mean(value), exp_change = mean(exp_change)) 

true_values <- tibble(
  parameter = c("halfLife", "sigma2Stateless",
                "theta_state0", "theta_state1"),
  value = c(0.35, 0.0625*2*log(2)/0.35, 0.5, 2.0))


cuteTetrisPlot <- mean_parameters %>% 
  filter(parameter %in% c("halfLife", "sigma2Stateless",
                          "theta_state0", "theta_state1")) %>%
  ggplot(aes(x = mean, alpha = 0.5, fill = factor(exp_change))) +
  geom_histogram(bins = 15) +
  facet_grid(cols = vars(parameter), rows = vars(exp_change), scales = "free_x") +
  geom_vline(aes(xintercept = value, alpha = 0.5), data = true_values) +
  theme_bw()
ggsave("figures/sim3_tetrisHist.pdf", cuteTetrisPlot, width = 200, height = 150, units = "mm")
