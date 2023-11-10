library(ggridges)
library(ggplot2)
library(tidyverse)

cleanDf <- function(path, simID, treeSize) {
  messy_df <- read.csv(path, sep = "\t")
  clean_df <- messy_df %>% 
    filter(Iteration >= 240) %>% 
    select(halfLife, alpha.1., sigma2.1., stationaryVariance.1., theta.1., theta.2.) %>% 
    rename(alpha = alpha.1.,
           sigma2 = sigma2.1.,
           stationaryVariance = stationaryVariance.1.,
           theta_state0 = theta.1.,
           theta_state1 = theta.2.) %>% 
    mutate(simulationID = simID) %>%
    pivot_longer(cols = c("halfLife", "alpha", "sigma2", "stationaryVariance",
                          "theta_state0", "theta_state1"),
                 names_to = "parameter",
                 values_to = "value") %>% 
    mutate(tree_size = treeSize)
}



num_tips   = c(100, 250, 500)
tree_reps   = 1
num_dtraits = 4
num_ctraits = 4

grid = expand.grid(num_tips=num_tips, tree=1:tree_reps,
                   ctraits=1:num_ctraits, dtraits = 1:num_dtraits,
                   parameters = 1:6,
                   stringsAsFactors=FALSE)

this_dir = "output/"
looongPosteriorDf <- data.frame()


for(i in 1:nrow(grid)) {
  
  this_row = grid[i,]
  this_num_tips   = this_row[[1]]
  this_tree       = this_row[[2]]
  this_num_ctraits = this_row[[3]]
  this_num_dtraits = this_row[[4]]
  
  
  this_simID = paste0("n", this_num_tips,
                      "t", this_tree,
                      "d", this_num_dtraits,
                      "c", this_num_ctraits)
  this_path = paste0(this_dir, "sdOU_simulation2_", this_simID, ".log")
  clean_df <- cleanDf(this_path, this_simID, this_num_tips)
  looongPosteriorDf <- looongPosteriorDf %>% 
    bind_rows(clean_df)
}






# basic example
mean_parameters <- looongPosteriorDf %>% 
  group_by(parameter, simulationID) %>% 
  summarise(mean = mean(value), tree_size = mean(tree_size)) 

true_values <- tibble(
  parameter = c("halfLife", "stationaryVariance",
                "theta_state0", "theta_state1"),
  value = c(0.35, 0.0625, 0.5, 2.0))

#cuteLegoPlot <- mean_parameters %>% ggplot(aes(x = mean, alpha = 0.5, fill = factor(tree_size))) +
#  geom_histogram(bins = 10) +
#  facet_wrap(vars(parameter) , scales = "free_x") +
#  geom_vline(aes(xintercept = value, alpha = 0.5), data = true_values) +
#  theme_classic()


cuteTetrisPlot <- mean_parameters %>% 
  filter(parameter %in% c("halfLife", "stationaryVariance",
                          "theta_state0", "theta_state1")) %>%
  ggplot(aes(x = mean, alpha = 0.5, fill = factor(tree_size))) +
  geom_histogram(bins = 10) +
  facet_grid(cols = vars(parameter), rows = vars(tree_size), scales = "free_x") +
  geom_vline(aes(xintercept = value, alpha = 0.5), data = true_values) +
  theme_bw()
ggsave("figures/sim2_tetrisHist.pdf", cuteTetrisPlot, width = 200, height = 150, units = "mm")







mean_parameters_wide <- mean_parameters %>% pivot_wider(id_cols = ,
                                                        names_from = "tree_size",
                                                        values_from = "mean")
  


mean_n100 <- mean_parameters %>% filter(tree_size == 100)
mean_n250 <- mean_parameters %>% filter(tree_size == 250)
mean_n500 <- mean_parameters %>% filter(tree_size == 500)


mean_parameters %>% ggplot() +
  geom_histogram(data = mean_n100, aes(x = mean,
                 fill = "n100"), alpha = 0.2) +
  geom_histogram(data = mean_n250, aes(x = mean,
                 fill = "n250"), alpha = 0.2) +
  geom_histogram(data = mean_n100, aes(x = mean,
                fill = "n500"), alpha = 0.2) +
  geom_vline(aes(xintercept = value, alpha = 0.5), data = true_values) +
  scale_fill_manual(values = c("n100" = "red", "n250" = "green", "n500" = "blue")) +
  facet_wrap(vars(parameter) , scales = "free_x") +
  theme_classic()



ggsave("figures/sim2_histogram.pdf", cuteLegoPlot, width = 200, height = 150, units = "mm")







looongPosteriorDf %>% 
  filter(parameter == "halfLife") %>% 
  group_by(simulationID) %>% 
  ggplot(aes(x = value, y = simulationID, fill = simulationID, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(0, 2) +
  theme_ridges() + 
  theme(legend.position = "none")

looongPosteriorDf %>% 
  filter(parameter == "theta_state1") %>% 
  group_by(simulationID) %>% 
  ggplot(aes(x = value, y = simulationID, fill = simulationID, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(0, 5) +
  theme_ridges() + 
  theme(legend.position = "none")

looongPosteriorDf %>% 
  filter(parameter == "stationaryVariance") %>% 
  group_by(simulationID) %>% 
  ggplot(aes(x = value, y = simulationID, fill = simulationID, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(0, 0.25) +
  theme_ridges() + 
  theme(legend.position = "none")

looongPosteriorDf %>% 
  filter(parameter == "theta_state0") %>% 
  group_by(simulationID) %>% 
  ggplot(aes(x = value, y = simulationID, fill = simulationID, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(-0.1, 1.5) +
  theme_ridges() + 
  theme(legend.position = "none")





