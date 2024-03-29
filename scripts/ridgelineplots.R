library(ggridges)
library(ggplot2)
library(tidyverse)

cleanDf <- function(path, simID) {
  messy_df <- read.csv(path, sep = "\t")
  clean_df <- messy_df %>% 
    filter(Iteration >= 240) %>% 
    select(halfLife, stationaryVariance.1., theta.1., theta.2.) %>% 
    rename(stationaryVariance = stationaryVariance.1.,
           theta_state0 = theta.1.,
           theta_state1 = theta.2.) %>% 
    mutate(simulationID = simID) %>%
    pivot_longer(cols = c("halfLife", "stationaryVariance", "theta_state0", "theta_state1"),
                 names_to = "parameter",
                 values_to = "value")
}



num_tips   = c(100, 250, 500)
tree_reps   = 1
num_dtraits = 4
num_ctraits = 4

grid = expand.grid(num_tips=num_tips, tree=1:tree_reps,
                   ctraits=1:num_ctraits, dtraits = 1:num_dtraits,
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
  clean_df <- cleanDf(this_path, this_simID)
  looongPosteriorDf <- looongPosteriorDf %>% 
    bind_rows(clean_df)
}





# basic example

mean_parameters <- looongPosteriorDf %>% 
  group_by(parameter, simulationID) %>% 
  summarise(mean = mean(value))
  
true_values <- tibble(
  parameter = c("halfLife", "stationaryVariance", "theta_state0", "theta_state1"),
  value = c(0.35, 0.0625, 0.5, 2.0))

mean_parameters %>% ggplot(aes(x = mean, alpha = 0.5)) +
  geom_histogram() +
  facet_wrap(vars(parameter) , scales = "free_x") +
  geom_vline(aes(xintercept = value), data = true_values)


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






