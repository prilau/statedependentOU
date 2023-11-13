library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(ggridges)
library(tidyverse)

cleanDf <- function(path, simID) {
  messy_df <- read.csv(path, sep = "\t")
  clean_df <- messy_df %>% 
    filter(Iteration >= 240) %>% 
    select(halfLife.1., halfLife.2., sigma2.1., sigma2.2., theta.1., theta.2.) %>% 
    rename(theta_state0 = theta.1.,
           theta_state1 = theta.2.) %>% 
    mutate(trait = cont_trait_labels$trait[which(cont_trait_labels$index == simID)]) %>%
    pivot_longer(cols = c("halfLife.1.", "halfLife.2.", "sigma2.1.", "sigma2.2.",
                          "theta_state0", "theta_state1"),
                 names_to = "parameter",
                 values_to = "value")
}

cont_trait_labels <- tibble(
  index = c(1, 2, 5, 6, 7, 9, 10, 11, 12),
  trait = c("body size",
            "adductor mass",
            "ascending process",
            "raker length",
            "eye width",
            "buccal length",
            "buccal width",
            "head height",
            "head length"))


looongPosteriorDf <- data.frame()


for(i in c(1,2,5,6,7,9,10,11,12)) {
  this_path = paste0("output/4_empirical_haemulidae_", i, ".log")
  clean_df <- cleanDf(this_path, i)
  looongPosteriorDf <- looongPosteriorDf %>% 
    bind_rows(clean_df)
}






cuteRidgePlot_halfLife <- looongPosteriorDf %>% 
  filter(parameter %in% c("halfLife.1.", "halfLife.2.")) %>% 
  group_by(trait) %>% 
  ggplot(aes(x = value, y = trait, fill = parameter, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(-0.5, 4) +
  theme_ridges() + 
  theme(legend.position = "none") 
ggsave("figures/haemulidae_halfLife.pdf", cuteRidgePlot_halfLife, width = 200, height = 150, units = "mm")

cuteRidgePlot_sigma2 <- looongPosteriorDf %>% 
  filter(parameter %in% c("sigma2.1.", "sigma2.2.")) %>% 
  group_by(trait) %>% 
  ggplot(aes(x = value, y = trait, fill = parameter, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(-0.1, 4) +
  theme_ridges() + 
  theme(legend.position = "none") 
ggsave("figures/haemulidae_sigma2.pdf", cuteRidgePlot_sigma2, width = 200, height = 150, units = "mm")

cuteRidgePlot_theta <- looongPosteriorDf %>% 
  filter(parameter %in% c("theta_state0", "theta_state1")) %>% 
  group_by(trait) %>% 
  ggplot(aes(x = value, y = trait, fill = parameter, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(-0.4, 80) +
  theme_ridges() + 
  theme(legend.position = "none") 
ggsave("figures/haemulidae_theta.pdf", cuteRidgePlot_theta, width = 200, height = 150, units = "mm")
