library(ggplot2)
library(ggridges)
library(tidyverse)
library(patchwork)

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
  xlim(0, 3) +
  theme_ridges() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("half life")

cuteRidgePlot_sigma2 <- looongPosteriorDf %>% 
  filter(parameter %in% c("sigma2.1.", "sigma2.2.")) %>% 
  group_by(trait) %>% 
  ggplot(aes(x = value, y = trait, fill = parameter, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(0, 3) +
  theme_ridges() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("sigma2")

cuteRidgePlot_theta <- looongPosteriorDf %>% 
  filter(parameter %in% c("theta_state0", "theta_state1")) %>% 
  group_by(trait) %>% 
  ggplot(aes(x = value, y = trait, fill = parameter, alpha = 0.5)) +
  geom_density_ridges() +
  xlim(-0.4, 80) +
  theme_ridges() + 
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("theta")


nested <- (cuteRidgePlot_sigma2|cuteRidgePlot_halfLife|cuteRidgePlot_theta)+
  plot_annotation(tag_levels = 'A') #add figure labels

ggsave("figures/haemulidae_all_parameters.pdf", nested, width = 400, height = 200, units = "mm")
