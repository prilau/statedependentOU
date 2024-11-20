library(RevGadgets)
library(tidyverse)
library(grid)
library(ggplot2)
library(kableExtra)
library(gridExtra)

fp <- tibble(sim = NA,
             theta = NA,
             alpha = NA,
             sigma2 = NA,
             halflife = NA,
             stv = NA)

num_sim = 1000
dir_in = "output/2_simulation/false_positive/"
dir_out = "figures/2_simulation/false_positive/"
pars <- c("theta", "alpha", "sigma2", "halflife", "stv")

bar = txtProgressBar(style=3, width=40)
for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_", i, "_run_1.log")
  df <- readTrace(file_in)[[1]] %>% 
    select(all_of(contains(c("theta_compare", "alpha_compare", "sigma2_compare", "halflife_compare", "stv_compare"))))
  fp[i,] = c(i, lapply(df, mean))
  setTxtProgressBar(bar, i / num_sim)
}
write.csv(fp, "output/2_simulation/false_positive/grid.csv")


thresholds <- tibble(sim=1:num_sim, upper=0.975, lower=0.025)
p <- list()
for (par in pars){
  fp_par <- fp %>% 
    pivot_longer(contains(par), names_to = "par")
  
  p[[par]] <- ggplot(fp_par) +
    geom_histogram(aes(x=value), binwidth = 0.05) +
    theme_classic() +
    geom_vline(mapping = aes(xintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
    geom_vline(mapping = aes(xintercept = upper), data=thresholds, linetype = "dashed", color = "brown") +
    ggtitle(par) +
    #scale_y_continuous(breaks = c(0, 20, 40)) +
    theme(axis.title = element_blank())
}

# plot table
fp_bin <- fp %>% 
  pivot_longer(contains("_"), names_to = "par") %>% 
  mutate(value = ifelse(value > 0.975 | value < 0.025, 1, 0)) %>% 
  group_by(par) %>%
  summarise(`False positive rate` = signif(mean(value),2)) %>% 
  rename(Parameter = par)

fpr <- tableGrob(fp_bin, theme=ttheme_minimal(), rows=NULL)

plot_layout <- 
  "ABC
   DEF"  

p_all <- wrap_plots(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], fpr, design = plot_layout)

file_out <- paste0(dir_out, "fpr.pdf")
ggsave(file_out, p_all, width = 300, height = 180, units = "mm")


