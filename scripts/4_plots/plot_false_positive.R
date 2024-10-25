library(RevGadgets)
library(tidyverse)
library(grid)
library(ggplot2)

fp <- tibble(sim = NA,
             theta_12 = 0, theta_13 = 0, theta_23 = 0,
             alpha_12 = 0, alpha_13 = 0, alpha_23 = 0,
             sigma2_12 = 0, sigma2_13 = 0, sigma2_23 = 0)

num_sim = 500
dir_in = "output/2_simulation/sx_logs/"
dir_out = "figures/2_simulation/triState/"
pars <- c("theta", "alpha", "sigma2")

for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_", i, "_sx_run_1_2mv_linkedPrior.log")
  df <- readTrace(file_in)[[1]] %>% 
    select(all_of(contains(c("theta_compare", "alpha_compare", "sigma2_compare"))))
  fp[i,] = c(i, lapply(df, mean))
}

fp_bin <- fp %>% 
  pivot_longer(contains("_"), names_to = "par")

fp_bin$par <- rep(c(rep("theta", 3), rep("alpha", 3), rep("sigma2", 3)), num_sim)

fp_bin$value <- ifelse(fp_bin$value > 0.975 | fp_bin$value < 0.025, 1, 0)
fp_bin <- fp_bin %>% 
  group_by(par) %>% 
  summarise(mean = mean(value))

thresholds <- tibble(sim=1:num_sim, lower=0.05)
p <- ggplot(fp_bin, aes(x=par)) +
  geom_point(aes(x=par, y=mean)) +
  theme_classic() +
  geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
  ylim(c(0, 0.10)) +
  ggtitle("False positive rates")
  #theme(axis.text.x = element_text(angle = 40, vjust = 0.8, hjust=0.9)) 

ggsave("figures/2_simulation/triState/fp_rates.pdf", p, width = 200, height = 120, units = "mm")


thresholds <- tibble(sim=1:num_sim, upper=0.975, lower=0.025)
for (par in pars){
  fp_par <- fp %>% 
    pivot_longer(contains(par), names_to = "par")
  
  p <- ggplot(fp_par) +
    geom_histogram(aes(x=value), binwidth = 0.05) +
    theme_classic() +
    geom_vline(mapping = aes(xintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
    geom_vline(mapping = aes(xintercept = upper), data=thresholds, linetype = "dashed", color = "brown") +
    ggtitle(par) +
    #scale_y_continuous(breaks = c(0, 20, 40)) +
    theme(axis.title = element_blank())
  
  file_out <- paste0(dir_out, "fp_hist_", par, ".pdf")
  ggsave(file_out, p, width = 100, height = 60, units = "mm")
}



