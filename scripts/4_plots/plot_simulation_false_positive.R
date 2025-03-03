library(RevGadgets)
library(tidyverse)
library(grid)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(patchwork)
library(latex2exp)

fp <- tibble(sim = NA,
             theta = NA,
             alpha = NA,
             sigma2 = NA,
             halflife = NA,
             stv = NA)

num_sim = 1000
dir_in = "output/2_simulation/false_positive/"
dir_out = "figures/2_simulation/"
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
pars <- c("halflife", "stv", "theta")
pars_title <- c(TeX("Phylogenetic half-life $t_{0.5}$"), TeX("Stationary variance $V_y$"), TeX("Optimum $\\theta$"))
pars_ylab <- c(TeX("P($t_{0.5\\,0} > t_{0.5\\,1}$)"), TeX("P($V_{y\\,0} > V_{y\\,1}$)"), TeX("P($\\theta_0 > \\theta_1$)"))

i=1
for (par in pars){
  fp_par <- fp %>% 
    pivot_longer(contains(par), names_to = "par")
  
  p[[par]] <- ggplot(fp_par) +
    geom_histogram(aes(x=value), binwidth = 0.05) +
    theme_classic() +
    geom_vline(mapping = aes(xintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
    geom_vline(mapping = aes(xintercept = upper), data=thresholds, linetype = "dashed", color = "brown") +
    ggtitle(pars_title[i]) +
    scale_y_continuous(breaks=c(0, 30, 60, 90)) + 
    coord_cartesian(ylim=c(0, 105)) +
    xlab(pars_ylab[i]) +
    ylab("") +
    #scale_y_continuous(breaks = c(0, 20, 40)) +
    theme(plot.title = element_text(hjust = 0.5, size=14),
          legend.position = "bottom",
          axis.text = element_text(size=12),
          legend.title = element_text(size=12),
          legend.text = element_text(size=12),
          axis.title.x = element_text(size=12),
          axis.title.y = element_blank())
  i=i+1
}

# plot table
fp_bin <- fp %>% 
  pivot_longer(!c("sim"), names_to = "par") %>% 
  filter(par %in% c("halflife", "stv", "theta")) %>% 
  mutate(value = ifelse(value > 0.975 | value < 0.025, 1, 0)) %>% 
  group_by(par) %>%
  summarise(`False positive rate` = signif(mean(value),2)) %>% 
  rename(Parameter = par)

p_all <- cowplot::plot_grid(p[[1]], p[[2]], p[[3]], nrow = 1)

file_out <- paste0(dir_out, "/fpr.pdf")
ggsave(file_out, p_all, width = 180, height = 70, units = "mm")

