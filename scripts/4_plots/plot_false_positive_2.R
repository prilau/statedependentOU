library(RevGadgets)
library(tidyverse)
library(grid)
library(ggplot2)

fp <- tibble(sim = NA,
             alpha_12 = 0, alpha_13 = 0, alpha_23 = 0,
             halflife_12 = 0, halflife_13 = 0, halflife_23 = 0,
             rho_12 = 0, rho_13 = 0, rho_23 = 0,
             sigma2_12 = 0, sigma2_13 = 0, sigma2_23 = 0,
             stv_12 = 0, stv_13 = 0, stv_23 = 0,
             theta_12 = 0, theta_13 = 0, theta_23 = 0)

num_sim = 50
dir_in = "output/2_simulation/sx_logs/"

for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_", i, "_sx_run_1_2mv_linkedPrior.log")
  df <- readTrace(file_in)[[1]] %>% 
  select(contains("compare"))
  fp[i,] = c(i, lapply(df, mean))
}

fp <- fp %>% 
  pivot_longer(contains("_"), names_to = "par") %>% 
  group_by(par)

fp_bin <- fp
fp_bin$value <- ifelse(fp$value > 0.975 | fp$value < 0.025, 1, 0)
fp_bin <- fp_bin %>% 
  group_by(par) %>% 
  summarise(mean = mean(value))

thresholds <- tibble(sim=1:num_sim, lower=0.05)

p <- ggplot(fp_bin, aes(x=par)) +
  geom_point(aes(x=par, y=mean)) +
  theme_classic() +
  geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
  ylim(c(0, 1)) +
  ggtitle("False positive rates") +
  theme(axis.text.x = element_text(angle = 40, vjust = 0.8, hjust=0.9)) 

ggsave("figures/2_simulation/triState/fp_rates.pdf", p, width = 200, height = 120, units = "mm")



thresholds <- tibble(sim=1:num_sim, upper=0.975, lower=0.025)

p <- ggplot(fp) +
  geom_histogram(aes(x=value), binwidth = 0.05) +
  theme_classic() +
  geom_vline(mapping = aes(xintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
  geom_vline(mapping = aes(xintercept = upper), data=thresholds, linetype = "dashed", color = "brown") +
  ggtitle("False positive rates") +
  facet_wrap(facets = vars(par), nrow = 6) +
  scale_y_continuous(breaks = c(0, 5, 10))
  
ggsave("figures/2_simulation/triState/fp_hist.pdf", p, width = 120, height = 200, units = "mm")






fn <- tibble(sim = NA,
             alpha_12 = 0, alpha_13 = 0, alpha_23 = 0,
             halflife_12 = 0, halflife_13 = 0, halflife_23 = 0,
             rho_12 = 0, rho_13 = 0, rho_23 = 0,
             sigma2_12 = 0, sigma2_13 = 0, sigma2_23 = 0,
             stv_12 = 0, stv_13 = 0, stv_23 = 0,
             theta_12 = 0, theta_13 = 0, theta_23 = 0)

num_sim = 50
dir_in = "output/2_simulation/sd_logs/"

for (i in 1:num_sim){
  file_in <- paste0(dir_in, "sim_", i, "_sd_run_1_2mv_linkedPrior.log")
  df <- readTrace(file_in)[[1]] %>% 
    select(contains("compare"))
  fn[i,] = c(i, lapply(df, mean))
}



# plot true dthetas against P(state_0 > state_1)
num_sim=50
dir_in="output/2_simulation/sd_logs/"
dir_out="figures/2_simulation/triState/"
par_values <- read.csv2("data/2_simulation/triState/pars.csv", sep=",")
pars <- c("theta", "alpha", "sigma2", "rho", "stv", "halflife")
df <- list()
for (i in 1:num_sim){
  prefix_sd <- paste0("sim_", i, "_sd")
  filename_sd <- paste0(dir_in, list.files(dir_in)[grep(prefix_sd, list.files(dir_in))])
  df[[i]] <- readTrace(filename_sd, burnin = 0.1)[[1]]
  df[[i]]$sim = i
}



dfx <- bind_rows(df[[1]], df[[2]], df[[3]], df[[4]], df[[5]]) %>%
  select(all_of(contains(c("compare", "sim"))))

ratios <- tibble(sim = 1:5)
ratios$theta =    pars_sd$theta_0    - pars_sd$theta_1
ratios$alpha =    log10(pars_sd$alpha_0    / pars_sd$alpha_1)
ratios$halflife = log10(pars_sd$halflife_0 / pars_sd$halflife_1)
ratios$sigma2 =   log10(pars_sd$sigma2_0   / pars_sd$sigma2_1)
ratios$stv =      log10(pars_sd$stv_0      / pars_sd$stv_1)
ratios$rho =      log10(pars_sd$rho_0      / pars_sd$rho_1)
thresholds <- tibble(sim=1:5, upper=0.975, lower=0.025)

for (par in pars){
  dfx_par <- pivot_longer(dfx, starts_with(par), names_to = "par")
  dfx_prob <- dfx_par %>%
    group_by(par, sim) %>%
    summarize(prob = mean(value))
  dfx_delta_prob <- tibble(par = ratios[[par]],
                           prob = dfx_prob$prob)
  
  
  p <- ggplot(dfx_delta_prob, aes(x=par, y=prob)) +
    geom_point(aes(x = par, y=prob)) +
    theme_classic() +
    labs(x = "log10(state_0 / state_1)", y = "P(state_0 > state_1)") + 
    geom_hline(mapping = aes(yintercept = upper), data=thresholds, linetype = "dashed", color = "brown") + 
    geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
    ylim(c(0, 1)) +
    ggtitle(par)
  filename <- paste0(dir_out, "p_greater_vs_delta_", par, ".pdf")
  ggsave(filename, p, width = 200, height = 120, units = "mm")
}








