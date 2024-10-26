library(RevGadgets)
library(tidyverse)
library(ggplot2)

# plot true dthetas against P(state_0 > state_1)
#load(output/2_simulation/power_theta/grid.Rda")
num_sim_per_combo=200
num_combo=9
dir_in="output/2_simulation/power_theta/logs/"
dir_out="figures/2_simulation/power_theta/"
par_values <- read.csv("data/2_simulation/power_theta/pars.csv") %>% 
  mutate(stv = paste0(round(stv / 12.39784, 1), " v"),
         halflife = paste0(halflife, " TL"))
par_values$halflife = as.character(par_values$halflife)
grid <- expand.grid(sim=1:num_sim_per_combo, combo=1:num_combo)
grid$stv=par_values$stv
grid$halflife=par_values$halflife
grid$dtheta=par_values$theta_1 - par_values$theta_2

bar = txtProgressBar(style=3, width=40)
for (i in 1:nrow(grid)){
  filename <- paste0(dir_in, "sim_", i, "_run_1.log")
  df <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select("theta_compare_12") %>% 
    summarise(prob = mean(theta_compare_12))
  grid$prob[i] <- df$prob
  setTxtProgressBar(bar, i / nrow(grid))
}

thresholds <- tibble(sim=1:nrow(grid), upper=0.975, lower=0.025)

med_prob_theta <- grid %>%
  mutate(bin = round(dtheta, digits = 0)) %>% 
  group_by(combo, bin) %>% 
  summarise(med_prob=median(prob)) %>% 
  mutate(stv = paste0(combo / 2, " v"))
med_prob_theta$halflife = ifelse(med_prob_theta$combo %in% 1:3, "0.1 TL",
                                 ifelse(med_prob_theta$combo %in% 4:6, "0.3 TL",
                                        "0.6 TL"))

# plot individual reps as points
p1 <- ggplot(grid, aes(x=dtheta, y=prob)) +
  geom_point(aes(x = dtheta, y=prob), color="grey", alpha=0.5) +
  geom_point(data=med_prob_theta, aes(x = bin, y=med_prob, group = combo), color="black") +
  geom_line(data=med_prob_theta, aes(x = bin, y=med_prob, group=combo), color="black", linetype="dashed") +
  theme_bw() +
  labs(x = "theta1 - theta2", y = "P(theta1 > theta2)") + 
  geom_hline(mapping = aes(yintercept = upper), data=thresholds, linetype = "dashed", color = "brown") + 
  geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
  ylim(c(0, 1)) +
  facet_grid(rows = vars(stv), col = vars(halflife)) +
  ggtitle("Power of theta")
filename <- paste0(dir_out, "power_theta.pdf")
ggsave(filename, p1, width = 100, height = 100, units = "mm")

save(grid, file="output/2_simulation/power_theta/grid.Rda")

