library(RevGadgets)
library(tidyverse)
library(ggplot2)

# set up grid to summarise 
num_sim_per_combo=200
num_combo=9


##########################
# testing power of theta #
##########################

# uncomment the following line if there is an existing grid
#load("output/2_simulation/power_theta/grid.Rda")

grid <- expand.grid(sim=1:num_sim_per_combo, combo=1:num_combo)


# specify input and output directories
dir_in="output/2_simulation/power_theta/"
dir_out="figures/2_simulation/power_theta/"

# read 'true' parameters used for simulations
# v = empirical variance of log body size (kg) in mammals; ~= 12.39784
par_values <- read.csv("data/2_simulation/power_theta/pars.csv") %>% 
  mutate(stv = paste0(round(stv / 12.39784, 1), " v"),
         halflife = paste0(halflife, " TH"))
# add background OU parameters to grid
grid$stv=par_values$stv
grid$halflife=par_values$halflife
# add 'true' dtheta to grid
grid$dtheta=par_values$theta_1 - par_values$theta_2
# probability of theta_1 > theta_2
grid$prob=NA

bar = txtProgressBar(style=3, width=40)
# input P(theta_1 > theta_2) for each replicate into grid
for (i in 1:nrow(grid)){
  filename <- paste0(dir_in, "sim_", i, "_run_2.log")
  if (!file.exists(filename)) next
  grid$prob[i] <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select("theta_compare_12") %>% 
    summarise(prob = mean(theta_compare_12)) %>% 
    unlist() %>% 
    unname()
  setTxtProgressBar(bar, i / nrow(grid))
}
save(grid, file="output/2_simulation/power_theta/grid.Rda")

# remove the rows without input in prob
# this row is used when not all runs have been completed but you want plot a figure
grid <- grid  %>% 
  filter(!is.na(prob))

# for each combo of background parameters
# bin replicates by 'true' dthetas
# and calculate the mean P(theta_1 > theta_2) for each bin
mean_prob <- grid %>%
  mutate(bin = round(dtheta*2, digits = 0)/2) %>%  # bin size = 0.5
  group_by(combo, bin) %>% 
  summarise(mean_prob=mean(prob)) %>% 
  mutate(halflife = ifelse(combo %in% 1:3, "0.1 TH",
                           ifelse(combo %in% 4:6, "0.3 TH",
                                  "0.6 TH")),
         stv = ifelse(combo %in% c(1,4,7), "0.5 v",
                      ifelse(combo %in% c(2, 5, 8), "1 v",
                             "2 v")))


# for each combo
# bins between which probability intercepts with 0.025 or 0.975
intercepts <- mean_prob %>%
  group_by(combo) %>% 
  summarise(i025_1 = min(bin[which(mean_prob > 0.025)]),
            i025_2 = min(bin[which(mean_prob > 0.025)-1]),
            i975_1 = max(bin[which(mean_prob < 0.975)]),
            i975_2 = max(bin[which(mean_prob < 0.975)+1])) %>% 
  mutate(i025 = NA,
         i975 = NA)

for (i in 1:num_combo){
  x = c(intercepts$i025_1[i], intercepts$i025_2[i])
  y = mean_prob %>% filter(combo == i, bin %in% x) %>%
    ungroup() %>% 
    select(mean_prob) %>% 
    unlist()
  # skip if one of the bin exceed the boundary
  if (length(y)==2){
    eq = unname(lm(y~x)[[1]])
    intercepts$i025[i] = (0.025 - eq[1]) / eq[2]
  }
  
  x = c(intercepts$i975_1[i], intercepts$i975_2[i])
  y = mean_prob %>% filter(combo == i, bin %in% x) %>%
    ungroup() %>% 
    select(mean_prob) %>% 
    unlist()
  if (length(y)==2){
  eq = unname(lm(y~x)[[1]])
  intercepts$i975[i] = (0.975 - eq[1]) / eq[2]
  }
}

# calculate the width of dtheta within which 0.025 < mean(prob) < 0.975
intercepts <- intercepts %>% 
  mutate(width = i975-i025,
         halflife = ifelse(combo %in% 1:3, "0.1 TH",
                           ifelse(combo %in% 4:6, "0.3 TH",
                                  "0.6 TH")),
         stv = ifelse(combo %in% c(1,4,7), "0.5 v",
                      ifelse(combo %in% c(2, 5, 8), "1 v",
                             "2 v")))

# regression of prob~dtheta
grid$pred = NA
for (i in 1:num_combo){
  if (nrow(grid[which(grid$combo==i),]) == 0) next
  grid$pred[which(grid$combo==i)] <- predict(loess(prob~dtheta, grid[which(grid$combo==i),], span = 0.5))
}

# plot 1: elaborated plot where each point represents one replicate
thresholds <- tibble(sim=1:nrow(grid), upper=0.975, lower=0.025)
p1 <- ggplot(grid, aes(x=dtheta, y=prob)) +
  geom_point(aes(x = dtheta, y=prob), color="grey", alpha=0.5) +
  geom_point(data=mean_prob, aes(x = bin, y=mean_prob, group = combo), color="black") +
  #geom_line(data=mean_prob, aes(x = bin, y=mean_prob, group=combo), color="black", linetype="dashed") +
  geom_line(aes(y=pred), alpha=1.0) +
  theme_bw() +
  labs(x = "theta1 - theta2", y = "P(theta1 > theta2)") + 
  geom_hline(mapping = aes(yintercept = upper), data=thresholds, linetype = "dashed", color = "brown") + 
  geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
  ylim(c(0, 1)) +
  facet_grid(rows = vars(stv), col = vars(halflife)) +
  ggtitle("Power of theta")
filename <- paste0(dir_out, "power_theta_1.pdf")
ggsave(filename, p1, width = 100, height = 100, units = "mm")




##########################
# testing power of alpha #
##########################

# uncomment the following line if there is an existing grid
#load("output/2_simulation/power_alpha/grid.Rda")

grid <- expand.grid(sim=1:num_sim_per_combo, combo=1:num_combo)

# specify input and output directories
dir_in="output/2_simulation/power_alpha/"
dir_out="figures/2_simulation/power_alpha/"

# read 'true' parameters used for simulations
# v = empirical variance of log body size (kg) in mammals; ~= 12.39784
par_values <- read.csv("data/2_simulation/power_alpha/pars.csv") %>% 
  mutate(stv = paste0(round(stv / 12.39784, 1), " v"),
         dtheta_T = theta_1-theta_2)
# add background OU parameters to grid
grid$stv=par_values$stv
grid$dtheta_T=par_values$dtheta_T
# add 'true' dalpha and dhalflife to grid
grid$dalpha=par_values$alpha_1 - par_values$alpha_2
grid$dhalflife=log(2) / par_values$alpha_1 - log(2) / par_values$alpha_2
# add 'true' alphas ratio to grid
grid$ralpha=par_values$alpha_1 / par_values$alpha_2
grid$rhalflife= 1 / grid$ralpha
# probability of alpha_1 > alpha_2
grid$prob_a=NA
grid$prob_hl=NA

bar = txtProgressBar(style=3, width=40)
# input P(alpha_1 > alpha_2) for each replicate into grid
for (i in 1:nrow(grid)){
  filename <- paste0(dir_in, "sim_", i, "_run_1.log")
  if (!file.exists(filename)) next
  probs <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select(all_of(c("alpha_compare", "halflife_compare"))) %>% 
    summarise(prob_a = mean(alpha_compare),
              prob_hl = mean(halflife_compare)) %>% 
    unlist() %>% 
    unname()
  grid$prob_a[i] <- probs[1]
  grid$prob_hl[i] <- probs[2]
  setTxtProgressBar(bar, i / nrow(grid))
}
save(grid, file="output/2_simulation/power_alpha/grid.Rda")

# remove the rows without input in prob
# this row is used when not all runs have been completed but you want plot a figure
grid <- grid  %>% 
  filter(!is.na(prob))



###########################
# testing power of sigma2 #
###########################

grid <- expand.grid(sim=1:num_sim_per_combo, combo=1:num_combo)

# specify input and output directories
dir_in="output/2_simulation/power_sigma2/"
dir_out="figures/2_simulation/power_sigma2/"

par_values <- read.csv("data/2_simulation/power_sigma2/pars.csv") %>% 
  mutate(dstv_T = stv_1 - stv_2,
         dtheta_T = theta_1-theta_2,
         dsigma2_T = sigma2_1 - sigma2_2)
# add background OU d_stv to grid
grid$dstv_T = par_values$dstv_T
grid$halflife_T = par_values$halflife
grid$dtheta_T = par_values$dtheta_T
grid$dsigma2_T = par_values$dsigma2_T
grid$prob_stv = NA
grid$prob_sigma2 = NA

bar = txtProgressBar(style=3, width=40)
# input P(alpha_1 > alpha_2) for each replicate into grid
for (i in 1:nrow(grid)){
  filename <- paste0(dir_in, "sim_", i, "_run_1.log")
  if (!file.exists(filename)) next
  probs <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select(all_of(c("stv_compare", "sigma2_compare"))) %>% 
    summarise(prob_stv = mean(stv_compare),
              prob_sigma2 = mean(sigma2_compare)) %>% 
    unlist() %>% 
    unname()
  grid$prob_stv[i] <- probs[1]
  grid$prob_sigma2[i] <- probs[2]
  setTxtProgressBar(bar, i / nrow(grid))
}
save(grid, file="output/2_simulation/power_sigma2/grid.Rda")



#############
# new plots #
#############
powers_th <- grid  %>%
  mutate(sig = ifelse(prob >= 0.975 | prob <= 0.025, 1, 0),
         sig_T = ifelse(prob >= 0.975 & dtheta > 0 | prob <= 0.025 & dtheta < 0, 1, 0)) %>% 
  group_by(combo) %>% 
  summarise(power = mean(sig_T),
            ppv = sum(sig_T)/sum(sig),
            num_ppv = sum(sig)) %>% 
  mutate(halflife = ifelse(combo %in% 1:3, 0.1,
                           ifelse(combo %in% 4:6, 0.3,
                                  0.6)),
         stv = ifelse(combo %in% c(1,4,7), "0.5 v",
                      ifelse(combo %in% c(2, 5, 8), "1 v",
                             "2 v")),
         q025_power = qnorm(0.025, mean=power, sd = sqrt(power*(1-power)/200)), # sqrt(p(1-p)/n) is 1 standard error
         q975_power = qnorm(0.975, mean=power, sd = sqrt(power*(1-power)/200)),
         q025_ppv = qnorm(0.025, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)),
         q975_ppv = qnorm(0.975, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)))



powers_hl <- grid  %>%
  mutate(sig = ifelse(prob_hl >= 0.975 | prob_hl <= 0.025, 1, 0),
         sig_T = ifelse(prob_hl >= 0.975 & dhalflife > 0 | prob_hl <= 0.025 & dhalflife < 0, 1, 0)) %>% 
  group_by(combo) %>% 
  summarise(power = mean(sig_T),
            ppv = sum(sig_T)/sum(sig),
            num_ppv = sum(sig)) %>% 
  mutate(dtheta_T = ifelse(combo %in% 1:3, 2,
                           ifelse(combo %in% 4:6, 6,
                                  10)),
         stv = ifelse(combo %in% c(1,4,7), 0.5,
                      ifelse(combo %in% c(2, 5, 8), 1,
                             2)),
         q025_power = qnorm(0.025, mean=power, sd = sqrt(power*(1-power)/200)), # sqrt(p(1-p)/n) is 1 standard error
         q975_power = qnorm(0.975, mean=power, sd = sqrt(power*(1-power)/200)),
         q025_ppv = qnorm(0.025, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)),
         q975_ppv = qnorm(0.975, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)))

powers_stv <- grid  %>%
  mutate(sig = ifelse(prob_stv >= 0.975 | prob_stv <= 0.025, 1, 0),
         sig_T = ifelse(prob_stv >= 0.975 & dstv_T > 0 | prob_stv <= 0.025 & dstv_T < 0, 1, 0)) %>% 
  group_by(combo) %>% 
  summarise(power = mean(sig_T),
            ppv = sum(sig_T)/sum(sig),
            num_ppv = sum(sig)) %>% 
  mutate(dtheta_T = ifelse(combo %in% 1:3, 2,
                           ifelse(combo %in% 4:6, 6,
                                  10)),
         halflife = ifelse(combo %in% c(1,4,7), 0.1,
                      ifelse(combo %in% c(2, 5, 8), 0.3,
                             0.6)),
         q025_power = qnorm(0.025, mean=power, sd = sqrt(power*(1-power)/200)), # sqrt(p(1-p)/n) is 1 standard error
         q975_power = qnorm(0.975, mean=power, sd = sqrt(power*(1-power)/200)),
         q025_ppv = qnorm(0.025, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)),
         q975_ppv = qnorm(0.975, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)))


p_th_power <- ggplot(powers_th) +
  geom_point(aes(x = halflife, y=power, group=stv, shape=stv, color=stv), size=2) +
  geom_line(aes(x = halflife, y=power, group=stv, color=stv)) +
  geom_errorbar(aes(x = halflife, ymin = q025_power, ymax=q975_power, group=stv, color=stv))+
  scale_shape_manual(values=c(0, 1, 5))+
  theme_bw() +
  scale_x_continuous(breaks=c(0.1, 0.3, 0.6), label=c("0.1TH", "0.3TH", "0.6TH")) + 
  #labs(x = "", y = "theta") + 
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle("Power of theta")

p_th_ppv <- ggplot(powers_th) +
  geom_point(aes(x = halflife, y=ppv, group=stv, shape=stv, color=stv), size=2) +
  geom_line(aes(x = halflife, y=ppv, group=stv, color=stv)) +
  geom_errorbar(aes(x = halflife, ymin = q025_ppv, ymax=q975_ppv, group=stv, color=stv))+
  scale_shape_manual(values=c(0, 1, 5))+
  scale_x_continuous(breaks=c(0.1, 0.3, 0.6), label=c("0.1TH", "0.3TH", "0.6TH")) + 
  theme_bw() +
  theme(legend) +
  coord_cartesian(ylim=c(0, 1)) +
  #labs(x = "", y = "theta") + 
  ggtitle("PPV of theta")


p_hl_power <- ggplot(powers_hl) +
  geom_point(aes(x = stv, y=power, group=stv,shape=factor(dtheta_T), color=factor(dtheta_T)), size=2) +
  geom_line(aes(x = stv, y=power, group=factor(dtheta_T), color=factor(dtheta_T))) +
  theme_bw() +
  scale_shape_manual(values=c(0, 1, 5))+
  geom_errorbar(aes(x = stv, ymin = q025_power, ymax=q975_power, group=factor(dtheta_T), color=factor(dtheta_T)))+
  #labs(x = "", y = "theta") + 
  scale_x_continuous(breaks=c(0.5, 1, 2), label=c("0.5v", "v", "2v")) + 
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle("Power of halflife")


p_hl_ppv <- ggplot(powers_hl) +
  geom_point(aes(x = stv, y=ppv, group=stv,shape=factor(dtheta_T), color=factor(dtheta_T)), size=2) +
  geom_line(aes(x = stv, y=ppv, group=factor(dtheta_T), color=factor(dtheta_T))) +
  theme_bw() +
  geom_errorbar(aes(x = stv, ymin = q025_ppv, ymax=q975_ppv, group=factor(dtheta_T), color=factor(dtheta_T)))+
  scale_shape_manual(values=c(0, 1, 5)) +
  scale_x_continuous(breaks=c(0.5, 1, 2), label=c("0.5v", "v", "2v")) + 
  coord_cartesian(ylim=c(0, 1)) +
  #labs(x = "", y = "theta") + 
  ggtitle("PPV of halflife")

p_stv_power <- ggplot(powers_stv) +
  geom_point(aes(x = dtheta_T, y=power, group=factor(halflife),shape=factor(halflife), color=factor(halflife)), size=2) +
  geom_line(aes(x = dtheta_T, y=power, group=factor(halflife), color=factor(halflife))) +
  theme_bw() +
  scale_shape_manual(values=c(0, 1, 5))+
  geom_errorbar(aes(x = dtheta_T, ymin = q025_power, ymax=q975_power, group=halflife, color=factor(halflife)))+
  #labs(x = "", y = "theta") + 
  scale_x_continuous(breaks=c(2, 6, 10)) + 
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle("Power of stationary variance")

p_stv_ppv <- ggplot(powers_stv) +
  geom_point(aes(x = dtheta_T, y=ppv, group=factor(halflife), shape=factor(halflife), color=factor(halflife)), size=2) +
  geom_line(aes(x = dtheta_T, y=ppv, group=factor(halflife), color=factor(halflife))) +
  theme_bw() +
  geom_errorbar(aes(x = dtheta_T, ymin = q025_ppv, ymax=q975_ppv, group=factor(halflife), color=factor(halflife)))+
  scale_shape_manual(values=c(0, 1, 5)) +
  
  scale_x_continuous(breaks=c(2, 6, 10)) + 
  coord_cartesian(ylim=c(0, 1)) +
  #labs(x = "", y = "theta") + 
  ggtitle("PPV of stationary variance")

grid.arrange(p_th_power, p_hl_power, p_stv_power, p_th_ppv, p_hl_ppv, p_stv_ppv, nrow = 2)
ggsave("figures/2_simulation/power_theta/power_ppv.pdf", width = 200, height = 180, unit = "mm")
