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



# plot 2: concise plot where each point represents a background parameter combo
# width of dtheta where 0.025 < P(i>j) < 0.975 is used as proxy
# of false negative rate
p2 <- ggplot(intercepts, aes(x=halflife, y=width)) +
  geom_point(aes(x = halflife, y=width, group=stv, color=stv)) +
  geom_line(aes(x = halflife, y=width, group=stv, color=stv)) +
  theme_bw() +
  #labs(x = "dtheta", y = "Phylogenetic half-life") + 
  ggtitle("Width of dtheta")
filename <- paste0(dir_out, "power_theta_2.pdf")
ggsave(filename, p2, width = 100, height = 100, units = "mm")




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
         theta = paste0("{", theta_1, theta_2, "}"))
# add background OU parameters to grid
grid$stv=par_values$stv
grid$theta=par_values$theta
# add 'true' dalpha and dhalflife to grid
grid$dalpha=par_values$alpha_1 - par_values$alpha_2
grid$dhalflife=log(2) / grid$dalpha
# add 'true' alphas ratio to grid
grid$ralpha=par_values$alpha_1 / par_values$alpha_2
grid$rhalflife= 1 / grid$ralpha
# probability of alpha_1 > alpha_2
grid$prob=NA

bar = txtProgressBar(style=3, width=40)
# input P(alpha_1 > alpha_2) for each replicate into grid
for (i in 1:nrow(grid)){
  filename <- paste0(dir_in, "sim_", i, "_run_1.log")
  if (!file.exists(filename)) next
  grid$prob[i] <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select("theta_compare_12") %>% 
    summarise(prob = mean(theta_compare_12)) %>% 
    unlist() %>% 
    unname()
  setTxtProgressBar(bar, i / nrow(grid))
}
save(grid, file="output/2_simulation/power_alpha/grid.Rda")

# remove the rows without input in prob
# this row is used when not all runs have been completed but you want plot a figure
grid <- grid  %>% 
  filter(!is.na(prob))

# for each combo of background parameters
# bin replicates by 'true' dthetas
# and calculate the mean P(theta_1 > theta_2) for each bin
mean_prob <- grid %>%
  mutate(bin = round(dalpha*2, digits = 0)/2) %>%  # bin size = 0.5
  group_by(combo, bin) %>% 
  summarise(mean_prob=mean(prob)) %>% 
  mutate(theta = ifelse(combo %in% 1:3, "{1,-1}",
                           ifelse(combo %in% 4:6, "{3,-3}",
                                  "{5,-5}")),
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
  y = mean_prob_theta %>% filter(combo == i, bin %in% x) %>%
    ungroup() %>% 
    select(mean_prob) %>% 
    unlist()
  # skip if one of the bin exceed the boundary
  if (length(y)==2){
    eq = unname(lm(y~x)[[1]])
    intercepts$i025[i] = (0.025 - eq[1]) / eq[2]
  }
  
  x = c(intercepts$i975_1[i], intercepts$i975_2[i])
  y = mean_prob_theta %>% filter(combo == i, bin %in% x) %>%
    ungroup() %>% 
    select(mean_prob) %>% 
    unlist()
  if (length(y)==2){
    eq = unname(lm(y~x)[[1]])
    intercepts$i975[i] = (0.975 - eq[1]) / eq[2]
  }
}

# calculate the width of dalpha within which 0.025 < mean(prob) < 0.975
intercepts <- intercepts %>% 
  mutate(width = i975-i025,
         theta = ifelse(combo %in% 1:3, "{1,-1}",
                           ifelse(combo %in% 4:6, "{3,-3}",
                                  "{5,-5}")),
         stv = ifelse(combo %in% c(1,4,7), "0.5 v",
                      ifelse(combo %in% c(2, 5, 8), "1 v",
                             "2 v")))

# plot 1: elaborated plot where each point represents one replicate
thresholds <- tibble(sim=1:nrow(grid), upper=0.975, lower=0.025)
p1 <- ggplot(grid, aes(x=dalpha, y=prob)) +
  geom_point(aes(x = dalpha, y=prob), color="grey", alpha=0.5) +
  geom_point(data=mean_prob, aes(x = bin, y=mean_prob, group = combo), color="black") +
  geom_line(data=mean_prob, aes(x = bin, y=mean_prob, group=combo), color="black", linetype="dashed") +
  theme_bw() +
  labs(x = "alpha1 - alpha2", y = "P(alpha1 > alpha2)") + 
  geom_hline(mapping = aes(yintercept = upper), data=thresholds, linetype = "dashed", color = "brown") + 
  geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
  ylim(c(0, 1)) +
  facet_grid(rows = vars(stv), col = vars(theta)) +
  ggtitle("Power of alpha")
filename <- paste0(dir_out, "power_alpha_1.pdf")
ggsave(filename, p1, width = 100, height = 100, units = "mm")

# plot 2: concise plot where each point represents a background parameter combo
# width of dalpha where 0.025 < P(i>j) < 0.975 is used as proxy
# of false negative rate
p2 <- ggplot(intercepts, aes(x=alpha, y=width)) +
  geom_point(aes(x = stv, y=width, group=theta, color=theta)) +
  geom_line(aes(x = stv, y=width, group=theta, color=theta)) +
  theme_bw() +
  #labs(x = "", y = "theta") + 
  ggtitle("Width of dalpha")
filename <- paste0(dir_out, "power_alpha_2.pdf")
ggsave(filename, p2, width = 100, height = 100, units = "mm")


###########################
# testing power of sigma2 #
###########################
