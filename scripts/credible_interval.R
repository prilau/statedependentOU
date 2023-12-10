library(RevGadgets)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(tikzDevice)


# true values
true_halfLife <- 0.35
true_alpha <- log(2)/true_halfLife
true_stV <- 0.0625
true_sigma2 <- true_stV * 2 * true_alpha
true_theta0 <- 0.5
true_theta1 <- 2.0


##################
## Simulation 2 ##
##################

num_tips = c(100, 250, 500)
num_dtraits = 4
num_ctraits = 4

sim2_grid = expand.grid(tips = num_tips, dtrait=1:num_dtraits,
                   ctrait=1:num_ctraits, stringsAsFactors=FALSE)

p_CI_100 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

p_CI_250 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

p_CI_500 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_100 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_250 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_500 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))


bar = txtProgressBar(style=3, width=40)
for(i in 1:nrow(sim2_grid)) {
  
  this_row = sim2_grid[i,]
  this_num_tips = this_row[[1]]
  this_num_dtraits = this_row[[2]]
  this_num_ctraits = this_row[[3]]
  
  this_log_file = paste0('output/sdOU_simulation2_n', this_num_tips,
                         't1', 'd', this_num_dtraits, 'c', this_num_ctraits,
                         '.log')
  
  trace_quant <- readTrace(path = this_log_file, burnin = 1/7)

  CI_halfLife <- quantile(trace_quant[[1]]$halfLife, c(0.025, 0.975))
  CI_alpha <- quantile(trace_quant[[1]]$`alpha[1]`, c(0.025, 0.975))
  CI_stV <- quantile(trace_quant[[1]]$`stationaryVariance[1]`, c(0.025, 0.975))
  CI_sigma2 <- quantile(trace_quant[[1]]$`sigma2[1]`, c(0.025, 0.975))
  CI_theta0 <- quantile(trace_quant[[1]]$`theta[1]`, c(0.025, 0.975))
  CI_theta1 <- quantile(trace_quant[[1]]$`theta[2]`, c(0.025, 0.975))
  if (this_num_tips == 100) {
    r_CI_100$value[1] <- r_CI_100$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_100$value[2] <- r_CI_100$value[2] + CI_alpha[2] - CI_alpha[1]
    r_CI_100$value[3] <- r_CI_100$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_100$value[4] <- r_CI_100$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_100$value[5] <- r_CI_100$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_100$value[6] <- r_CI_100$value[6] + CI_theta1[2] - CI_theta1[1]

    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_100$value[1] <- p_CI_100$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_100$value[2] <- p_CI_100$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_100$value[3] <- p_CI_100$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_100$value[4] <- p_CI_100$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_100$value[5] <- p_CI_100$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_100$value[6] <- p_CI_100$value[6] + 1
    }
  }
  else if (this_num_tips == 250) {
    r_CI_250$value[1] <- r_CI_250$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_250$value[2] <- r_CI_250$value[2] + CI_halfLife[2] - CI_alpha[1]
    r_CI_250$value[3] <- r_CI_250$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_250$value[4] <- r_CI_250$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_250$value[5] <- r_CI_250$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_250$value[6] <- r_CI_250$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_250$value[1] <- p_CI_250$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_250$value[2] <- p_CI_250$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_250$value[3] <- p_CI_250$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_250$value[4] <- p_CI_250$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_250$value[5] <- p_CI_250$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_250$value[6] <- p_CI_250$value[6] + 1
    }
  }
  else if (this_num_tips == 500) {
    r_CI_500$value[1] <- r_CI_500$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_500$value[2] <- r_CI_500$value[2] + CI_halfLife[2] - CI_alpha[1]
    r_CI_500$value[3] <- r_CI_500$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_500$value[4] <- r_CI_500$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_500$value[5] <- r_CI_500$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_500$value[6] <- r_CI_500$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_500$value[1] <- p_CI_500$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_500$value[2] <- p_CI_500$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_500$value[3] <- p_CI_500$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_500$value[4] <- p_CI_500$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_500$value[5] <- p_CI_500$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_500$value[6] <- p_CI_500$value[6] + 1
    }
  }
  setTxtProgressBar(bar, i / nrow(sim2_grid))
}

p_CI_100 <- p_CI_100 %>% 
  mutate(value = value/16)

p_CI_250 <- p_CI_250 %>% 
  mutate(value = value/16)

p_CI_500 <- p_CI_500 %>% 
  mutate(value = value/16)

r_CI_100 <- r_CI_100 %>% 
  mutate(value = value/16)

r_CI_250 <- r_CI_250 %>% 
  mutate(value = value/16)

r_CI_500 <- r_CI_500 %>% 
  mutate(value = value/16)




##################
## Simulation 3 ##
##################
rate = c(5, 10, 20, 50, 500)
num_dtraits = 4
num_ctraits = 4

sim3_grid = expand.grid(rate = rate, dtrait=1:num_dtraits,
                   ctrait=1:num_ctraits, stringsAsFactors=FALSE)



p_CI_r5 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))
p_CI_r10 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))
p_CI_r20 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))
p_CI_r50 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))
p_CI_r500 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_r5 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_r10 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_r20 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_r50 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

r_CI_r500 <- tibble(
  parameter = c("halfLife", "alpha", "stV",
                "sigma2", "theta0", "theta1"),
  value = c(0, 0, 0, 0, 0, 0))

bar = txtProgressBar(style=3, width=40)

for(i in 1:nrow(sim3_grid)) {
  this_row = sim3_grid[i,]
  this_rate = this_row[[1]]
  this_num_dtraits = this_row[[2]]
  this_num_ctraits = this_row[[3]]
  
  this_log_file = paste0('output/sdOU_simulation3_n500t1r', this_rate,
                         'd', this_num_dtraits, 'c', this_num_ctraits,
                         '.log')
  
  trace_quant <- readTrace(path = this_log_file, burnin = 1/7)
  
  CI_halfLife <- quantile(trace_quant[[1]]$halfLife, c(0.025, 0.975))
  CI_alpha <- quantile(trace_quant[[1]]$`alpha[1]`, c(0.025, 0.975))
  CI_stV <- quantile(trace_quant[[1]]$`stationaryVariance[1]`, c(0.025, 0.975))
  CI_sigma2 <- quantile(trace_quant[[1]]$`sigma2[1]`, c(0.025, 0.975))
  CI_theta0 <- quantile(trace_quant[[1]]$`theta[1]`, c(0.025, 0.975))
  CI_theta1 <- quantile(trace_quant[[1]]$`theta[2]`, c(0.025, 0.975))

  if (this_rate == 5) {
    r_CI_r5$value[1] <- r_CI_r5$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_r5$value[2] <- r_CI_r5$value[2] + CI_alpha[2] - CI_alpha[1]
    r_CI_r5$value[3] <- r_CI_r5$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_r5$value[4] <- r_CI_r5$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_r5$value[5] <- r_CI_r5$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_r5$value[6] <- r_CI_r5$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_r5$value[1] <- p_CI_r5$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_r5$value[2] <- p_CI_r5$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_r5$value[3] <- p_CI_r5$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_r5$value[4] <- p_CI_r5$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_r5$value[5] <- p_CI_r5$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_r5$value[6] <- p_CI_r5$value[6] + 1
    }
  }
  else if (this_rate == 10) {
    r_CI_r10$value[1] <- r_CI_r10$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_r10$value[2] <- r_CI_r10$value[2] + CI_alpha[2] - CI_alpha[1]
    r_CI_r10$value[3] <- r_CI_r10$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_r10$value[4] <- r_CI_r10$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_r10$value[5] <- r_CI_r10$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_r10$value[6] <- r_CI_r10$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_r10$value[1] <- p_CI_r10$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_r10$value[2] <- p_CI_r10$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_r10$value[3] <- p_CI_r10$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_r10$value[4] <- p_CI_r10$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_r10$value[5] <- p_CI_r10$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_r10$value[6] <- p_CI_r10$value[6] + 1
    }
  }
  else if (this_rate == 20) {
    r_CI_r20$value[1] <- r_CI_r20$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_r20$value[2] <- r_CI_r20$value[2] + CI_alpha[2] - CI_alpha[1]
    r_CI_r20$value[3] <- r_CI_r20$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_r20$value[4] <- r_CI_r20$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_r20$value[5] <- r_CI_r20$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_r20$value[6] <- r_CI_r20$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_r20$value[1] <- p_CI_r20$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_r20$value[2] <- p_CI_r20$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_r20$value[3] <- p_CI_r20$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_r20$value[4] <- p_CI_r20$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_r20$value[5] <- p_CI_r20$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_r20$value[6] <- p_CI_r20$value[6] + 1
    }
  }
  else if (this_rate == 50) {
    r_CI_r50$value[1] <- r_CI_r50$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_r50$value[2] <- r_CI_r50$value[2] + CI_alpha[2] - CI_alpha[1]
    r_CI_r50$value[3] <- r_CI_r50$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_r50$value[4] <- r_CI_r50$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_r50$value[5] <- r_CI_r50$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_r50$value[6] <- r_CI_r50$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_r50$value[1] <- p_CI_r50$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_r50$value[2] <- p_CI_r50$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_r50$value[3] <- p_CI_r50$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_r50$value[4] <- p_CI_r50$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_r50$value[5] <- p_CI_r50$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_r50$value[6] <- p_CI_r50$value[6] + 1
    }
  }
  else if (this_rate == 500) {
    r_CI_r500$value[1] <- r_CI_r500$value[1] + CI_halfLife[2] - CI_halfLife[1]
    r_CI_r500$value[2] <- r_CI_r500$value[2] + CI_alpha[2] - CI_alpha[1]
    r_CI_r500$value[3] <- r_CI_r500$value[3] + CI_stV[2] - CI_stV[1]
    r_CI_r500$value[4] <- r_CI_r500$value[4] + CI_sigma2[2] - CI_sigma2[1]
    r_CI_r500$value[5] <- r_CI_r500$value[5] + CI_theta0[2] - CI_theta0[1]
    r_CI_r500$value[6] <- r_CI_r500$value[6] + CI_theta1[2] - CI_theta1[1]
    if (true_halfLife > CI_halfLife[1] && true_halfLife < CI_halfLife[2]) {
      p_CI_r500$value[1] <- p_CI_r500$value[1] + 1
    }
    if (true_alpha > CI_alpha[1] && true_alpha < CI_alpha[2]) {
      p_CI_r500$value[2] <- p_CI_r500$value[2] + 1
    }
    if (true_stV > CI_stV[1] && true_stV < CI_stV[2]) {
      p_CI_r500$value[3] <- p_CI_r500$value[3] + 1
    }
    if (true_sigma2 > CI_sigma2[1] && true_sigma2 < CI_sigma2[2]) {
      p_CI_r500$value[4] <- p_CI_r500$value[4] + 1
    }
    if (true_theta0 > CI_theta0[1] && true_theta0 < CI_theta0[2]) {
      p_CI_r500$value[5] <- p_CI_r500$value[5] + 1
    }
    if (true_theta1 > CI_theta1[1] && true_theta1 < CI_theta1[2]) {
      p_CI_r500$value[6] <- p_CI_r500$value[6] + 1
    }
  }
  setTxtProgressBar(bar, i / nrow(sim3_grid))
}

p_CI_r5 <- p_CI_r5 %>% 
  mutate(value = value/16)
p_CI_r10 <- p_CI_r10 %>% 
  mutate(value = value/16)
p_CI_r20 <- p_CI_r20 %>% 
  mutate(value = value/16)
p_CI_r50 <- p_CI_r50 %>% 
  mutate(value = value/16)
p_CI_r500 <- p_CI_r500 %>% 
  mutate(value = value/16)

r_CI_r5 <- r_CI_r5 %>% 
  mutate(value = value/16)
r_CI_r10 <- r_CI_r10 %>% 
  mutate(value = value/16)
r_CI_r20 <- r_CI_r20 %>% 
  mutate(value = value/16)
r_CI_r50 <- r_CI_r50 %>% 
  mutate(value = value/16)
r_CI_r500 <- r_CI_r500 %>% 
  mutate(value = value/16)

#################
## figure time ##
#################

sim3_r <- bind_rows(list(r5 = r_CI_r5, r10 = r_CI_r10, r20 = r_CI_r20,
                         r50 = r_CI_r50, r500 = r_CI_r500), .id = "rate")
sim3_p <- bind_rows(list(r5 = p_CI_r5, r10 = p_CI_r10, r20 = p_CI_r20,
                         r50 = p_CI_r50, r500 = p_CI_r500), .id = "rate")
sim2_r <- bind_rows(list(n100 = r_CI_100, n250 = r_CI_250, n500 = r_CI_500),
                    .id = "num_tips")
sim2_p <- bind_rows(list(n100 = p_CI_100, n250 = p_CI_250, n500 = p_CI_500),
                    .id = "num_tips")
parameter_color <- c('#6699CC', '#004488', '#EECC66', '#994455', '#997700', '#EE99AA')
names(parameter_color) <- levels(sim3_r$parameter)
parameter_shape <- c(0, 1, 2, 15, 16, 17)
names(parameter_shape) <- levels(sim3_r$parameter)


sim3_level_order <- c('r5', 'r10', 'r20', 'r50', 'r500') 

sim3_rangePlot <- sim3_r %>%
  ggplot(aes(x=factor(rate, level = sim3_level_order), y=log(value),
             group = parameter, col = parameter, shape = parameter)) +
  geom_line() +
  geom_point() +
  xlab("Rate") +
  ylab("$\\ln$ (size)") +
  scale_colour_manual(name = "Parameter",
                      labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                 "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                      values = parameter_color) +
  scale_shape_manual(name = "Parameter",
                     labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                     values = parameter_shape) +
  theme_bw() +
  ggtitle("Size of 95% credible interval")




sim3_percPlot <- sim3_p %>%
  ggplot(aes(x=factor(rate, level = sim3_level_order), y=value,
             group = parameter, col = parameter, shape = parameter)) +
  geom_line() +
  geom_point() +
  xlab("Rate") +
  ylab("") +
  scale_colour_manual(name = "Parameter",
                      labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                 "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                      values = parameter_color) +
  scale_shape_manual(name = "Parameter",
                     labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                     values = parameter_shape) + 
  theme_bw() +
  ggtitle("Frequency of true parameter values falling within $\\newline$ 95% credible interval")


nested_sim3 <- (sim3_percPlot + sim3_rangePlot) +
  plot_annotation(tag_levels = 'A')

tikzDevice::tikz(file = "figures/sim3.tex", width = 5, height = 3)
nested_sim3
dev.off()



sim2_level_order <- c('n100', 'n250', 'n500') 


sim2_rangePlot <- sim2_r %>%
  ggplot(aes(x=factor(num_tips, level = sim2_level_order), y=log(value),
             group = parameter, col = parameter, shape = parameter)) +
  geom_line() +
  geom_point() +
  xlab("Number of tips") +
  ylab("$\\ln$ (size)") +
  scale_colour_manual(name = "Parameter",
                      labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                 "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                      values = parameter_color) +
  scale_shape_manual(name = "Parameter",
                     labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                     values = parameter_shape) +   theme_bw() +
  ggtitle("Size of 95% credible interval")

sim2_percPlot <- sim2_p %>%
  ggplot(aes(x=factor(num_tips, level = sim2_level_order), y=value,
             group = parameter, col = parameter, shape = parameter)) +
  geom_line() +
  geom_point() +
  xlab("Number of tips") +
  ylab("") +
  scale_colour_manual(name = "Parameter",
                      labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                 "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                      values = parameter_color) +
  scale_shape_manual(name = "Parameter",
                     labels = c("$\\alpha$", "$t_{1/2}$", "$\\sigma^2$",
                                "$V_{st}$", "$\\theta_0$", "$\\theta_1$"),
                     values = parameter_shape) +   ggtitle("Frequency of true values falling within $\\newline$ 95% credible interval") +
  theme_bw()


nested_sim2 <- (sim2_percPlot + sim2_rangePlot) +
  plot_annotation(tag_levels = 'A')

tikzDevice::tikz(file = "figures/sim2.tex", width = 5, height = 3)
nested_sim2
dev.off()


ggsave("figures/sims_posterior_summary.pdf", nested, width = 400, height = 300, units = "mm")

