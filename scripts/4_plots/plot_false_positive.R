library(RevGadgets)
library(tidyverse)
library(grid)

plot_dthetas <- function(df, med_true, dir_out, num_sim){
  p <- ggplot(df, aes(x = dthetas, fill = state)) +
    geom_density(aes(x = dthetas, after_stat(density), fill=state), color = "black", alpha = 0.4) +
    theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.8, 0.8)) +
    geom_vline(mapping = aes(xintercept = med, color=state), data = med_true, linetype = "dashed") + 
    geom_vline(mapping = aes(xintercept = true), data = med_true, color = c("darkred", "darkblue")) + 
    labs(x = "dthetas", y = "posterior density")
  filename <- paste0(dir_out, "sim_", i, "_dthetas.pdf")
  ggsave(filename, p, width = 200, height = 120, units = "mm")
}




# plot true dthetas, sd vs stateless
num_sim=5
dir_in="output/2_simulation/binState/logs/"
dir_out="figures/2_simulation/binState/"
load("data/2_simulation/binState/par_values.Rda")

for (i in 1:num_sim){
  prefix_sd <- paste0("sim_", i, "_sd")
  prefix_sx <- paste0("sim_", i, "_stateless")
  filename_sd <- paste0(dir_in, list.files(dir_in)[grep(prefix_sd, list.files(dir_in))])
  filename_sx <- paste0(dir_in, list.files(dir_in)[grep(prefix_sx, list.files(dir_in))])
  
  df <- readTrace(filename_sd, burnin = 0.1)
  df[[2]] <- readTrace(filename_sx, burnin = 0.1)[[1]]
  
  df[[1]]$state <- "sd"
  df[[2]]$state <- "sx"
  
  dfx <- bind_rows(df[[1]], df[[2]]) %>%
    select(c(dthetas, state))
  
  df_med_true <- dfx %>%
    group_by(state) %>%
    summarize(med = median(dthetas))
  df_med_true$true <- c(par_values$theta_0[i] - par_values$theta_1[i], 0)
  
  if (!dir.exists(dir_out)){
    dir.create(dir_out, showWarnings = FALSE)
  }
  plot_dthetas(dfx, df_med_true, dir_out, num_sim = i)
}



# plot true dthetas against P(state_0 > state_1)
num_sim=100
dir_in="output/2_simulation/sd_logs/"
dir_out="figures/2_simulation/triState/"
par_values <- read.csv("data/2_simulation/triState/pars.csv") %>% 
  filter(sim < 101,
         state == "sd")



pars <- c("theta", "alpha", "sigma2", "rho", "stv", "halflife")


for (i in 1:num_sim){
  filename_sd <- paste0(dir_in, "sim_", i, "_sd_run_1_2mv_linkedPrior.log")
  #filename_sd <- paste0(dir_in, list.files(dir_in)[grep(prefix_sd, list.files(dir_in))])
  df <- readTrace(filename_sd, burnin = 0.1)[[1]]
  df$sim = i
  if (i == 1){
    dfx <- df
  } else {
    dfx <- bind_rows(dfx, df)
  }
}



dfx <- dfx %>%
  select(all_of(contains(c("compare", "sim"))))

ratios <- tibble(sim = 1:100)
ratios$theta12 =    par_values$theta_1    - par_values$theta_2
ratios$alpha12 =    log10(par_values$alpha_1    / par_values$alpha_2)
ratios$halflife12 = log10(par_values$halflife_1 / par_values$halflife_2)
ratios$sigma212 =   log10(par_values$sigma2_1   / par_values$sigma2_2)
ratios$stv12 =      log10(par_values$stv_1      / par_values$stv_2)
ratios$rho12 =      log10(par_values$rho_1      / par_values$rho_2)

ratios$theta13 =    par_values$theta_1    - par_values$theta_3
ratios$alpha13 =    log10(par_values$alpha_1    / par_values$alpha_3)
ratios$halflife13 = log10(par_values$halflife_1 / par_values$halflife_3)
ratios$sigma213 =   log10(par_values$sigma2_1   / par_values$sigma2_3)
ratios$stv13 =      log10(par_values$stv_1      / par_values$stv_3)
ratios$rho13 =      log10(par_values$rho_1      / par_values$rho_3)

ratios$theta23 =    par_values$theta_2    - par_values$theta_3
ratios$alpha23 =    log10(par_values$alpha_2    / par_values$alpha_3)
ratios$halflife23 = log10(par_values$halflife_2 / par_values$halflife_3)
ratios$sigma223 =   log10(par_values$sigma2_2   / par_values$sigma2_3)
ratios$stv23 =      log10(par_values$stv_2      / par_values$stv_3)
ratios$rho23 =      log10(par_values$rho_2      / par_values$rho_3)


thresholds <- tibble(sim=1:100, upper=0.975, lower=0.025)

for (par in pars){
  dfx_par <- pivot_longer(dfx, starts_with(par), names_to = "par")
  dfx_prob <- dfx_par %>%
    group_by(par, sim) %>%
    summarize(prob = mean(value))
  
  par_12 <- paste0(par, "12")
  par_13 <- paste0(par, "13")
  par_23 <- paste0(par, "23")
  dfx_delta_prob <- tibble(par = c(ratios[[par_12]], ratios[[par_13]], ratios[[par_23]]),
                           prob = dfx_prob$prob)
  
  # plot individual reps as points
  p1 <- ggplot(dfx_delta_prob, aes(x=par, y=prob)) +
    geom_point(aes(x = par, y=prob)) +
    theme_classic() +
    labs(x = "log10(state_i / state_j)", y = "P(state_i > state_j)") + 
    geom_hline(mapping = aes(yintercept = upper), data=thresholds, linetype = "dashed", color = "brown") + 
    geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
    ylim(c(0, 1)) +
    ggtitle(par)
  filename <- paste0(dir_out, "power_points_", par, ".pdf")
  ggsave(filename, p1, width = 200, height = 120, units = "mm")
  
  
  # plot rolling medians
  dfx_bin <- dfx_delta_prob %>% 
    mutate(bin = round(par/2, digits = 1)*2) %>% 
    group_by(bin)
  
  dfx_mean_bin <- dfx_bin %>% 
    summarise(mean_prob = mean(prob))
  
  p2 <- ggplot(dfx_bin, aes(x=bin, y=prob)) +
    geom_boxplot(aes(x=bin, y=prob, group=bin)) +
    #geom_line(data=dfx_mean_bin, aes(x=bin, y=mean_prob), color="blue") +
    theme_classic() +
    labs(x = "log10(state_i / state_j)", y = "P(state_i > state_j)") + 
    geom_hline(mapping = aes(yintercept = upper), data=thresholds, linetype = "dashed", color = "brown") + 
    geom_hline(mapping = aes(yintercept = lower), data=thresholds, linetype = "dashed", color = "brown") +
    ylim(c(0, 1)) +
    ggtitle(par)
  
  filename <- paste0(dir_out, "power_binned_", par, ".pdf")
  ggsave(filename, p2, width = 200, height = 120, units = "mm")
  
}






