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
load("data/2_simulation/binState/pars_sd.Rda")

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
  df_med_true$true <- c(pars_sd$theta_0[i] - pars_sd$theta_1[i], 0)
  
  if (!dir.exists(dir_out)){
    dir.create(dir_out, showWarnings = FALSE)
  }
  plot_dthetas(dfx, df_med_true, dir_out, num_sim = i)
}



# plot true dthetas against P(state_0 > state_1)
num_sim=5
dir_in="output/2_simulation/binState/logs/"
dir_out="figures/2_simulation/binState/"
load("data/2_simulation/binState/pars_sd.Rda")
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
