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
}s