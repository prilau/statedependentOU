library(RevGadgets)
library(tidyverse)
library(grid)

process_par <- function(par, sdORsx, num_sim){
  df_par <- pivot_longer(dfx, starts_with(par), names_to = "par")
  df_med_true <- df_par %>%
    group_by(par, run) %>%
    summarize(med = median(value))
  par_0 <- paste0(par, "_0")
  par_1 <- paste0(par, "_1")
  if (sdORsx == "sd"){
    df_med_true$true <- c(rep(pars_sd[[par_0]][num_sim], 2), rep(pars_sd[[par_1]][num_sim], 2))
  } else {
    df_med_true$true <- c(rep(pars_stateless[[par_0]][num_sim], 4))
  }
  return(list(df_par, df_med_true))
}

plot_par <- function(par_name, df, med_true, dir_out, combo){
  p <- ggplot(df, aes(x = value, fill = par)) +
    geom_density(aes(x = value, after_stat(density), fill=run), color = "black", alpha = 0.4) +
    theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.7, 0.8)) +
    facet_grid(cols = vars(par), scales = "free_y") +
    geom_vline(mapping = aes(xintercept = med), data = med_true, linetype = "dashed", color = "darkgrey") + 
    geom_vline(mapping = aes(xintercept = true), data = med_true, linetype = "dashed", color = "brown") + 
    #scale_x_log10(breaks = c(p_breaks[[par_name]])) +
    #coord_cartesian(xlim = p_xlim[[par_name]]) +
    labs(x = combo, y = "posterior density")
  
  filename <- paste0(dir_out, combo, ".pdf")
  ggsave(filename, p, width = 200, height = 120, units = "mm")
}



load("data/2_simulation/convergence/pars_sd.Rda")
load("data/2_simulation/convergence/pars_stateless.Rda")
pars <- c("theta", "alpha", "sigma2", "rho", "stv", "halflife")
dir_in = "output/2_simulation/convergence/logs/"
dir_out = "figures/2_simulation/convergence/"

# plot non-DPP
num_sim = 3
sdORsx   = c("sd", "stateless")
moves       = c("mv", "mvAVM", "2mv")
link = c("linkedPrior", "unlinkedPrior")

grid = expand.grid(num_sim=1:num_sim, sdORsx=sdORsx, moves=moves, link=link, stringsAsFactors=FALSE)

for (i in 1:nrow(grid)){
  this_row = grid[i,]
  this_num_sim = this_row[[1]]
  this_state = this_row[[2]]
  this_move = this_row[[3]]
  this_link = this_row[[4]]
  this_combo = paste0(c("sim", this_num_sim, this_state, this_move, this_link), collapse = "_")
  read_run1 <- paste0(dir_in, paste0(c("sim", this_num_sim, this_state, "run_1", this_move, this_link), collapse = "_"), ".log")
  read_run2 <- paste0(dir_in, paste0(c("sim", this_num_sim, this_state, "run_2", this_move, this_link), collapse = "_"), ".log")
  
  if ( isFALSE(file.exists(read_run1)) | isFALSE(file.exists(read_run2))) {
    next
  }
  
  df <- readTrace(path = read_run1, burnin = 0.1)
  df[[2]] <- readTrace(path = read_run2, burnin = 0.1)[[1]]
  
  df[[1]]$run <- "run 1"
  df[[2]]$run <- "run 2"
  
  dfx <- bind_rows(df[[1]], df[[2]]) %>%
    select(-starts_with(c("theta_compare", "alpha_compare", "rho_compare", "sigma2_compare", "halflife_compare", "stv_compare")))
  
  for (par in pars){
    if (this_state == "sd"){
      dfs <- process_par(par, "sd", this_num_sim)
    } else {
      dfs <- process_par(par, "stateless", this_num_sim)
    }
    df_par <- dfs[[1]]
    df_par_medtrue <- dfs[[2]]
    this_combo_par <- paste0(this_combo, "_", par)
    plot_par(par_name=par, df_par, df_par_medtrue, dir_out, this_combo_par)
  }
}



# plot DPP as well
num_sim = 3
sdORsx   = c("sd", "stateless")

grid = expand.grid(num_sim=1:num_sim, sdORsx=sdORsx, stringsAsFactors=FALSE)

for (i in 1:nrow(grid)){
  this_row = grid[i,]
  this_num_sim = this_row[[1]]
  this_state = this_row[[2]]
  this_combo = paste0(c("sim", this_num_sim, this_state, "DPP"), collapse = "_")
  read_run1 <- paste0(dir_in, paste0(c("sim", this_num_sim, this_state, "run_1_DPP.log"), collapse = "_"))
  read_run2 <- paste0(dir_in, paste0(c("sim", this_num_sim, this_state, "run_2_DPP.log"), collapse = "_"))
  
  if ( isFALSE(file.exists(read_run1)) | isFALSE(file.exists(read_run2))) {
    next
  }
  
  df <- readTrace(path = read_run1, burnin = 0.1)
  df[[2]] <- readTrace(path = read_run2, burnin = 0.1)[[1]]
  
  df[[1]]$run <- "run 1"
  df[[2]]$run <- "run 2"
  
  dfx <- bind_rows(df[[1]], df[[2]]) %>%
    select(-starts_with(c("theta_compare", "alpha_compare", "rho_compare", "sigma2_compare", "halflife_compare", "stv_compare")))
  
  for (par in pars){
    if (this_state == "sd"){
      dfs <- process_par(par, "sd", this_num_sim)
    } else {
      dfs <- process_par(par, "stateless", this_num_sim)
    }
    df_par <- dfs[[1]]
    df_par_medtrue <- dfs[[2]]
    this_combo_par <- paste0(this_combo, "_", par)
    plot_par(par_name=par, df_par, df_par_medtrue, dir_out, this_combo_par)
  }
}



# Plot set 2
num_sim = 3
sdORsx   = c("sd", "stateless")
link = c("linkedPrior", "unlinkedPrior")

grid = expand.grid(num_sim=1:num_sim, sdORsx=sdORsx, link=link, stringsAsFactors=FALSE)

for (i in 1:nrow(grid)){
  this_row = grid[i,]
  this_num_sim = this_row[[1]]
  this_state = this_row[[2]]
  this_link = this_row[[3]]
  this_combo = paste0(c("sim", this_num_sim, this_state, this_link, "set_2"), collapse = "_")
  read_run1 <- paste0(dir_in, paste0(c("sim", this_num_sim, this_state, "run_1_2mv", this_link, "set_2"), collapse = "_"), ".log")
  read_run2 <- paste0(dir_in, paste0(c("sim", this_num_sim, this_state, "run_2_2mv", this_link, "set_2"), collapse = "_"), ".log")
  
  if ( isFALSE(file.exists(read_run1)) | isFALSE(file.exists(read_run2))) {
    next
  }
  
  df <- readTrace(path = read_run1, burnin = 0.1)
  df[[2]] <- readTrace(path = read_run2, burnin = 0.1)[[1]]
  
  df[[1]]$run <- "run 1"
  df[[2]]$run <- "run 2"
  
  dfx <- bind_rows(df[[1]], df[[2]]) %>%
    select(-starts_with(c("theta_compare", "alpha_compare", "rho_compare", "sigma2_compare", "halflife_compare", "stv_compare")))
  
  for (par in pars){
    if (this_state == "sd"){
      dfs <- process_par(par, "sd", this_num_sim)
    } else {
      dfs <- process_par(par, "stateless", this_num_sim)
    }
    df_par <- dfs[[1]]
    df_par_medtrue <- dfs[[2]]
    this_combo_par <- paste0(this_combo, "_", par)
    plot_par(par_name=par, df_par, df_par_medtrue, dir_out, this_combo_par)
  }
}


