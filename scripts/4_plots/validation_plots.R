library(RevGadgets)
library(grid)
library(ggplot2)
library(tidyverse)
library(tikzDevice)


bm1 <- readTrace(path = "output/1_validation/state_less_BM/state_less_BM_run_1.log", burnin = 0.05)
sdou1 <-  readTrace(path = "output/1_validation/state_less_BM/state_dependent_OU_run_1.log", burnin = 0.05)
sdbm1 <- readTrace(path = "output/1_validation/state_less_BM/state_dependent_BM_run_1.log", burnin = 0.05)
ou2 <- readTrace(path = "output/1_validation/state_less_BM/state_less_OU_run_1.log", burnin = 0.05)
bm1[[1]]$sigma2_sdou <- sdou1[[1]]$sigma2
bm1[[1]]$sigma2_ou <- ou2[[1]]$sigma2
bm1[[1]]$sigma2_sdbm <- sdbm1[[1]]$sigma2


colors_bm <- c("sigma2" = "#ABC3C9", "sigma2_sdbm" = "#CCBe9F",
               "sigma2_ou"="#FFCCCC", "sigma2_sdou"="#663333")
plot_bm <- plotTrace(trace = bm1, vars = c("sigma2", "sigma2_sdbm",
                                          "sigma2_ou", "sigma2_sdou"),
                     color=colors_bm)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("s")

pdf("figures/stateless_bm.pdf")
plot_bm
dev.off()




ou <- readTrace(path = "output/1_validation/state_less_OU/state_less_OU_run_2.log", burnin = 0.05)
#ou[[1]] <- ou[[1]] %>% 
#  mutate(stv = sigma2 / (2 * alpha))
sdou2 <-  readTrace(path = "output/1_validation/state_less_OU/state_dependent_OU_run_2.log", burnin = 0.05)
#sdou2[[1]] <- sdou2[[1]] %>% 
#  mutate(stv = sigma2 / (2 * alpha))


ou[[1]]$sigma2_sdou <- sdou2[[1]]$sigma2
ou[[1]]$theta_sdou <- sdou2[[1]]$theta
ou[[1]]$alpha_sdou <- sdou2[[1]]$alpha
ou[[1]]$stv_sdou <- sdou2[[1]]$`stVs[1]`

colors_ou <- c("stv_sdou" = "#663333", "stV" = "#FFCCCC")
ou_stv <- plotTrace(trace = ou, vars = c("stV", "stv_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("v")
colors_ou <- c("sigma2_sdou" = "#663333", "sigma2" = "#FFCCCC")
ou_sigma2 <- plotTrace(trace = ou, vars = c("sigma2", "sigma2_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("s")
colors_ou <- c("theta_sdou" = "#663333", "theta" = "#FFCCCC")
ou_theta <- plotTrace(trace = ou, vars = c("theta", "theta_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("t")
colors_ou <- c("alpha_sdou" = "#663333", "alpha" = "#FFCCCC")
ou_alpha <- plotTrace(trace = ou, vars = c("alpha", "alpha_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("a")


pdf("figures/stateless_ou.pdf")
ou_stv
ou_alpha
ou_sigma2
ou_theta
dev.off()

sdbm <- readTrace(path = "output/1_validation/state_dependent_BM/state_dependent_BM_run_1.log", burnin = 0.05)
sdbm2020 <- readTrace(path = "output/1_validation/state_dependent_BM/state_dependent_BM_MayMoore_run_1.log", burnin = 0.05)
sdou <-  readTrace(path = "output/1_validation/state_dependent_BM/state_dependent_OU_run_1.log", burnin = 0.05)
sdbm2020[[1]]$`sigma2_sdou_0` <- sdou[[1]]$`sigma2s[1]`
sdbm2020[[1]]$`sigma2_sdou_1` <- sdou[[1]]$`sigma2s[2]`
sdbm2020[[1]]$`sigma2_0` <- sdbm[[1]]$`sigma2s[1]`
sdbm2020[[1]]$`sigma2_1` <- sdbm[[1]]$`sigma2s[2]`

colors_sdbm <- c("sigma2_sdou_0" = "#663333", "sigma2s[1]" = "#CCBe9F", "sigma2_0" = "#ABC3C9",
                 "sigma2_sdou_1" = "#663333", "sigma2s[2]" = "#CCBe9F", "sigma2_1" = "#ABC3C9")

plot_sdbm0 <- plotTrace(trace = sdbm2020,
          color = colors_sdbm,
          vars = c("sigma2_0", "sigma2_1", "sigma2s[1]", "sigma2s[2]", "sigma2_sdou_0", "sigma2_sdou_1"))[[1]] +
  ggtitle("s") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
  )

pdf("figures/state-dependent_bm.pdf")
plot_sdbm0
dev.off()






load("output/likelihood_comparison_across_methods.Rda")
logL <- logL %>% 
  filter(rb > -1500, rb < -800)
#logL <- logL[sample(1:length(logL$rb), 150, replace=FALSE),]
logL = logL[order(logL$rb),]
logL$Replicates = 1:length(logL$rb)

lnl_comparison <- ggplot(logL) +
  geom_point(shape=15, alpha=0.7, size=3, aes(x=Replicates, y=R_vcv),col="grey") +
  geom_point(shape=4,alpha=1.0,size=3, aes(x=Replicates, y=R_pruning),col="#0077bb") +
  geom_point(shape=3,alpha=1.0,size=2, aes(x=Replicates, y=rb),col="#ee7733") +
  ylab("Likelihood") +
  xlab("Replicate number") +
  theme_bw()



pdf("figures/likelihoods.pdf")
lnl_comparison
dev.off()


processValidation <- function(analysis_name, n_reps = 1000, n_bins = 50) {
  if (missing(analysis_name)) {
    stop("Missing required argument: Analysis Name")
  }
  
  # Define directories
  results_dir = paste0("results_", analysis_name)
  output_dir = paste0("output_", analysis_name)
  
  # create directories
  dir.create(results_dir, showWarnings = FALSE)
  
  # function to read data from the file 
  read_data <- function(file_path) {
    if (file.exists(file_path)) {
      return(read.table(file_path, 
                        sep = "\t", 
                        header = TRUE, 
                        skip = 0, 
                        check.names = FALSE))
    } else {
      return(NULL)
    }
  }
  
  # get the list of parameters 
  parameters <- colnames(read_data(paste0(output_dir, 
                                          "/Validation_Sim_0/posterior_samples.var")))
  parameters <- parameters[-1]   # remove the first column ("Iteration")
  
  # exclude "branch_rates" if present 
  parameters <- parameters[parameters != "branch_rates"]
  
  # Print list of parameters 
  cat("parameters:\n", parameters, "\n\n")
  
  # Iterate over each parameter
  for (param in parameters) {
    # initialize variables 
    coverage_probs <- data.frame(total_count = numeric(0), 
                                 in_count = numeric(0), 
                                 hpd_width = numeric(0), 
                                 stringsAsFactors = FALSE)
    hpd_width <- seq(from = 0.0, to = 1.0, by = 1/n_bins)
    
    for (i in 1:(n_bins+1)) {
      coverage_probs[i,] = c(total_count = 0, in_count = 0, hpd_width = hpd_width[i])
    }
    
    # initialize progress bar 
    pb <- txtProgressBar(min = 0, max = n_reps, char = "*", style = 3)
    
    # iterate over each replication 
    for (i in 1:n_reps) {
      setTxtProgressBar(pb, i)
      
      # read in the data 
      data <- read_data(paste0(output_dir, 
                               "/Validation_Sim_", 
                               i-1, 
                               "/posterior_samples.var"))
      if (is.null(data)) next
      
      # extract samples 
      num_samples = length(data[,1])
      
      x <- as.mcmc(data[round(0.25*num_samples):num_samples, param])
      
      true_val_ext <- ifelse(param == "branch_rates", ".out", ".txt")
      true_val <- read.table(file=paste0(output_dir, 
                                         "/Validation_Sim_", 
                                         i-1, 
                                         "/", 
                                         param, 
                                         true_val_ext))[1,1]
      
      # calculate coverage probabilities 
      for (k in 1:(n_bins + 1)) {
        hpd <- HPDinterval(x, prob = hpd_width[k])
        if (true_val >= hpd[1,1] && true_val <= hpd[1,2]){
          coverage_probs[k, "in_count"] <- coverage_probs[k, "in_count"] + 1
        }
        coverage_probs[k, "total_count"] <- coverage_probs[k, "total_count"] + 1
      }      
    }
    close(pb)
    
    # calculate frequency of coverage
    coverage_probs$freq = coverage_probs$in_count / coverage_probs$total_count
    
    # save coverage probabilities 
    saveRDS(coverage_probs, file = paste0(results_dir, "/", param, ".rds"))
    
    # print results to the screen 
    cat(param,"\n")
    cat("HPD-width:\t\t",hpd_width,"\n")
    cat("Coverage-freq:\t\t",coverage_probs$freq,"\n")
  }
  
  return(coverage_probs)
}


generate_coverage_plots <- function(analysis_name, 
                                    num_reps = 10000, 
                                    num_bins = 50, 
                                    results_dir, 
                                    figs_dir) {
  library(coda)
  library(ggplot2)
  
  # Create directories if they don't exist
  dir.create(figs_dir, showWarnings = FALSE)
  dir.create(results_dir, showWarnings = FALSE)
  
  # Read data
  in_file <- paste0("output/1_validation/", analysis_name, "/Validation_Sim_1/posterior_samples.var")
  data <- read.table(in_file, sep="\t", header=TRUE, skip=0, check.names=FALSE)
  parameters <- colnames(data)[-1]
  
  # Iterate over each parameter
  for (param in parameters) {
    # Read coverage probabilities
    coverage_probs <- readRDS(file = paste0(results_dir, "/", param, ".rds"))
    
    # Generate plot
    p <- ggplot(coverage_probs) +
      geom_bar(stat="identity", aes(x=hpd_width, y=freq), colour="lightgray", fill="lightgray") +
      theme_classic() +
      xlab("HPD width") + ylab("coverage probability") + ggtitle(param) +
      geom_segment(aes(x=0, y=0, xend=1, yend=1), linetype="dashed", size=1.5, show.legend=FALSE) +
      theme(legend.position="none", plot.title = element_text(hjust = 0.5))
    
    # Save plot
    ggsave(paste0(figs_dir, "/hpd_width_vs_coverage_", param, ".pdf"), plot=p, width=10, height=10, units="cm")
  }
}