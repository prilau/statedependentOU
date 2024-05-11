library(RevGadgets)
library(grid)
library(ggplot2)
library(tidyverse)
library(tikzDevice)


#trace_quant_bm1 <- readTrace(path = "state_less_BM/state_less_BM_run_1.log", burnin = 0.05)
#trace_quant_sdou1 <-  readTrace(path = "state_less_BM/state_dependent_OU_run_1.log", burnin = 0.05)
#trace_quant_bm1[[1]]$sigma2_sdou <- trace_quant_sdou1[[1]]$sigma2
#
#plotTrace(trace = trace_quant_bm1, vars = c("sigma2", "sigma2_sdou"))[[1]]

#plots_ou <- list()

ou <- readTrace(path = "state_less_OU/state_less_OU_run_2.log", burnin = 0.05)
ou[[1]] <- ou[[1]] %>% 
  mutate(stv = sigma2 / (2 * alpha))
sdou2 <-  readTrace(path = "state_less_OU/state_dependent_OU_run_2.log", burnin = 0.05)
sdou2[[1]] <- sdou2[[1]] %>% 
  mutate(stv = sigma2 / (2 * alpha))


#ou[[1]]$sigma2_sdou <- sdou2[[1]]$sigma2
#ou[[1]]$theta_sdou <- sdou2[[1]]$theta
#ou[[1]]$alpha_sdou <- sdou2[[1]]$alpha
ou[[1]]$stv_sdou <- sdou2[[1]]$stv

#plots_ou[[1]] <- plotTrace(trace = ou, vars = c("sigma2", "sigma2_sdou"))[[1]] +
#  theme(panel.grid = element_blank(),
#        legend.position = "none") +
#  ggtitle("")
#plots_ou[[2]] <- plotTrace(trace = ou, vars = c("theta", "theta_sdou"))[[1]] +
#  theme(panel.grid = element_blank(),
#        legend.position = "none") +
#  ggtitle("")
#plots_ou[[3]] <- plotTrace(trace = ou, vars = c("alpha", "alpha_sdou"))[[1]] +
#  theme(panel.grid = element_blank(),
#        legend.position = "none") +
#  ggtitle("")

colors_ou <- c("stv_sdou" = "#663333", "stv" = "#FFCCCC")
plots_ou <- plotTrace(trace = ou, color = colors_ou, vars = c("stv", "stv_sdou"))[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("")


pdf("ou_stv.pdf", width=5, height=6)
plots_ou
dev.off()

#pdf("testing4.pdf", width=5, height=6)
#plots_ou[[2]]
#dev.off()
#
#pdf("testing5.pdf", width=5, height=6)
#plots_ou[[3]]
#dev.off()


















#sdbm <- readTrace(path = "state_dependent_BM/state_dependent_BM_run_1.log", burnin = 0.05)
sdbm2020 <- readTrace(path = "state_dependent_BM/state_dependent_BM_MayMoore_run_1.log", burnin = 0.05)
sdou <-  readTrace(path = "state_dependent_BM/state_dependent_OU_run_1.log", burnin = 0.05)
sdbm2020[[1]]$`sigma2_sdou_0` <- sdou[[1]]$`sigma2s[1]`
sdbm2020[[1]]$`sigma2_sdou_1` <- sdou[[1]]$`sigma2s[2]`
#sdbm[[1]]$`sigma2_sdbm2020_0` <- sdbm2020[[1]]$`sigma2s[1]`
#sdbm[[1]]$`sigma2_sdbm2020_1` <- sdbm2020[[1]]$`sigma2s[2]`

#sdbm_state0 <- list()
#sdbm_state1 <- list()
#sdbm_state0[[1]] <- sdbm[[1]] %>% 
#  select(-c(`sigma2s_sdou[2]`, `sigma2s[2]`, `sigma2s_sdbm2020[2]`)) %>% 
#  rename(`state-dependent OU`=`sigma2s_sdou[1]`,
#         `state-dependent BM`=`sigma2s[1]`,
#         `state-dependent BM (May & Moore 2020)`=`sigma2s_sdbm2020[1]`)
#sdbm_state1[[1]] <- sdbm[[1]] %>% 
#  select(-c(`sigma2s_sdou[1]`, `sigma2s[1]`, `sigma2s_sdbm2020[1]`)) %>% 
#  rename(`state-dependent OU`=`sigma2s_sdou[2]`,
#         `state-dependent BM`=`sigma2s[2]`,
#         `state-dependent BM (May & Moore 2020)`=`sigma2s_sdbm2020[2]`)
#004488', '#DDAA33', '#BB5566'.

colors_sdbm0 <- c("sigma2_sdou_0" = "#666633", "sigma2s[1]" = "#EEEEBB")
colors_sdbm1 <- c("sigma2_sdou_1" = "#ee7733", "sigma2s[2]" = "#ffccaa")

plots_sdbm0 <- plotTrace(trace = sdbm2020,
          color = colors_sdbm0,
          vars = c("sigma2s[1]", "sigma2_sdou_0"))[[1]] +
  ggtitle("") +
  theme(panel.grid = element_blank(),
        #axis.text = element_text(size = 16),
        #legend.text = element_text(size = 16),
        #legend.title = element_text(size = 16),
        legend.position = "none",
        #legend.position.inside = c(0.9, 0.9),
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #legend.margin = margin(4, 4, 4, 4),
        axis.title.x = element_text(size = 10),
        #axis.text.y = element_text(size = 8)
        axis.title.y = element_text(size = 10)
  ) +
  xlim(0,0.10)
  
plots_sdbm1 <- plotTrace(trace = sdbm2020,
                         color = colors_sdbm1,
                         vars = c("sigma2s[2]", "sigma2_sdou_1"))[[1]] +
  ggtitle("") +
  theme(panel.grid = element_blank(),
        #axis.text = element_text(size = 16),
        #legend.text = element_text(size = 16),
        #legend.title = element_text(size = 16),
        legend.position = "none",
        #legend.position.inside = c(0.9, 0.9),
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #legend.margin = margin(4, 4, 4, 4),
        axis.title.x = element_text(size = 10),
        #axis.text.y = element_text(size = 8)
        axis.title.y = element_text(size = 10)
  )

pdf("sdbm_sigma0.pdf", width=5, height=6)
plots_sdbm0
dev.off()

pdf("testing1.pdf", width=5, height=6)
plots_sdbm1
dev.off()





load("output/likelihood_comparison_across_methods.Rda")
logL <- logL %>% 
  filter(rb > -1500, rb < -800)
logL <- logL[sample(1:length(logL$rb), 150, replace=FALSE),]
logL = logL[order(logL$rb),]
logL$Replicates = 1:length(logL$rb)

lnl_comparison <- ggplot(logL) +
  geom_point(shape=15, alpha=0.7, aes(x=Replicates, y=R_vcv),col="black") +
  geom_point(shape=3,alpha=0.7,aes(x=Replicates, y=R_pruning),col="#ee7733") +
  geom_point(shape=4,alpha=0.7,aes(x=Replicates, y=rb),col="#0077bb") +
  ylab("Likelihood") +
  xlab("Replicate number") +
  theme_bw()



pdf("figures/likelihoods.pdf")
lnl_comparison
dev.off()
