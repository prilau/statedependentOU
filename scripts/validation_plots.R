library(RevGadgets)
library(grid)
library(ggplot2)
library(tidyverse)
library(tikzDevice)


bm1 <- readTrace(path = "output/state_less_BM/state_less_BM_run_1.log", burnin = 0.05)
sdou1 <-  readTrace(path = "output/state_less_BM/state_dependent_OU_run_1.log", burnin = 0.05)
sdbm1 <- readTrace(path = "output/state_less_BM/state_dependent_BM_run_1.log", burnin = 0.05)
ou2 <- readTrace(path = "output/state_less_BM/state_less_OU_run_1.log", burnin = 0.05)
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




ou <- readTrace(path = "output/state_less_OU/state_less_OU_run_2.log", burnin = 0.05)
ou[[1]] <- ou[[1]] %>% 
  mutate(stv = sigma2 / (2 * alpha))
sdou2 <-  readTrace(path = "output/state_less_OU/state_dependent_OU_run_2.log", burnin = 0.05)
sdou2[[1]] <- sdou2[[1]] %>% 
  mutate(stv = sigma2 / (2 * alpha))


ou[[1]]$sigma2_sdou <- sdou2[[1]]$sigma2
ou[[1]]$theta_sdou <- sdou2[[1]]$theta
ou[[1]]$alpha_sdou <- sdou2[[1]]$alpha
ou[[1]]$stv_sdou <- sdou2[[1]]$stv

colors_ou <- c("stv_sdou" = "#663333", "stv" = "#FFCCCC")
ou_stv <- plotTrace(trace = ou, vars = c("stv", "stv_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("")
colors_ou <- c("sigma2_sdou" = "#663333", "sigma2" = "#FFCCCC")
ou_sigma2 <- plotTrace(trace = ou, vars = c("sigma2", "sigma2_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("")
colors_ou <- c("theta_sdou" = "#663333", "theta" = "#FFCCCC")
ou_theta <- plotTrace(trace = ou, vars = c("theta", "theta_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("")
colors_ou <- c("alpha_sdou" = "#663333", "alpha" = "#FFCCCC")
ou_alpha <- plotTrace(trace = ou, vars = c("alpha", "alpha_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  ggtitle("")


pdf("figures/stateless_ou.pdf")
ou_stv
ou_alpha
ou_sigma2
ou_theta
dev.off()





sdbm <- readTrace(path = "output/state_dependent_BM/state_dependent_BM_run_1.log", burnin = 0.05)
sdbm2020 <- readTrace(path = "output/state_dependent_BM/state_dependent_BM_MayMoore_run_1.log", burnin = 0.05)
sdou <-  readTrace(path = "output/state_dependent_BM/state_dependent_OU_run_1.log", burnin = 0.05)
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

