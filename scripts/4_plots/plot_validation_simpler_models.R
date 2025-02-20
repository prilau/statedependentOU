library(RevGadgets)
library(grid)
library(ggplot2)
library(tidyverse)
library(latex2exp)
source("scripts/5_miscellaneous/functions.R")

bm1 <- readTrace(path = "output/1_validation/state_less_BM/state_less_BM_run_1.log", burnin = 0.05)
sdou1 <-  readTrace(path = "output/1_validation/state_less_BM/state_dependent_OU_run_1.log", burnin = 0.05)
sdbm1 <- readTrace(path = "output/1_validation/state_less_BM/state_dependent_BM_run_1.log", burnin = 0.05)
ou2 <- readTrace(path = "output/1_validation/state_less_BM/state_less_OU_run_1.log", burnin = 0.05)
bm1[[1]]$sigma2_sdou <- sdou1[[1]]$sigma2
bm1[[1]]$sigma2_ou <- ou2[[1]]$sigma2
bm1[[1]]$sigma2_sdbm <- sdbm1[[1]]$sigma2

colors_all <- c('#66CCEE', '#EE6677', '#228833', '#4477AA', '#CCBB44') 
names(colors_all) <- c("bm", "sdbm", "ou", "sdou", "musscrat")

colors_bm <- colors_all[1:4]
names(colors_bm) <- paste0("sigma2_", names(colors_bm))
plot_bm <- plotTrace(trace = bm1, vars = c("sigma2", "sigma2_sdbm",
                                          "sigma2_ou", "sigma2_sdou"),
                     color=colors_bm)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank()) +
  ggtitle(TeX("$\\sigma^2$"))
  



ou <- readTrace(path = "output/1_validation/state_less_OU/state_less_OU_run_2.log", burnin = 0.05)
#ou[[1]] <- ou[[1]] %>% 
#  mutate(stv = sigma2 / (2 * alpha))
sdou2 <-  readTrace(path = "output/1_validation/state_less_OU/state_dependent_OU_run_2.log", burnin = 0.05)
#sdou2[[1]] <- sdou2[[1]] %>% 
#  mutate(stv = sigma2 / (2 * alpha))


ou[[1]]$sigma2_sdou <- sdou2[[1]]$sigma2
ou[[1]]$theta_sdou <- sdou2[[1]]$theta
ou[[1]]$alpha_sdou <- sdou2[[1]]$alpha
ou[[1]]$halflife_sdou <- sdou2[[1]]$halflife
ou[[1]]$stv_sdou <- sdou2[[1]]$`stVs[1]`

colors_ou <- colors_all[which(names(colors_all) %in% c("ou", "sdou"))]
names(colors_ou) <- c("halflife", "halflife_sdou")
ou_halflife <- plotTrace(trace = ou, vars = c("halflife", "halflife_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank()) +
  ggtitle(TeX("$t_{0.5}$"))

names(colors_ou) <- c("stV", "stv_sdou")
ou_stv <- plotTrace(trace = ou, vars = c("stV", "stv_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank()) +
  ggtitle(TeX("$V_y$"))

#names(colors_ou) <- c("sigma2_sdou" = "#663333", "sigma2" = "#FFCCCC")
#ou_sigma2 <- plotTrace(trace = ou, vars = c("sigma2", "sigma2_sdou"), color=colors_ou)[[1]] +
#  theme(panel.grid = element_blank(),
#        legend.position = "none") +
#  ggtitle("s")
names(colors_ou) <- c("theta", "theta_sdou")
ou_theta <- plotTrace(trace = ou, vars = c("theta", "theta_sdou"), color=colors_ou)[[1]] +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  ggtitle(TeX("$\\theta$"))

#names(colors_ou) <- c("alpha_sdou" = "#663333", "alpha" = "#FFCCCC")
#ou_alpha <- plotTrace(trace = ou, vars = c("alpha", "alpha_sdou"), color=colors_ou)[[1]] +
#  theme(panel.grid = element_blank(),
#        legend.position = "none") +
#  ggtitle("a")



sdbm <- readTrace(path = "output/1_validation/state_dependent_BM/state_dependent_BM_run_1.log", burnin = 0.05)
sdbm2020 <- readTrace(path = "output/1_validation/state_dependent_BM/state_dependent_BM_MayMoore_run_1.log", burnin = 0.05)
sdou <-  readTrace(path = "output/1_validation/state_dependent_BM/state_dependent_OU_run_1.log", burnin = 0.05)
sdbm2020[[1]]$`sigma2_sdou_0` <- sdou[[1]]$`sigma2s[1]`
sdbm2020[[1]]$`sigma2_sdou_1` <- sdou[[1]]$`sigma2s[2]`
sdbm2020[[1]]$`sigma2_0` <- sdbm[[1]]$`sigma2s[1]`
sdbm2020[[1]]$`sigma2_1` <- sdbm[[1]]$`sigma2s[2]`

colors_sdbm <- colors_all[which(names(colors_all) %in% c("sdbm", "sdou", "musscrat"))]
names(colors_sdbm) <- c("sigma2_0", "sigma2_sdou_0", "sigma2s[1]")
plot_sdbm0 <- plotTrace(trace = sdbm2020,
          color = colors_sdbm,
          vars = c("sigma2_0", "sigma2s[1]", "sigma2_sdou_0"))[[1]] +
  ggtitle(TeX("$\\sigma^2_0$")) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )

names(colors_sdbm) <- c("sigma2_1", "sigma2_sdou_1", "sigma2s[2]")
plot_sdbm1 <- plotTrace(trace = sdbm2020,
                        color = colors_sdbm,
                        vars = c("sigma2_1", "sigma2s[2]", "sigma2_sdou_1"))[[1]] +
  ggtitle(TeX("$\\sigma^2_1$")) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank()
  )


# dummy plot with all color variables
bm1[[1]]$bm <- bm1[[1]]$sigma2
bm1[[1]]$sdbm <- bm1[[1]]$sigma2_sdou
bm1[[1]]$ou <- bm1[[1]]$sigma2_ou
bm1[[1]]$sdou <- bm1[[1]]$sigma2_sdbm
bm1[[1]]$musscrat <- sdbm2020[[1]]$`sigma2s[1]`


colors_all <- c('#66CCEE', '#EE6677', '#228833', '#4477AA', '#CCBB44') 
names(colors_all) <- c("bm", "sdbm", "ou", "sdou", "musscrat")

plot_dummy <- plotTrace(trace = bm1,
                        color = colors_all,
                        vars = c("bm", "sdbm", "ou", "sdou", "musscrat"))[[1]]

legend <- get_legend2(plot_dummy + theme(legend.position = "left",
                       legend.box.margin = margin(0, 0, 0, 12))
                      + scale_color_manual(values=c('#66CCEE', '#EE6677', '#228833', '#4477AA', '#CCBB44'), 
                                           name="Model",
                                           labels=c("State-independent BM", "State-dependent BM",
                                                    "State-independent OU", "State-dependent OU", "MuSSCRat"))
                      + guides(size = "none",
                                color = guide_legend(override.aes = list(size = 3),
                                                     title='Model'),
                                fill=guide_legend(title='Model')))



val_col1row1 <- cowplot::plot_grid(plot_bm, plot_sdbm0, plot_sdbm1, ncol=3, labels=c("(a)", "(b)", "(c)"))
val_col1row2 <- cowplot::plot_grid(ou_halflife, ou_stv, ou_theta, ncol=3, labels=c("(d)", "(e)", "(f)"))
val_col2 <- cowplot::plot_grid(legend)
val_col1 <- cowplot::plot_grid(val_col1row1, val_col1row2, ncol=1)
val_all <- cowplot::plot_grid(val_col1, val_col2, ncol=2, rel_widths = c(0.75, 0.25))

val_all

ggsave("figures/1_validation/compare_simple_models.pdf", val_all, width=8, height=6, units="in")