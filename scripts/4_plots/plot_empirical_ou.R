library(cowplot)
library(RevGadgets)
library(ggplot2)
library(tidyverse)
library(latex2exp)
library(grid)
library(gridExtra)
source("scripts/5_miscellaneous/functions.R")

tree <- read.tree("data/3_empirical/mammal_perMY_r500.tre")
trait <- read.csv("data/3_empirical/mammal_traits.csv") %>% filter(Binomial.1.2 %in% tree$tip.label) %>% 
  mutate(diet3 = as.factor(diet3),
         diet4 = as.factor(diet4))
root_age <- max(node.depth.edgelength(tree))

###########
# 3-state #
###########
trace_sdou3 <- readTrace(c("output/3_empirical/sdOU_r500_3StateOrderedModel/trace_run_1.log", "output/3_empirical/sdOU_r500_3StateOrderedModel/trace_run_4.log"), burnin = 0.1)
trace_sdou3[[1]]$`halflife[1]` <- trace_sdou3[[1]]$`halflife[1]` / root_age
trace_sdou3[[1]]$`halflife[2]` <- trace_sdou3[[1]]$`halflife[2]` / root_age
trace_sdou3[[1]]$`halflife[3]` <- trace_sdou3[[1]]$`halflife[3]` / root_age
trace_sdou3[[2]]$`halflife[1]` <- trace_sdou3[[2]]$`halflife[1]` / root_age
trace_sdou3[[2]]$`halflife[2]` <- trace_sdou3[[2]]$`halflife[2]` / root_age
trace_sdou3[[2]]$`halflife[3]` <- trace_sdou3[[2]]$`halflife[3]` / root_age

trait_tmp <- trait %>%
  mutate(diet3 = ifelse(diet3==0, "theta[1]", ifelse(diet3==1, "theta[2]", "theta[3]")))

color_sdou3 <- c("#44aa99", "#ddcc77", "#882255")

names(color_sdou3) <- c("halflife[1]", "halflife[2]", "halflife[3]")
p1 <- plotTrace(trace_sdou3, vars = c("halflife[1]", "halflife[2]", "halflife[3]"), color = color_sdou3)[[1]] +
  ggtitle(TeX("Phylogenetic half-life $t_{0.5}$")) +
  theme(legend.position = "none") +
  xlab("Time (tree height)") +
  ylab("Posterior probability density") +
  xlim(0,10)

names(color_sdou3) <- c("stv[1]", "stv[2]", "stv[3]")
p2 <- plotTrace(trace_sdou3, vars = c("stv[1]", "stv[2]", "stv[3]"), color = color_sdou3)[[1]] +
  ggtitle(TeX("Stationary variance $V_y$")) +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  xlab("ln(body mass (kg))^2") +
  xlim(0, 100)

names(color_sdou3) <- c("theta[1]", "theta[2]", "theta[3]")
p3 <- plotTrace(trace_sdou3, vars = c("theta[1]", "theta[2]", "theta[3]"), color = color_sdou3)[[1]] +
  ggtitle(TeX("Optimum $\\theta$")) +
  theme(axis.title.y = element_blank()) +
  scale_color_manual(values=c("#44aa99", "#ddcc77", "#882255"), 
                    name="Diet",
                    labels=c("Herbivore", "Omnivore", "Carnivore")) +
  xlab("ln(body mass (kg))") +
  annotate("rect", xmin = -8, xmax = 12, ymin = 0.25, ymax = 0.32, fill="grey80", alpha=0.5, size=0.1) +
  annotate("text", x=-2.5, y=0.31, label= "Empirical distribution", size=2.5) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet3=="theta[1]"),],
               aes(x=log_mass_kg, y=0.29),
               color="#44aa99",
               width=0.01,
               linewidth=0.2) +
  geom_point(data=trait_tmp[which(trait_tmp$diet3=="theta[1]"),], aes(x=log_mass_kg, y=0.29),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#44aa99",
             position=position_jitter(height=0.0025)) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet3=="theta[2]"),],
               aes(x=log_mass_kg, y=0.275),
               color="#ddcc77",
               width=0.01,
               linewidth=0.2) +
  geom_point(data=trait_tmp[which(trait_tmp$diet3=="theta[2]"),], aes(x=log_mass_kg, y=0.275),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#ddcc77",
             position=position_jitter(height=0.0025)) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet3=="theta[3]"),],
               aes(x=log_mass_kg, y=0.26),
               color="#882255",
               width=0.01,
               linewidth=0.2,
               outlier.shape=NA) +
  geom_point(data=trait_tmp[which(trait_tmp$diet3=="theta[3]"),], aes(x=log_mass_kg, y=0.26),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#882255",
             position=position_jitter(height=0.0025)) +
  theme(legend.position = "none")
p3  


post_sdou3 <- cowplot::plot_grid(p1, p2, p3, ncol=3)
legend3 <- get_legend2(p3 + theme(legend.position = "right",
                                 legend.box.margin = margin(0, 0, 0, 12)))

post_sdou3_wlegend <- plot_grid(post_sdou3, legend3, rel_widths = c(3, .4))

ggsave("figures/3_empirical/sdOU_r500_3StateOrderedModel/ou_posterior.pdf", post_sdou3_wlegend, width = 10, height = 4, units = "in")



###########
# 4-state #
###########

trace_sdou4 <- readTrace("output/3_empirical/sdOU_r500_4StateModel/trace_run_1.log", burnin = 0.1)
trace_sdou4[[1]]$`halflife[1]` <- trace_sdou4[[1]]$`halflife[1]` / root_age
trace_sdou4[[1]]$`halflife[2]` <- trace_sdou4[[1]]$`halflife[2]` / root_age
trace_sdou4[[1]]$`halflife[3]` <- trace_sdou4[[1]]$`halflife[3]` / root_age
trace_sdou4[[1]]$`halflife[4]` <- trace_sdou4[[1]]$`halflife[4]` / root_age

trait_tmp <- trait %>%
  mutate(diet4 = ifelse(diet4==0, "theta[1]", ifelse(diet4==1, "theta[2]",ifelse(diet4==2, "theta[3]", "theta[4]"))))

color_sdou4 <- c("#44aa99", "#999933", "#CC6677", "#882255")

names(color_sdou4) <- c("halflife[1]", "halflife[2]", "halflife[3]", "halflife[4]")
p7 <- plotTrace(trace_sdou4, vars = c("halflife[1]", "halflife[2]", "halflife[3]", "halflife[4]"), color = color_sdou4)[[1]] +
  ggtitle(TeX("Phylogenetic half-life $t_{0.5}$")) +
  theme(legend.position = "none") +
  ylab("Posterior probability density") +
  xlim(0,15) +
  xlab("Time (tree height)")
names(color_sdou4) <- c("stv[1]", "stv[2]", "stv[3]", "stv[4]")
p8 <- plotTrace(trace_sdou4, vars = c("stv[1]", "stv[2]", "stv[3]", "stv[4]"), color = color_sdou4)[[1]] +
  ggtitle(TeX("Stationary variance $V_y$")) +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  xlim(0,100) +
  xlab("ln(body mass (kg))^2")
names(color_sdou4) <- c("theta[1]", "theta[2]", "theta[3]", "theta[4]")
p9 <- plotTrace(trace_sdou4, vars = c("theta[1]", "theta[2]", "theta[3]", "theta[4]"), color = color_sdou4)[[1]] +
  ggtitle(TeX("Optimum $\\theta$")) +
  theme(axis.title.y = element_blank(),
        legend.position = "none") +
  xlab("ln(body mass (kg))") +
  scale_color_manual(values=c("#44aa99", "#999933", "#CC6677", "#882255"), 
                     name="Diet",
                     labels=c("Herbivore", TeX("Omnivore ($>50\\%$ plants)"), TeX("Omnivore ($\\leq 50\\%$ plants)"), "Carnivore")) +
  annotate("rect", xmin = -8, xmax = 12, ymin = 0.25, ymax = 0.33, fill="grey80", alpha=0.5, size=0.1) +
  annotate("text", x=0, y=0.32, label= "Empirical distribution", size=2.5) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet4=="theta[1]"),],
               aes(x=log_mass_kg, y=0.305),
               color="#44aa99",
               width=0.01,
               linewidth=0.2) +
  geom_point(data=trait_tmp[which(trait_tmp$diet4=="theta[1]"),], aes(x=log_mass_kg, y=0.305),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#44aa99",
             position=position_jitter(height=0.0025)) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet4=="theta[2]"),],
               aes(x=log_mass_kg, y=0.29),
               color="#999933",
               width=0.01,
               linewidth=0.2,
               outlier.shape=NA) +
  geom_point(data=trait_tmp[which(trait_tmp$diet4=="theta[2]"),], aes(x=log_mass_kg, y=0.29),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#999933",
             position=position_jitter(height=0.0025)) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet4=="theta[3]"),],
               aes(x=log_mass_kg, y=0.275),
               color="#CC6677",
               width=0.01,
               linewidth=0.2) +
  geom_point(data=trait_tmp[which(trait_tmp$diet4=="theta[3]"),], aes(x=log_mass_kg, y=0.275),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#CC6677",
             position=position_jitter(height=0.0025)) +
  geom_boxplot(data=trait_tmp[which(trait_tmp$diet4=="theta[3]"),],
               aes(x=log_mass_kg, y=0.26),
               color="#882255",
               width=0.01,
               linewidth=0.2) +
  geom_point(data=trait_tmp[which(trait_tmp$diet4=="theta[3]"),], aes(x=log_mass_kg, y=0.26),
             shape = 16,
             size = 0.3,
             alpha = 0.3,
             color="#882255",
             position=position_jitter(height=0.0025))
p9

post_sdou4 <- cowplot::plot_grid(p7, p8, p9, ncol=3)
legend4 <- get_legend2(p9 + theme(legend.position = "right",
                                 legend.box.margin = margin(0, 0, 0, 12)))

post_sdou4_wlegend <- plot_grid(post_sdou4, legend4, rel_widths = c(3, .7))

ggsave("figures/3_empirical/sdOU_r500_4StateModel/ou_posterior.pdf", post_sdou4_wlegend, width = 11.5, height = 4, units = "in")


###########
# st dist #
###########
trace_tmp3 <- rbind(trace_sdou3[[1]], trace_sdou3[[2]]) %>% 
  summarize(theta1 = median(`theta[1]`),
            theta2 = median(`theta[2]`),
            theta3 = median(`theta[3]`),
            stv1 = median(`stv[1]`),
            stv2 = median(`stv[2]`),
            stv3 = median(`stv[3]`))

st_dist_3 <- 
  ggplot() +
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp3$theta1, sd = sqrt(trace_tmp3$stv1)), color="#44aa99") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp3$theta1, sd=sqrt(trace_tmp3$stv1)), 
                geom = "area", fill = "#44aa99", alpha = .5) +
  geom_vline(xintercept = trace_tmp3$theta1, color="#44aa99", linetype="dashed", linewidth=0.6) +
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp3$theta2, sd = sqrt(trace_tmp3$stv2)), color="#ddcc77") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp3$theta2, sd=sqrt(trace_tmp3$stv2)), 
                geom = "area", fill = "#ddcc77", alpha = .5) +
  geom_vline(xintercept = trace_tmp3$theta2, color="#ddcc77", linetype="dashed", linewidth=0.6) +
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp3$theta3, sd = sqrt(trace_tmp3$stv3)), color="#882255") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp3$theta3, sd=sqrt(trace_tmp3$stv3)), 
                geom = "area", fill = "#882255", alpha = .5) +
  geom_vline(xintercept = trace_tmp3$theta3, color="#882255", linetype="dashed", linewidth=0.6) +
  theme_classic() +
  xlab("ln(body mass (kg))") +
  ylab("Density") +
  ggtitle("3-state model") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15), limits=c(0,0.15))
st_dist_3

trace_tmp4 <- trace_sdou4[[1]] %>% 
  summarize(theta1 = median(`theta[1]`),
            theta2 = median(`theta[2]`),
            theta3 = median(`theta[3]`),
            theta4 = median(`theta[4]`),
            stv1 = median(`stv[1]`),
            stv2 = median(`stv[2]`),
            stv3 = median(`stv[3]`),
            stv4 = median(`stv[4]`))

st_dist_4 <- 
  ggplot() +
  # herb
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp4$theta1, sd = sqrt(trace_tmp4$stv1)), color="#44aa99") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp4$theta1, sd=sqrt(trace_tmp4$stv1)), 
                geom = "area", fill = "#44aa99", alpha = .5) +
  geom_vline(xintercept = trace_tmp4$theta1, color="#44aa99", linetype="dashed", linewidth=0.6) +
  # omn1
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp4$theta2, sd = sqrt(trace_tmp4$stv2)), color="#999933") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp4$theta2, sd=sqrt(trace_tmp4$stv2)), 
                geom = "area", fill = "#999933", alpha = .5) +
  geom_vline(xintercept = trace_tmp4$theta2, color="#999933", linetype="dashed", linewidth=0.6) +
  # omn2
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp4$theta3, sd = sqrt(trace_tmp4$stv3)), color="#CC6677") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp4$theta3, sd=sqrt(trace_tmp4$stv3)), 
                geom = "area", fill = "#CC6677", alpha = .5) +
  geom_vline(xintercept = trace_tmp4$theta3, color="#CC6677", linetype="dashed", linewidth=0.6) +
  # carn
  stat_function(data = data.frame(x = c(-20, 30)), aes(x), fun = dnorm, n = 1001, args = list(mean = trace_tmp4$theta4, sd = sqrt(trace_tmp4$stv4)), color="#882255") + ylab("") +
  stat_function(fun = funcShaded, args = list(qLower=0.025, qUpper=0.975, mean=trace_tmp4$theta4, sd=sqrt(trace_tmp4$stv4)), 
                geom = "area", fill = "#882255", alpha = .5) +
  geom_vline(xintercept = trace_tmp4$theta4, color="#882255", linetype="dashed", linewidth=0.6) +
  theme_classic() +
  xlab("ln(body mass (kg))") +
  ylab("Density") +
  ggtitle("4-state model") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_continuous(breaks = c(0,0.05,0.1,0.15), limits=c(0,0.15))
st_dist_4




color_sdou0 <- c("#44aa99", "#ddcc77", "#999933", "#CC6677", "#882255")
p0 <- plotTrace(trace_sdou3, vars = c("theta[1]", "theta[2]", "theta[3]", "stv[1]", "stv[2]"), color = color_sdou0)[[1]] +
  scale_color_manual(values=c("#44aa99", "#ddcc77", "#999933", "#CC6677", "#882255"), 
                     name="Diet",
                     labels=c("Herbivore",
                              "Omnivore (all)",
                              TeX("Omnivore ($>50\\%$ plants)"),
                              TeX("Omnivore ($\\leq 50\\%$ plants)"),
                              "Carnivore"))
legend0 <- get_legend2(p0)

stdist <- cowplot::plot_grid(st_dist_3, st_dist_4, ncol=2)
stdist_wlegend <- plot_grid(stdist, legend0, rel_widths = c(2, .7))

ggsave("figures/3_empirical/sdOU_r500_4StateModel/stationary_distributions.pdf", stdist_wlegend, width = 9, height = 4, units = "in")

