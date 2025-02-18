library(phytools)
library(RevGadgets)
library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)
library(latex2exp)
library(ggpmisc)
source("scripts/5_miscellaneous/functions.R")

tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")
t <- read.tree("data/3_empirical/mammal_perMY_r500.tre")
trait <- read.csv("data/3_empirical/mammal_traits.csv") %>% filter(Binomial.1.2 %in% t$tip.label)
#############
#  2-state  #
#############
trace <- readTrace("output/3_empirical/sdOU_r500_missingStateModel/trace_nstate_2_run_6.log", burnin=0.1)
ase <- processAncStates("output/3_empirical/sdOU_r500_missingStateModel/anc_states_nstate_2_run_6.log",
                        state_labels=c("0"="Large", "1"="Small"))
# Plot p0 to for its legend only
p0 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Large"="#364B9a",
                                      "Small"="#c2e4ef"),
                       node_size = c(0.3, 1.2),
                       tip_states = TRUE,
                       tip_states_size = 0.3,
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25)
legend2 <- get_legend2(p0 + theme(legend.position = "left",
                                  legend.box.margin = margin(0, 0, 0, 12))
                       + scale_color_manual(values=c("#364B9a", "#c2e4ef"), 
                                            name="Body size optimum",
                                            labels=c("Large", "Small"))
                       + guides(size = "none",
                                color = guide_legend(override.aes = list(size = 3),
                                                     title="Body size optimum"),
                                fill=guide_legend(title="Body size optimum")))


p1 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Large"="#364B9a", "Small"="#c2e4ef", "2"="grey0", "3"="grey1"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.3,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  #theme(legend.position.inside = c(0.6,0.81))
  theme(legend.position = "none")

est <- as.data.frame(ase@data$anc_state_1[1:500]) %>% rename(est=`ase@data$anc_state_1[1:500]`)
est$Binomial.1.2 <- ase@phylo$tip.label

compare_tips <- trait %>% select(Binomial.1.2, diet4) %>% 
  rename(tru=diet4)
compare_tips <- merge(compare_tips, est, by="Binomial.1.2")

compare_tips <- compare_tips %>%
  group_by(tru, est) %>% 
  summarise(count = n()) %>% 
  pivot_wider(names_from = est, values_from = count) %>% 
  summarise(large = Large/(Large+Small),
            small = Small/(Large+Small)) %>% 
  pivot_longer(cols=2:3, names_to = "est", values_to = "ratio")
#c("Herbivore", TeX("Omnivore ($>50\\%$ plants)"),
#                      TeX("Omnivore ($\\leq 50\\%$ plants)"), "Carnivore")

p2 <- ggplot(compare_tips) +
  geom_col(aes(y=as.factor(tru), x=ratio, group=est, fill=est)) +
  scale_fill_manual(values=c("#364B9a", "#c2e4ef")) +
  scale_y_discrete(labels = c("0"="H", "1"="pO",
                             "2"="npO", "3"="C"), name="True tip state") +
  #ggtitle("Ratio of inferred tip states per diet group") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5)) +
  xlab("Relative frequency of inferred tip states")

p2

post <- as.data.frame(ase@data$anc_state_1_pp[1:500]) %>%
  rename(p=`ase@data$anc_state_1_pp[1:500]`) %>% 
  mutate(p=as.numeric(p))
p3 <- ggplot(post) +
  geom_histogram(aes(x=p), bins = 30, fill="grey50") +
  geom_vline(xintercept = 0.975, linetype="dashed", color="darkred") +
  theme_classic() +
  ylab("") +
  xlab("Posterior probability of MAP tip state")

tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)),
                     rowhead=list(fg_params=list(col=NA)))
ptable <- trace[[1]] %>% mutate(dhalflife12=ifelse(`halflife[1]`>`halflife[2]`, 1, 0)) %>% 
  summarise(Prob=sprintf('%.3f',sum(dhalflife12)/length(dhalflife12))) %>% 
  mutate(#Parameter = "Vy",
    `ij` = "Large; small") %>% 
  select(`ij`, Prob)

colnames(ptable) <- c("Optima i;j", TeX("\\textbf{P(}$t_{0.5_i}>t_{0.5_j}$\\textbf{)}"))

ptable <- tibble(x = 1500, y = 0.015,
                       tb = list(ptable))

color2 <- c("#364B9a", "#c2e4ef")
names(color2) <- c("halflife[1]", "halflife[2]")
p4 <- plotTrace(trace, vars = c("halflife[1]", "halflife[2]"), color = color2)[[1]] +
  theme(legend.position="none") +
  xlim(0,1500) +
  geom_table(data = ptable,
             aes(x = x, y = y, label = tb),
             vjust = 1,
             table.theme = tt) +
  ggtitle(TeX("$t_{0.5}$")) +
  xlab("Time (million years)") +
  ylab("Posterior density")

p4


ptable <- trace[[1]] %>% mutate(dstv12=ifelse(`stv[1]`>`stv[2]`, 1, 0)) %>% 
  summarise(Prob=sprintf('%.3f',sum(dstv12)/length(dstv12))) %>% 
  mutate(#Parameter = "Vy",
    `State i; state j` = "Large; small") %>% 
  select(`State i; state j`, Prob)

colnames(ptable) <- c("Optima i;j", TeX("\\textbf{P(}$V_{y_i}> V_{y_j}$\\textbf{)}"))

ptable <- tibble(x = 30, y = 0.3,
                 tb = list(ptable))

names(color2) <- c("stv[1]", "stv[2]")
p5 <- plotTrace(trace, vars = c("stv[1]", "stv[2]"), color = color2)[[1]] +
  theme(legend.position="none") +
  xlim(0,30) +
  geom_table(data = ptable,
             aes(x = x, y = y, label = tb),
             vjust = 1,
             table.theme = tt) +
  ggtitle(TeX("$V_y$")) +
  xlab("ln(body size (kg))^2") +
  ylab("Posterior density")
p5

names(color2) <- c("theta[1]", "theta[2]")
p6 <- plotTrace(trace, vars = c("theta[1]", "theta[2]"), color = color2)[[1]] +
  #scale_color_manual(values=c("#364B9a", "#4a7bb7", "#c2e4ef"), 
  #                   name="Diet",
  #                   labels=c("Large optimum",
  #                            "Medium optimum",
  #                            "Small optimum")) +
  theme(legend.position="none") +
  ggtitle(TeX("$\\theta$"))  +
  xlab("ln(body size (kg))") +
  ylab("Posterior density")



row1 <- cowplot::plot_grid(p1, legend2, ncol=2, labels = c("(a)",""), rel_widths = c(0.7,0.2), label_size = 12)
row2col1 <- cowplot::plot_grid(p2, p3, ncol=1, labels = c("(b)", "(c)"), rel_heights = c(1,1),
                               vjust = -0.5, label_size = 12)
row3col1 <- cowplot::plot_grid(p6, ncol=1, labels = c("(d)"),
                               vjust = -0.5, label_size = 12)
col1 <- cowplot::plot_grid(row2col1, row3col1, ncol=1, rel_heights = c(1,1))
row23col2 <- cowplot::plot_grid(p4, p5, ncol=1, labels = c("(e)", "(f)"), rel_heights = c(1,1),
                                vjust = -0.5, label_size = 12)
row23 <- cowplot::plot_grid(col1, row23col2, ncol=2, rel_widths = c(1,1))
row123 <- cowplot::plot_grid(row1, row23, ncol=1, rel_heights = c(1, 1))


ggsave("figures/3_empirical/sdOU_r500_MissingStateModel/nstate_2.pdf", row123, width = 7.4, height = 11, units = "in")




#############
#  3-state  #
#############
trace <- readTrace("output/3_empirical/sdOU_r500_missingStateModel/trace_nstate_3_run_5.log", burnin=0.1)
ase <- processAncStates("output/3_empirical/sdOU_r500_missingStateModel/anc_states_nstate_3_run_5.log",
                        state_labels=c("0"="Large", "1"="Medium", "2"="Small"))

# Plot p0 to for its legend only
p0 <- plotAncStatesMAP(t = ase,
                       #tip_labels_offset = 0.5,
                       tip_labels = FALSE,
                       node_color_as = "state",
                       node_color = c("Large"="#364B9a", "Medium"="#4a7bb7", "Small"="#c2e4ef"),
                       node_size = c(0.3, 1.2),
                       tip_states = TRUE,
                       tip_states_size = 0.3,
                       #tip_states_shape = 1,
                       state_transparency = 0.7,
                       tree_layout = "circular",
                       #tip_labels_size = 0.5
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25)
legend3 <- get_legend2(p0 + theme(legend.position = "left",
                                  legend.box.margin = margin(0, 0, 0, 12))
                       + scale_color_manual(values=c("#364B9a", "#4a7bb7", "#c2e4ef"), 
                                            name="Body size optimum",
                                            labels=c("Large", "Medium", "Small"))
                       + guides(size = "none",
                                color = guide_legend(override.aes = list(size = 3),
                                                     title="Body size optimum"),
                                fill=guide_legend(title="Body size optimum")))


p7 <- plotAncStatesPie(t = ase,
                       pie_colors = c("Large"="#364B9a", "Medium"="#4a7bb7", "Small"="#c2e4ef", "3"="grey1"),
                       #tip_labels_size = 1,
                       tip_pies = TRUE,
                       #tip_labels_offset = 0.5,
                       node_pie_size = 0.7,
                       tip_pie_size = 0.3,
                       tree_layout = "circular",
                       tip_labels = FALSE,
                       state_transparency = 0.7,
                       tree_color = "#bbbbbb",
                       tree_linewidth = 0.25) +
  # modify legend location using ggplot2
  #theme(legend.position.inside = c(0.6,0.81))
  theme(legend.position = "none")

est <- as.data.frame(ase@data$anc_state_1[1:500]) %>% rename(est=`ase@data$anc_state_1[1:500]`)
est$Binomial.1.2 <- ase@phylo$tip.label

compare_tips <- trait %>% select(Binomial.1.2, diet4) %>% 
  rename(tru=diet4)
compare_tips <- merge(compare_tips, est, by="Binomial.1.2")

compare_tips <- compare_tips %>%
  group_by(tru, est) %>% 
  summarise(count = n()) %>% 
  pivot_wider(names_from = est, values_from = count) %>% 
  summarise(large = Large/(Large+Medium+Small),
            medium = Medium/(Large+Medium+Small),
            small = Small/(Large+Medium+Small)) %>% 
  pivot_longer(cols=2:4, names_to = "est", values_to = "ratio")

p8 <- ggplot(compare_tips) +
  geom_col(aes(y=as.factor(tru), x=ratio, group=est, fill=est)) +
  scale_fill_manual(values=c("#364B9a", "#4a7bb7", "#c2e4ef")) +
  scale_y_discrete(labels = c("0"="H", "1"="pO",
                              "2"="npO", "3"="C"), name="True tip state") +
  #ggtitle("Ratio of inferred tip states per diet group") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(hjust=0.5)) +
  xlab("Relative frequency of inferred tip states")


post <- as.data.frame(ase@data$anc_state_1_pp[1:500]) %>%
  rename(p=`ase@data$anc_state_1_pp[1:500]`) %>% 
  mutate(p=as.numeric(p))
p9 <- ggplot(post) +
  geom_histogram(aes(x=p), bins = 30, fill="grey50") +
  geom_vline(xintercept = 0.975, linetype="dashed", color="darkred") +
  theme_classic() +
  ylab("") +
  xlab("Posterior probability of MAP tip state")

ptable <- trace[[1]] %>% mutate(dhalflife12=ifelse(`halflife[1]`>`halflife[2]`, 1, 0),
                                dhalflife13=ifelse(`halflife[1]`>`halflife[3]`, 1, 0),
                                dhalflife23=ifelse(`halflife[2]`>`halflife[3]`, 1, 0)) %>% 
  summarise(phalflife12=sprintf('%.3f',sum(dhalflife12)/length(dhalflife12)),
            phalflife13=sprintf('%.3f',sum(dhalflife13)/length(dhalflife13)),
            phalflife23=sprintf('%.3f',sum(dhalflife23)/length(dhalflife23)))
ptable <- data.frame(n=c("Large; medium", "Large; small", "Medium; small"), p=c(ptable$phalflife12, ptable$phalflife13, ptable$phalflife23))

colnames(ptable) <- c("Optima i;j", TeX("\\textbf{P(}$t_{0.5_i}>t_{0.5_j}$\\textbf{)}"))
ptable <- tibble(x = 2000, y = 0.018,
                 tb = list(ptable))

color3 <- c("#364B9a", "#4a7bb7", "#c2e4ef")
names(color3) <- c("halflife[1]", "halflife[2]", "halflife[3]")
p10 <- plotTrace(trace, vars = c("halflife[1]", "halflife[2]", "halflife[3]"), color = color3)[[1]] +
  theme(legend.position="none") +
  geom_table(data = ptable,
             aes(x = x, y = y, label = tb),
             vjust = 1,
             table.theme = tt) +
  ggtitle(TeX("$t_{0.5}$")) +
  xlab("Time (million years)") +
  ylab("Posterior density")

p10


ptable <- trace[[1]] %>% mutate(dstv12=ifelse(`stv[1]`>`stv[2]`, 1, 0),
                                dstv13=ifelse(`stv[1]`>`stv[3]`, 1, 0),
                                dstv23=ifelse(`stv[2]`>`stv[3]`, 1, 0)) %>% 
  summarise(pstv12=sprintf('%.3f',sum(dstv12)/length(dstv12)),
            pstv13=sprintf('%.3f',sum(dstv13)/length(dstv13)),
            pstv23=sprintf('%.3f',sum(dstv23)/length(dstv23)))
ptable <- data.frame(n=c("Large; medium", "Large; small", "Medium; small"), p=c(ptable$pstv12, ptable$pstv13, ptable$pstv23))
colnames(ptable) <- c("Optima i;j", TeX("\\textbf{P(}$t_{0.5_i}>t_{0.5_j}$\\textbf{)}"))
ptable <- tibble(x = 40, y = 0.45,
                 tb = list(ptable))

names(color3) <- c("stv[1]", "stv[2]", "stv[3]")
p11 <- plotTrace(trace, vars = c("stv[1]", "stv[2]", "stv[3]"), color = color3)[[1]] +
  theme(legend.position="none") +
  geom_table(data = ptable,
             aes(x = x, y = y, label = tb),
             vjust = 1,
             table.theme = tt) +
  ggtitle(TeX("$V_y$")) +
  xlab("ln(body size (kg))^2") +
  ylab("Posterior density")


names(color3) <- c("theta[1]", "theta[2]", "theta[3]")
p12 <- plotTrace(trace, vars = c("theta[1]", "theta[2]", "theta[3]"), color = color3)[[1]] +
  theme(legend.position="none") +
  ggtitle(TeX("$\\theta$"))  +
  xlab("ln(body size (kg))") +
  ylab("Posterior density")

p12


row1 <- cowplot::plot_grid(p7, legend3, ncol=2, labels = c("(a)",""), rel_widths = c(0.7,0.2), label_size = 12)
row2col1 <- cowplot::plot_grid(p8, p9, ncol=1, labels = c("(b)", "(c)"), rel_heights = c(1,1),
                               vjust = -0.5, label_size = 12)
row3col1 <- cowplot::plot_grid(p12, ncol=1, labels = c("(d)"),
                               vjust = -0.5, label_size = 12)
col1 <- cowplot::plot_grid(row2col1, row3col1, ncol=1, rel_heights = c(1,1))
row23col2 <- cowplot::plot_grid(p10, p11, ncol=1, labels = c("(e)", "(f)"), rel_heights = c(1,1),
                                vjust = -0.5, label_size = 12)
row23 <- cowplot::plot_grid(col1, row23col2, ncol=2, rel_widths = c(1,1))
row123 <- cowplot::plot_grid(row1, row23, ncol=1, rel_heights = c(1, 1))

ggsave("figures/3_empirical/sdOU_r500_MissingStateModel/nstate_3.pdf", row123, width = 7.4, height = 11, units = "in")


