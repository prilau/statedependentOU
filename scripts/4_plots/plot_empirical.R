library(RevGadgets)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(latex2exp)


######################
# sdOU 4-state model #
######################

# specify input and output directories
dir_in="output/3_empirical/sdOU_r500_4StateModel/"
dir_out="figures/3_empirical/sdOU_r500_4StateModel/"

runs = c(1:10)
probs_4 <- tibble(run = runs,
                ptheta.12 = 0, ptheta.13 = 0, ptheta.14 = 0,  ptheta.23 = 0, ptheta.24 = 0, ptheta.34 = 0,   
                pstv.12 = 0, pstv.13 = 0, pstv.14 = 0, pstv.23 = 0, pstv.24 = 0, pstv.34 = 0,
                phalflife.12 = 0, phalflife.13 = 0, phalflife.14 = 0, phalflife.23 = 0, phalflife.24 = 0, phalflife.34 = 0,
                psigma2.12 = 0, psigma2.13 = 0, psigma2.14 = 0, psigma2.23 = 0, psigma2.24 = 0, psigma2.34 = 0)

bar = txtProgressBar(style=3, width=40)
for (i in 1:length(runs)){
  filename <- paste0(dir_in, list.files(dir_in, pattern=paste0("trace_run_", runs[i])))
  probs_4[i,2:25] <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select(`theta[1]`, `theta[2]`, `theta[3]`, `theta[4]`,
           `stv[1]`, `stv[2]`, `stv[3]`, `stv[4]`,
           `halflife[1]`, `halflife[2]`, `halflife[3]`, `halflife[4]`,
           `sigma2[1]`, `sigma2[2]`, `sigma2[3]`, `sigma2[4]`) %>% 
    mutate(dtheta.12 = ifelse(`theta[1]` > `theta[2]`, 1, 0),
           dtheta.13 = ifelse(`theta[1]` > `theta[3]`, 1, 0),
           dtheta.14 = ifelse(`theta[1]` > `theta[4]`, 1, 0),
           dtheta.23 = ifelse(`theta[2]` > `theta[3]`, 1, 0),
           dtheta.24 = ifelse(`theta[2]` > `theta[4]`, 1, 0),
           dtheta.34 = ifelse(`theta[3]` > `theta[4]`, 1, 0),
           
           dstv.12 = ifelse(`stv[1]` > `stv[2]`, 1, 0),
           dstv.13 = ifelse(`stv[1]` > `stv[3]`, 1, 0),
           dstv.14 = ifelse(`stv[1]` > `stv[4]`, 1, 0),
           dstv.23 = ifelse(`stv[2]` > `stv[3]`, 1, 0),
           dstv.24 = ifelse(`stv[2]` > `stv[4]`, 1, 0),
           dstv.34 = ifelse(`stv[3]` > `stv[4]`, 1, 0),
           
           dhalflife.12 = ifelse(`halflife[1]` > `halflife[2]`, 1, 0),
           dhalflife.13 = ifelse(`halflife[1]` > `halflife[3]`, 1, 0),
           dhalflife.14 = ifelse(`halflife[1]` > `halflife[4]`, 1, 0),
           dhalflife.23 = ifelse(`halflife[2]` > `halflife[3]`, 1, 0),
           dhalflife.24 = ifelse(`halflife[2]` > `halflife[4]`, 1, 0),
           dhalflife.34 = ifelse(`halflife[3]` > `halflife[4]`, 1, 0),
           
           dsigma2.12 = ifelse(`sigma2[1]` > `sigma2[2]`, 1, 0),
           dsigma2.13 = ifelse(`sigma2[1]` > `sigma2[3]`, 1, 0),
           dsigma2.14 = ifelse(`sigma2[1]` > `sigma2[4]`, 1, 0),
           dsigma2.23 = ifelse(`sigma2[2]` > `sigma2[3]`, 1, 0),
           dsigma2.24 = ifelse(`sigma2[2]` > `sigma2[4]`, 1, 0),
           dsigma2.34 = ifelse(`sigma2[3]` > `sigma2[4]`, 1, 0)
           ) %>% 
    summarise(ptheta.12 = mean(dtheta.12),
              ptheta.13 = mean(dtheta.13),
              ptheta.14 = mean(dtheta.14),
              ptheta.23 = mean(dtheta.23),
              ptheta.24 = mean(dtheta.24),
              ptheta.34 = mean(dtheta.34),
              
              pstv.12 = mean(dstv.12),
              pstv.13 = mean(dstv.13),
              pstv.14 = mean(dstv.14),
              pstv.23 = mean(dstv.23),
              pstv.24 = mean(dstv.24),
              pstv.34 = mean(dstv.34),
              
              phalflife.12 = mean(dhalflife.12),
              phalflife.13 = mean(dhalflife.13),
              phalflife.14 = mean(dhalflife.14),
              phalflife.23 = mean(dhalflife.23),
              phalflife.24 = mean(dhalflife.24),
              phalflife.34 = mean(dhalflife.34),
              
              psigma2.12 = mean(dsigma2.12),
              psigma2.13 = mean(dsigma2.13),
              psigma2.14 = mean(dsigma2.14),
              psigma2.23 = mean(dsigma2.23),
              psigma2.24 = mean(dsigma2.24),
              psigma2.34 = mean(dsigma2.34)
              ) %>% 
    #unlist() %>% 
    unname()
  setTxtProgressBar(bar, i / length(runs))
}

probs_4_sum <- probs_4 %>% summarise(ptheta.12 = mean(ptheta.12),
                          ptheta.13 = mean(ptheta.13),
                          ptheta.14 = mean(ptheta.14),
                          ptheta.23 = mean(ptheta.23),
                          ptheta.24 = mean(ptheta.24),
                          ptheta.34 = mean(ptheta.34),
                          
                          pstv.12 = mean(pstv.12),
                          pstv.13 = mean(pstv.13),
                          pstv.14 = mean(pstv.14),
                          pstv.23 = mean(pstv.23),
                          pstv.24 = mean(pstv.24),
                          pstv.34 = mean(pstv.34),
                          
                          phalflife.12 = mean(phalflife.12),
                          phalflife.13 = mean(phalflife.13),
                          phalflife.14 = mean(phalflife.14),
                          phalflife.23 = mean(phalflife.23),
                          phalflife.24 = mean(phalflife.24),
                          phalflife.34 = mean(phalflife.34),
                          
                          psigma2.12 = mean(psigma2.12),
                          psigma2.13 = mean(psigma2.13),
                          psigma2.14 = mean(psigma2.14),
                          psigma2.23 = mean(psigma2.23),
                          psigma2.24 = mean(psigma2.24),
                          psigma2.34 = mean(psigma2.34))

######################
# sdOU 3-state model #
######################

# specify input and output directories
dir_in="output/3_empirical/sdOU_r500_3StateOrderedModel/"
dir_out="figures/3_empirical/sdOU_r500_3StateOrderedModel/"

#num_run = length(list.files(dir_in, pattern=".log"))
runs = c(1:10)
probs_3 <- tibble(run = runs,
                ptheta.12 = 0, ptheta.13 = 0, ptheta.23 = 0,
                pstv.12 = 0, pstv.13 = 0, pstv.23 = 0,
                phalflife.12 = 0, phalflife.13 = 0, phalflife.23 = 0,
                psigma2.12 = 0, psigma2.13 = 0, psigma2.23 = 0)

bar = txtProgressBar(style=3, width=40)
for (i in 1:length(runs)){
  filename <- paste0(dir_in, list.files(dir_in, pattern=paste0("trace_run_", runs[i])))
  probs_3[i,2:13] <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select(`theta[1]`, `theta[2]`, `theta[3]`,
           `stv[1]`, `stv[2]`, `stv[3]`,
           `halflife[1]`, `halflife[2]`, `halflife[3]`,
           `sigma2[1]`, `sigma2[2]`, `sigma2[3]`) %>% 
    mutate(dtheta.12 = ifelse(`theta[1]` > `theta[2]`, 1, 0),
           dtheta.13 = ifelse(`theta[1]` > `theta[3]`, 1, 0),
           dtheta.23 = ifelse(`theta[2]` > `theta[3]`, 1, 0),
           
           dstv.12 = ifelse(`stv[1]` > `stv[2]`, 1, 0),
           dstv.13 = ifelse(`stv[1]` > `stv[3]`, 1, 0),
           dstv.23 = ifelse(`stv[2]` > `stv[3]`, 1, 0),
           
           dhalflife.12 = ifelse(`halflife[1]` > `halflife[2]`, 1, 0),
           dhalflife.13 = ifelse(`halflife[1]` > `halflife[3]`, 1, 0),
           dhalflife.23 = ifelse(`halflife[2]` > `halflife[3]`, 1, 0),
           
           dsigma2.12 = ifelse(`sigma2[1]` > `sigma2[2]`, 1, 0),
           dsigma2.13 = ifelse(`sigma2[1]` > `sigma2[3]`, 1, 0),
           dsigma2.23 = ifelse(`sigma2[2]` > `sigma2[3]`, 1, 0)
    ) %>% 
    summarise(ptheta.12 = mean(dtheta.12),
              ptheta.13 = mean(dtheta.13),
              ptheta.23 = mean(dtheta.23),
              
              pstv.12 = mean(dstv.12),
              pstv.13 = mean(dstv.13),
              pstv.23 = mean(dstv.23),
              
              phalflife.12 = mean(dhalflife.12),
              phalflife.13 = mean(dhalflife.13),
              phalflife.23 = mean(dhalflife.23),
              
              psigma2.12 = mean(dsigma2.12),
              psigma2.13 = mean(dsigma2.13),
              psigma2.23 = mean(dsigma2.23)
    ) %>% 
    #unlist() %>% 
    unname()
  setTxtProgressBar(bar, i / length(runs))
}

probs_3_sum <- probs %>% summarise(ptheta.12 = mean(ptheta.12),
                                   ptheta.13 = mean(ptheta.13),
                                   ptheta.23 = mean(ptheta.23),
                                   
                                   pstv.12 = mean(pstv.12),
                                   pstv.13 = mean(pstv.13),
                                   pstv.23 = mean(pstv.23),
                                   
                                   phalflife.12 = mean(phalflife.12),
                                   phalflife.13 = mean(phalflife.13),
                                   phalflife.23 = mean(phalflife.23),
                                   
                                   psigma2.12 = mean(psigma2.12),
                                   psigma2.13 = mean(psigma2.13),
                                   psigma2.23 = mean(psigma2.23))


###########################
# sdOU hidden-state model #
###########################

# specify input and output directories
dir_in="output/3_empirical/sdOU_r500_hiddenStateModel/"
dir_out="figures/3_empirical/sdOU_r500_hiddenStateModel/"

runs = c(1:20)
probs_h <- tibble(run = runs,
                ptheta.12 = 0, ptheta.13 = 0, ptheta.15 = 0, ptheta.23 = 0, ptheta.25 = 0, ptheta.35 = 0,   
                pstv.12 = 0, pstv.13 = 0, pstv.15 = 0, pstv.23 = 0, pstv.25 = 0, pstv.35 = 0,
                phalflife.12 = 0, phalflife.13 = 0, phalflife.15 = 0, phalflife.23 = 0, phalflife.25 = 0, phalflife.35 = 0,
                psigma2.12 = 0, psigma2.13 = 0, psigma2.15 = 0, psigma2.23 = 0, psigma2.25 = 0, psigma2.35 = 0)

bar = txtProgressBar(style=3, width=40)
for (i in 1:length(runs)){
  filename <- paste0(dir_in, list.files(dir_in, pattern=paste0("trace_run_", runs[i])))
  probs_h[i,2:25] <- readTrace(filename, burnin = 0.0)[[1]] %>% 
    select(`theta[1]`, `theta[2]`, `theta[3]`, `theta[5]`,
           `stv[1]`, `stv[2]`, `stv[3]`, `stv[5]`,
           `halflife[1]`, `halflife[2]`, `halflife[3]`, `halflife[5]`,
           `sigma2[1]`, `sigma2[2]`, `sigma2[3]`, `sigma2[5]`) %>% 
    mutate(dtheta.12 = ifelse(`theta[1]` > `theta[2]`, 1, 0),
           dtheta.13 = ifelse(`theta[1]` > `theta[3]`, 1, 0),
           dtheta.15 = ifelse(`theta[1]` > `theta[5]`, 1, 0),
           dtheta.23 = ifelse(`theta[2]` > `theta[3]`, 1, 0),
           dtheta.25 = ifelse(`theta[2]` > `theta[5]`, 1, 0),
           dtheta.35 = ifelse(`theta[3]` > `theta[5]`, 1, 0),
           
           dstv.12 = ifelse(`stv[1]` > `stv[2]`, 1, 0),
           dstv.13 = ifelse(`stv[1]` > `stv[3]`, 1, 0),
           dstv.15 = ifelse(`stv[1]` > `stv[5]`, 1, 0),
           dstv.23 = ifelse(`stv[2]` > `stv[3]`, 1, 0),
           dstv.25 = ifelse(`stv[2]` > `stv[5]`, 1, 0),
           dstv.35 = ifelse(`stv[3]` > `stv[5]`, 1, 0),
           
           dhalflife.12 = ifelse(`halflife[1]` > `halflife[2]`, 1, 0),
           dhalflife.13 = ifelse(`halflife[1]` > `halflife[3]`, 1, 0),
           dhalflife.15 = ifelse(`halflife[1]` > `halflife[5]`, 1, 0),
           dhalflife.23 = ifelse(`halflife[2]` > `halflife[3]`, 1, 0),
           dhalflife.25 = ifelse(`halflife[2]` > `halflife[5]`, 1, 0),
           dhalflife.35 = ifelse(`halflife[3]` > `halflife[5]`, 1, 0),
           
           dsigma2.12 = ifelse(`sigma2[1]` > `sigma2[2]`, 1, 0),
           dsigma2.13 = ifelse(`sigma2[1]` > `sigma2[3]`, 1, 0),
           dsigma2.15 = ifelse(`sigma2[1]` > `sigma2[5]`, 1, 0),
           dsigma2.23 = ifelse(`sigma2[2]` > `sigma2[3]`, 1, 0),
           dsigma2.25 = ifelse(`sigma2[2]` > `sigma2[5]`, 1, 0),
           dsigma2.35 = ifelse(`sigma2[3]` > `sigma2[5]`, 1, 0)
    ) %>% 
    summarise(ptheta.12 = mean(dtheta.12),
              ptheta.13 = mean(dtheta.13),
              ptheta.15 = mean(dtheta.15),
              ptheta.23 = mean(dtheta.23),
              ptheta.25 = mean(dtheta.25),
              ptheta.35 = mean(dtheta.35),
              
              pstv.12 = mean(dstv.12),
              pstv.13 = mean(dstv.13),
              pstv.15 = mean(dstv.15),
              pstv.23 = mean(dstv.23),
              pstv.25 = mean(dstv.25),
              pstv.35 = mean(dstv.35),
              
              phalflife.12 = mean(dhalflife.12),
              phalflife.13 = mean(dhalflife.13),
              phalflife.15 = mean(dhalflife.15),
              phalflife.23 = mean(dhalflife.23),
              phalflife.25 = mean(dhalflife.25),
              phalflife.35 = mean(dhalflife.35),
              
              psigma2.12 = mean(dsigma2.12),
              psigma2.13 = mean(dsigma2.13),
              psigma2.15 = mean(dsigma2.15),
              psigma2.23 = mean(dsigma2.23),
              psigma2.25 = mean(dsigma2.25),
              psigma2.35 = mean(dsigma2.35)
    ) %>% 
    #unlist() %>% 
    unname()
  setTxtProgressBar(bar, i / length(runs))
}
probs_h_sum <- probs_h %>% summarise(ptheta.12 = mean(ptheta.12),
                              ptheta.13 = mean(ptheta.13),
                              ptheta.15 = mean(ptheta.15),
                              ptheta.23 = mean(ptheta.23),
                              ptheta.25 = mean(ptheta.25),
                              ptheta.35 = mean(ptheta.35),
                              
                              pstv.12 = mean(pstv.12),
                              pstv.13 = mean(pstv.13),
                              pstv.15 = mean(pstv.15),
                              pstv.23 = mean(pstv.23),
                              pstv.25 = mean(pstv.25),
                              pstv.35 = mean(pstv.35),
                              
                              phalflife.12 = mean(phalflife.12),
                              phalflife.13 = mean(phalflife.13),
                              phalflife.15 = mean(phalflife.15),
                              phalflife.23 = mean(phalflife.23),
                              phalflife.25 = mean(phalflife.25),
                              phalflife.35 = mean(phalflife.35),
                              
                              psigma2.12 = mean(psigma2.12),
                              psigma2.13 = mean(psigma2.13),
                              psigma2.15 = mean(psigma2.15),
                              psigma2.23 = mean(psigma2.23),
                              psigma2.25 = mean(psigma2.25),
                              psigma2.35 = mean(psigma2.35))

################
# df for plots #
################
powers_th <- grid_th  %>%
  mutate(sig = ifelse(prob >= 0.975 | prob <= 0.025, 1, 0),
         sig_T = ifelse(prob >= 0.975 & dtheta > 0 | prob <= 0.025 & dtheta < 0, 1, 0)) %>% 
  group_by(combo) %>% 
  summarise(power = mean(sig_T),
            ppv = sum(sig_T)/sum(sig),
            num_ppv = sum(sig)) %>% 
  mutate(halflife = ifelse(combo %in% 1:3, 0.1,
                           ifelse(combo %in% 4:6, 0.3,
                                  0.6)),
         stv = ifelse(combo %in% c(1,4,7), "0.5 v",
                      ifelse(combo %in% c(2, 5, 8), "1 v",
                             "2 v")),
         q025_power = qnorm(0.025, mean=power, sd = sqrt(power*(1-power)/200)), # sqrt(p(1-p)/n) is 1 standard error
         q975_power = qnorm(0.975, mean=power, sd = sqrt(power*(1-power)/200)),
         q025_ppv = qnorm(0.025, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)),
         q975_ppv = qnorm(0.975, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)))


powers_hl <- grid_hl  %>%
  mutate(sig = ifelse(prob_hl >= 0.975 | prob_hl <= 0.025, 1, 0),
         sig_T = ifelse(prob_hl >= 0.975 & dhalflife > 0 | prob_hl <= 0.025 & dhalflife < 0, 1, 0)) %>% 
  group_by(combo) %>% 
  summarise(power = mean(sig_T),
            ppv = sum(sig_T)/sum(sig),
            num_ppv = sum(sig)) %>% 
  mutate(dtheta_T = ifelse(combo %in% 1:3, 2,
                           ifelse(combo %in% 4:6, 6,
                                  10)),
         stv = ifelse(combo %in% c(1,4,7), 0.5,
                      ifelse(combo %in% c(2, 5, 8), 1,
                             2)),
         q025_power = qnorm(0.025, mean=power, sd = sqrt(power*(1-power)/200)), # sqrt(p(1-p)/n) is 1 standard error
         q975_power = qnorm(0.975, mean=power, sd = sqrt(power*(1-power)/200)),
         q025_ppv = qnorm(0.025, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)),
         q975_ppv = qnorm(0.975, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)))

powers_stv <- grid_stv  %>%
  mutate(sig = ifelse(prob_stv >= 0.975 | prob_stv <= 0.025, 1, 0),
         sig_T = ifelse(prob_stv >= 0.975 & dstv_T > 0 | prob_stv <= 0.025 & dstv_T < 0, 1, 0)) %>% 
  group_by(combo) %>% 
  summarise(power = mean(sig_T),
            ppv = sum(sig_T)/sum(sig),
            num_ppv = sum(sig)) %>% 
  mutate(dtheta_T = ifelse(combo %in% 1:3, 2,
                           ifelse(combo %in% 4:6, 6,
                                  10)),
         halflife = ifelse(combo %in% c(1,4,7), 0.1,
                      ifelse(combo %in% c(2, 5, 8), 0.3,
                             0.6)),
         q025_power = qnorm(0.025, mean=power, sd = sqrt(power*(1-power)/200)), # sqrt(p(1-p)/n) is 1 standard error
         q975_power = qnorm(0.975, mean=power, sd = sqrt(power*(1-power)/200)),
         q025_ppv = qnorm(0.025, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)),
         q975_ppv = qnorm(0.975, mean=ppv, sd = sqrt(power*(1-ppv)/num_ppv)))


#########
# plots #
#########
p_th_power <- ggplot(powers_th) +
  geom_point(aes(x = halflife, y=power, group=stv, shape=stv, color=stv),
             size=1.5, alpha=0.9) +
  geom_line(aes(x = halflife, y=power, group=stv, color=stv),
            alpha=0.9) +
  geom_errorbar(aes(x = halflife, ymin = q025_power, ymax=q975_power,
                    group=stv, color=stv, width=0.0125),
                linewidth = 0.25, alpha=0.9)+
  scale_shape_manual("V", values=c(15, 16, 17))+
  scale_color_manual("V",values=c('#004488', '#DDAA33', '#BB5566')) +
  theme_bw() +
  scale_x_continuous(breaks=c(0.1, 0.3, 0.6)) + 
  xlab(TeX("$t_{0.5}$")) +
  ylab("") +
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(TeX("Optimum $\\theta$")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

p_hl_power <- ggplot(powers_hl) +
  geom_point(aes(x = stv, y=power, group=stv,
                 shape=factor(dtheta_T), color=factor(dtheta_T)),
             size=1.5, alpha=0.9) +
  geom_line(aes(x = stv, y=power, group=factor(dtheta_T), color=factor(dtheta_T)),
            alpha=0.9) +
  theme_bw() +
  scale_shape_manual(TeX("$\\Delta \\theta$"), values=c(15, 16, 17)) +
  scale_color_manual(TeX("$\\Delta \\theta$"), values=c('#004488', '#DDAA33', '#BB5566')) +
  geom_errorbar(aes(x = stv, ymin = q025_power, ymax=q975_power,
                    group=factor(dtheta_T), color=factor(dtheta_T), width=0.0375),
                linewidth = 0.25, alpha=0.9)+
  ylab("Power") +
  xlab("V") +
  scale_x_continuous(breaks=c(0.5, 1, 2), label=c("0.5v", "v", "2v")) + 
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle(TeX("Phylogenetic half-life $t_{0.5}$")) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

p_stv_power <- ggplot(powers_stv) +
  geom_point(aes(x = dtheta_T, y=power, group=factor(halflife),
                 shape=factor(halflife), color=factor(halflife)),
             size=1.5, alpha=0.9) +
  geom_line(aes(x = dtheta_T, y=power, group=factor(halflife), color=factor(halflife)),
            alpha=0.9) +
  theme_bw() +
  scale_shape_manual(TeX("$t_{0.5}$"), values=c(15, 16, 17)) +
  scale_color_manual(TeX("$t_{0.5}$"), values=c('#004488', '#DDAA33', '#BB5566')) +
  geom_errorbar(aes(x = dtheta_T, ymin = q025_power, ymax=q975_power,
                    group=halflife, color=factor(halflife), width=0.2),
                linewidth = 0.25, alpha=0.9) +
  ylab("") +
  xlab(TeX("$\\Delta \\theta$")) + 
  scale_x_continuous(breaks=c(2, 6, 10)) + 
  coord_cartesian(ylim=c(0, 1)) +
  ggtitle("Stationary variance V") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
p_stv_power
p_power <- arrangeGrob(p_hl_power, p_stv_power, p_th_power, nrow = 1)
ggsave("figures/2_simulation/power.pdf", p_power, width = 200, height = 90, unit = "mm")

#p_th_ppv <- ggplot(powers_th) +
#  geom_point(aes(x = halflife, y=ppv, group=stv, shape=stv, color=stv), size=2) +
#  geom_line(aes(x = halflife, y=ppv, group=stv, color=stv)) +
#  geom_errorbar(aes(x = halflife, ymin = q025_ppv, ymax=q975_ppv, group=stv, color=stv))+
#  scale_shape_manual("V", values=c(0, 1, 5))+
#  scale_color_manual("V",values=c('#004488', '#DDAA33', '#BB5566')) +
#  scale_x_continuous(breaks=c(0.1, 0.3, 0.6)) + 
#  theme_bw() +
#  theme(legend) +
#  coord_cartesian(ylim=c(0, 1)) +
#  labs(x = TeX("$t_{0.5}$"), y = "#true significant / #significant") + 
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.position = "bottom")



#p_hl_ppv <- ggplot(powers_hl) +
#  geom_point(aes(x = stv, y=ppv, group=stv,shape=factor(dtheta_T), color=factor(dtheta_T)), size=2) +
#  geom_line(aes(x = stv, y=ppv, group=factor(dtheta_T), color=factor(dtheta_T))) +
#  theme_bw() +
#  geom_errorbar(aes(x = stv, ymin = q025_ppv, ymax=q975_ppv, group=factor(dtheta_T), color=factor(dtheta_T)))+
#  scale_shape_manual(TeX("$\\Delta \\theta$"), values=c(0, 1, 5)) +
#  scale_color_manual(TeX("$\\Delta \\theta$"), values=c('#004488', '#DDAA33', '#BB5566')) +
#  scale_x_continuous(breaks=c(0.5, 1, 2), label=c("0.5v", "v", "2v")) + 
#  coord_cartesian(ylim=c(0, 1)) +
#  labs(x = "V", y = "#true significant / #significant") + 
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.position = "bottom")

#p_stv_ppv <- ggplot(powers_stv) +
#  geom_point(aes(x = dtheta_T, y=ppv, group=factor(halflife), shape=factor(halflife), color=factor(halflife)), size=2) +
#  geom_line(aes(x = dtheta_T, y=ppv, group=factor(halflife), color=factor(halflife))) +
#  theme_bw() +
#  geom_errorbar(aes(x = dtheta_T, ymin = q025_ppv, ymax=q975_ppv, group=factor(halflife), color=factor(halflife)))+
#  scale_shape_manual(TeX("$t_{0.5}$"), values=c(0, 1, 5)) +
#  scale_color_manual(TeX("$t_{0.5}$"), values=c('#004488', '#DDAA33', '#BB5566')) +
#  geom_errorbar(aes(x = dtheta_T, ymin = q025_ppv, ymax=q975_ppv, group=halflife, color=factor(halflife)))+
#  scale_x_continuous(breaks=c(2, 6, 10)) + 
#  coord_cartesian(ylim=c(0, 1)) +
#  labs(x = TeX("$\\Delta \\theta$"), y = "#true significant / #significant") + 
#  theme(plot.title = element_text(hjust = 0.5),
#        legend.position = "bottom")

#p_all <- arrangeGrob(p_hl_power, p_stv_power, p_th_power, p_hl_ppv, p_stv_ppv, p_th_ppv, nrow = 2)
#ggsave("figures/2_simulation/power_theta/power_ppv.pdf", p_all, width = 240, height = 180, unit = "mm")


