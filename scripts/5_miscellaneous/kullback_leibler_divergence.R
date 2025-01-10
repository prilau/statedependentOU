library(seewave)
library(RevGadgets)
library(tidyverse)
library(ggplot2)

# state-less OU
# run 1
log1 <- readTrace("output/1_validation/artiodactyla_random_diets/1.log", burnin = 0.1)[[1]]
log2 <- readTrace("output/1_validation/artiodactyla_random_diets/2.log", burnin = 0.1)[[1]]
log3 <- readTrace("output/1_validation/artiodactyla_random_diets/3.log", burnin = 0.1)[[1]]
log4 <- readTrace("output/1_validation/artiodactyla_random_diets/4.log", burnin = 0.1)[[1]]
log5 <- readTrace("output/1_validation/artiodactyla_random_diets/5.log", burnin = 0.1)[[1]]
log6 <- readTrace("output/1_validation/artiodactyla_random_diets/6.log", burnin = 0.1)[[1]]
log7 <- readTrace("output/1_validation/artiodactyla_random_diets/7.log", burnin = 0.1)[[1]]
log8 <- readTrace("output/1_validation/artiodactyla_random_diets/8.log", burnin = 0.1)[[1]]
log9 <- readTrace("output/1_validation/artiodactyla_random_diets/9.log", burnin = 0.1)[[1]]
log10 <- readTrace("output/1_validation/artiodactyla_random_diets/10.log", burnin = 0.1)[[1]]
log_emp <- readTrace("output/1_validation/artiodactyla_random_diets/empirical.log", burnin = 0.1)[[1]]


#ou1 <- read.csv("output/1_validation/state_less_OU/state_less_OU_run_1.log", sep="\t") %>%
#  mutate(stationaryVariance = sigma2 / (2 * alpha))
#write_delim(ou1, file="output/1_validation/state_less_OU/state_less_OU_run_1_modified.log", delim="\t")


# sigma2
kl.dist(log_emp$`sigma2[1]`, log1$`sigma2[1]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (91+86)/(91*86))

kl.dist(log_emp$`sigma2[2]`, log1$`sigma2[2]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (47+963)/(47*963))

kl.dist(log_emp$`sigma2[3]`, log1$`sigma2[3]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (20+78)/(20*78))

# alpha
kl.dist(log_emp$`alpha[1]`, log1$`alpha[1]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (234+70)/(234*70))

kl.dist(log_emp$`alpha[2]`, log1$`alpha[2]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (75+110)/(75*110))

kl.dist(log_emp$`alpha[3]`, log1$`alpha[3]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (19+68)/(19*68))

# theta
log_emp$`theta[3]` <- log_emp$`theta[3]` - 12
log1$`theta[3]` <- log1$`theta[3]` - 12

kl.dist(log_emp$`theta[1]`, log1$`theta[1]`, exp(1))
sqrt(-log(0.05/2) * 0.5 * (10+62)/(10*62))

kl.dist(log_emp$`theta[2]`, log1$`theta[2]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (95+127)/(95*127))

kl.dist(log_emp$`theta[3]`, log1$`theta[3]`, exp(1))[[3]]
sqrt(-log(0.05/2) * 0.5 * (21+78)/(21*78))
