library(seewave)
library(RevGadgets)
library(tidyverse)
library(ggplot2)

# state-less OU
# run 1
sdou1 <- readTrace("output/1_validation/state_less_OU/state_dependent_OU_run_1.log")[[1]]
sdou1 <- sdou1 %>% mutate(stationaryVariance = sigma2 / (2 * alpha))
ou1 <- readTrace("output/1_validation/state_less_OU/state_less_OU_run_1.log")[[1]]
ou1 <- ou1 %>% mutate(stationaryVariance = sigma2 / (2 * alpha))

kl.dist(sdou1$sigma2, ou1$sigma2, exp(1))
kl.dist(sdou1$alpha, ou1$alpha, exp(1))
kl.dist(sdou1$theta, ou1$theta, exp(1))
kl.dist(sdou1$stationaryVariance, ou1$stationaryVariance, exp(1))

# run 2
sdou2 <- readTrace("output/1_validation/state_less_OU/state_dependent_OU_run_2.log")[[1]]
sdou2 <- sdou2 %>% mutate(stationaryVariance = sigma2 / (2 * alpha))
ou2 <- readTrace("output/1_validation/state_less_OU/state_less_OU_run_2.log")[[1]]
ou2 <- ou2 %>% mutate(stationaryVariance = sigma2 / (2 * alpha))

kl.dist(sdou2$sigma2, ou2$sigma2, exp(1))
kl.dist(sdou2$alpha, ou2$alpha, exp(1))
kl.dist(sdou2$theta, ou2$theta, exp(1))
kl.dist(sdou2$stationaryVariance, ou2$stationaryVariance, exp(1))


# state-less BM
# run 1
sdou3 <- readTrace("output/1_validation/state_less_BM/state_dependent_OU_run_1.log")[[1]]
bm3 <- readTrace("output/1_validation/state_less_BM/state_less_BM_run_1.log")[[1]]

kl.dist(sdou3$sigma2, bm3$sigma2, exp(1))

# run 2
sdou4 <- readTrace("output/1_validation/state_less_BM/state_dependent_OU_run_2.log")[[1]]
bm4 <- readTrace("output/1_validation/state_less_BM/state_less_BM_run_2.log")[[1]]

kl.dist(sdou4$sigma2, bm4$sigma2, exp(1))


# state-dependent BM
# run 1

kl_sdbm <- tibble(model=c("sdou", "sdbm", "sdbm_maymoore"), distance=c())

sdou5 <- readTrace("output/1_validation/state_dependent_BM/state_dependent_OU_run_1.log")[[1]]
sdbm5 <- readTrace("output/1_validation/state_dependent_BM/state_dependent_BM_run_1.log")[[1]]
sdbm_maymoore5 <- readTrace("output/1_validation/state_dependent_BM/state_dependent_BM_MayMoore_run_1.log")[[1]]

kl.dist(sdou5$`sigma2s[1]`, sdbm5$`sigma2s[1]`, exp(1))
kl.dist(sdou5$`sigma2s[1]`, sdbm_maymoore5$`sigma2s[1]`, exp(1))
kl.dist(sdbm_maymoore5$`sigma2s[1]`, sdbm5$`sigma2s[1]`, exp(1))

kl.dist(sdou5$`sigma2s[2]`, sdbm5$`sigma2s[2]`, exp(1))
kl.dist(sdou5$`sigma2s[2]`, sdbm_maymoore5$`sigma2s[2]`, exp(1))
kl.dist(sdbm_maymoore5$`sigma2s[2]`, sdbm5$`sigma2s[2]`, exp(1))

# run 2
sdou6 <- readTrace("output/1_validation/state_dependent_BM/state_dependent_OU_run_2.log")[[1]]
sdbm6 <- readTrace("output/1_validation/state_dependent_BM/state_dependent_BM_run_2.log")[[1]]
sdbm_maymoore6 <- readTrace("output/1_validation/state_dependent_BM/state_dependent_BM_MayMoore_run_2.log")[[1]]

kl.dist(sdou6$`sigma2s[1]`, sdbm6$`sigma2s[1]`, exp(1))
kl.dist(sdou6$`sigma2s[1]`, sdbm_maymoore6$`sigma2s[1]`, exp(1))
kl.dist(sdbm_maymoore6$`sigma2s[1]`, sdbm6$`sigma2s[1]`, exp(1))

kl.dist(sdou6$`sigma2s[2]`, sdbm6$`sigma2s[2]`, exp(1))
kl.dist(sdou6$`sigma2s[2]`, sdbm_maymoore6$`sigma2s[2]`, exp(1))
kl.dist(sdbm_maymoore6$`sigma2s[2]`, sdbm6$`sigma2s[2]`, exp(1))

ggplot(aes(x, y, fill=)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma", direction=-1, name="Whatever")
