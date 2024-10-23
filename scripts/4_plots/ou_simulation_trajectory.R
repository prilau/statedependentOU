library(ggplot2)
library(tidyverse)

OU <- function(y, sigma2, alpha, theta, dt) {
  dy <- alpha * ( theta - y ) * dt + sigma2 * rnorm(1, mean = 0, sd = dt)
  y_new <- y + dy
  return(y_new)
}


initialState_sx <- 0.0
sim_sx=list()
for (i in 1:17){
  initialState_sx <- 0.0
  sim_sx[[i]] <- tibble(it=1:10000,
                   value=0,
                   state="0")
  for (j in 1:9999) {
    sim_sx[[i]]$value[j+1] <- OU(y= initialState_sx,
                            sigma2 = 0.5,
                            alpha = 2,
                            theta = 0.2,
                            dt = 0.001)
    initialState_sx <- sim_sx[[i]]$value[j+1]
  }
}

p1 <- ggplot(sim_sx[[1]]) +
  geom_line(data = sim_sx[[2]], aes(x=it, y=value), color="#125A56", alpha=0.6) +
  geom_line(data = sim_sx[[3]], aes(x=it, y=value), color="#00767B", alpha=0.6) +
  geom_line(data = sim_sx[[4]], aes(x=it, y=value), color="#238F9D", alpha=0.6) +
  geom_line(data = sim_sx[[5]], aes(x=it, y=value), color="#42A7C6", alpha=0.6) +
  geom_line(data = sim_sx[[6]], aes(x=it, y=value), color="#60BCE9", alpha=0.6) +
  geom_line(data = sim_sx[[7]], aes(x=it, y=value), color="#9DCCEF", alpha=0.6) +
  geom_line(data = sim_sx[[8]], aes(x=it, y=value), color="#C6DBED", alpha=0.6) +
  geom_line(data = sim_sx[[9]], aes(x=it, y=value), color="#DEE6E7", alpha=0.6) +
  geom_line(data = sim_sx[[10]], aes(x=it, y=value), color="#ECEADA", alpha=0.6) +
  geom_line(data = sim_sx[[11]], aes(x=it, y=value), color="#F0E6B2", alpha=0.6) +
  geom_line(data = sim_sx[[12]], aes(x=it, y=value), color="#F9D576", alpha=0.6) +
  geom_line(data = sim_sx[[13]], aes(x=it, y=value), color="#FFB954", alpha=0.6) +
  geom_line(data = sim_sx[[14]], aes(x=it, y=value), color="#FD9A44", alpha=0.6) +
  geom_line(data = sim_sx[[15]], aes(x=it, y=value), color="#F57634", alpha=0.6) +
  geom_line(data = sim_sx[[16]], aes(x=it, y=value), color="#E94C1F", alpha=0.6) +
  geom_line(data = sim_sx[[17]], aes(x=it, y=value), color="#D11807", alpha=0.6) +
  geom_line(aes(x=it, y=value), color="#A01813", alpha=0.6) +
  theme_classic() +
  xlab("time") +
  ylab("continuous trait") +
  scale_y_continuous(breaks = c(0, 0.2)) +
  theme(axis.text.x = element_blank())
  
initialState_sx <- 0.0
sim_sx=list()
for (i in 1:17){
  initialState_sx <- 0.0
  sim_sx[[i]] <- tibble(it=1:10000,
                        value=0,
                        state="0")
  for (j in 1:9999) {
    sim_sx[[i]]$value[j+1] <- OU(y= initialState_sx,
                                 sigma2 = 2,
                                 alpha = 0.5,
                                 theta = 0.2,
                                 dt = 0.001)
    initialState_sx <- sim_sx[[i]]$value[j+1]
  }
} 
  
p2 <- ggplot(sim_sd[[1]]) +
  geom_line(data = sim_sd[[2]], aes(x=it, y=value), color="#125A56", alpha=0.6) +
  geom_line(data = sim_sd[[3]], aes(x=it, y=value), color="#00767B", alpha=0.6) +
  geom_line(data = sim_sd[[4]], aes(x=it, y=value), color="#238F9D", alpha=0.6) +
  geom_line(data = sim_sd[[5]], aes(x=it, y=value), color="#42A7C6", alpha=0.6) +
  geom_line(data = sim_sd[[6]], aes(x=it, y=value), color="#60BCE9", alpha=0.6) +
  geom_line(data = sim_sd[[7]], aes(x=it, y=value), color="#9DCCEF", alpha=0.6) +
  geom_line(data = sim_sd[[8]], aes(x=it, y=value), color="#C6DBED", alpha=0.6) +
  geom_line(data = sim_sd[[9]], aes(x=it, y=value), color="#DEE6E7", alpha=0.6) +
  geom_line(data = sim_sd[[10]], aes(x=it, y=value), color="#ECEADA", alpha=0.6) +
  geom_line(data = sim_sd[[11]], aes(x=it, y=value), color="#F0E6B2", alpha=0.6) +
  geom_line(data = sim_sd[[12]], aes(x=it, y=value), color="#F9D576", alpha=0.6) +
  geom_line(data = sim_sd[[13]], aes(x=it, y=value), color="#FFB954", alpha=0.6) +
  geom_line(data = sim_sd[[14]], aes(x=it, y=value), color="#FD9A44", alpha=0.6) +
  geom_line(data = sim_sd[[15]], aes(x=it, y=value), color="#F57634", alpha=0.6) +
  geom_line(data = sim_sd[[16]], aes(x=it, y=value), color="#E94C1F", alpha=0.6) +
  geom_line(data = sim_sd[[17]], aes(x=it, y=value), color="#D11807", alpha=0.6) +
  geom_line(aes(x=it, y=value), color="#A01813", alpha=0.6) +
  theme_classic() +
  xlab("time") +
  ylab("continuous trait") +
  scale_y_continuous(breaks = c(0, 0.2)) +
  theme(axis.text.x = element_blank()) 
  
  
  
initialState_sx <- 0.0
sim_sx=list()
for (i in 1:17){
  initialState_sx <- 0.0
  sim_sx[[i]] <- tibble(it=1:10000,
                        value=0,
                        state="0")
  for (j in 1:9999) {
    sim_sx[[i]]$value[j+1] <- OU(y= initialState_sx,
                                 sigma2 = 10,
                                 alpha = 0.1,
                                 theta = 0.2,
                                 dt = 0.001)
    initialState_sx <- sim_sx[[i]]$value[j+1]
  }
} 
  
p3 <- ggplot(sim_sd[[1]]) +
  geom_line(data = sim_sd[[2]], aes(x=it, y=value), color="#125A56", alpha=0.6) +
  geom_line(data = sim_sd[[3]], aes(x=it, y=value), color="#00767B", alpha=0.6) +
  geom_line(data = sim_sd[[4]], aes(x=it, y=value), color="#238F9D", alpha=0.6) +
  geom_line(data = sim_sd[[5]], aes(x=it, y=value), color="#42A7C6", alpha=0.6) +
  geom_line(data = sim_sd[[6]], aes(x=it, y=value), color="#60BCE9", alpha=0.6) +
  geom_line(data = sim_sd[[7]], aes(x=it, y=value), color="#9DCCEF", alpha=0.6) +
  geom_line(data = sim_sd[[8]], aes(x=it, y=value), color="#C6DBED", alpha=0.6) +
  geom_line(data = sim_sd[[9]], aes(x=it, y=value), color="#DEE6E7", alpha=0.6) +
  geom_line(data = sim_sd[[10]], aes(x=it, y=value), color="#ECEADA", alpha=0.6) +
  geom_line(data = sim_sd[[11]], aes(x=it, y=value), color="#F0E6B2", alpha=0.6) +
  geom_line(data = sim_sd[[12]], aes(x=it, y=value), color="#F9D576", alpha=0.6) +
  geom_line(data = sim_sd[[13]], aes(x=it, y=value), color="#FFB954", alpha=0.6) +
  geom_line(data = sim_sd[[14]], aes(x=it, y=value), color="#FD9A44", alpha=0.6) +
  geom_line(data = sim_sd[[15]], aes(x=it, y=value), color="#F57634", alpha=0.6) +
  geom_line(data = sim_sd[[16]], aes(x=it, y=value), color="#E94C1F", alpha=0.6) +
  geom_line(data = sim_sd[[17]], aes(x=it, y=value), color="#D11807", alpha=0.6) +
  geom_line(aes(x=it, y=value), color="#A01813", alpha=0.6) +
  theme_classic() +
  xlab("time") +
  ylab("continuous trait") +
  scale_y_continuous(breaks = c(0, 0.2)) +
  theme(axis.text.x = element_blank()) 
p3
  
  
  
  
  
initialState_sd <- 0.0
sim_sd=list()
for (i in 1:17){
  initialState_sd <- 0.0
  sim_sd[[i]] <- tibble(it=1:10000,
                        value=0,
                        state="0")
  for (j in 1:2999) {
    sim_sd[[i]]$value[j+1] <- OU(y= initialState_sd,
                                 sigma2 = 0.5,
                                 alpha = 0.5,
                                 theta = -0.2,
                                 dt = 0.001)
    sim_sd[[i]]$state[j+1] = "1"
    initialState_sd <- sim_sd[[i]]$value[j+1]
  }
  for (j in 3000:9999) {
    sim_sd[[i]]$value[j+1] <- OU(y= initialState_sd,
                                 sigma2 = 1.5,
                                 alpha = 1.0,
                                 theta = 0.2,
                                 dt = 0.001)
    initialState_sd <- sim_sd[[i]]$value[j+1]
  }
} 


p4 <- ggplot(sim_sd[[1]][1:3000,]) +
  geom_line(aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[2]][1:3000,], aes(x=it, y=value), color = "#882255", alpha=0.8) +
  geom_line(data = sim_sd[[3]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[4]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[5]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[6]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[7]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[8]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[9]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[10]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[11]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[12]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[13]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[14]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[15]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[16]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[17]][1:3000,], aes(x=it, y=value), color="#882255", alpha=0.8) +
  geom_line(data = sim_sd[[1]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[2]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[3]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[4]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[5]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[6]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[7]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[8]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[9]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[10]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[11]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[12]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[13]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[14]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[15]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[16]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  geom_line(data = sim_sd[[17]][3001:10000,], aes(x=it, y=value), color="#44aa99", alpha=0.8) +
  theme_classic() +
  xlab("time") +
  ylab("trait") +
  #scale_y_continuous(breaks = c(-0.1, 0.1)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) 

ggsave("sxou.pdf", p1, width = 100, height = 60, units = "mm")
ggsave("sxou2.pdf", p2, width = 100, height = 60, units = "mm")
ggsave("sxou3.pdf", p3, width = 100, height = 60, units = "mm")
ggsave("sdou.png", p4, width = 160, height = 90, units = "mm")
