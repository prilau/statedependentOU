library(RevGadgets)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(patchwork)
library(grid)

cont_trait_labels <- c("body size",
                       "adductor mass",
                       "ascending process",
                       "raker length",
                       "eye width",
                       "buccal length",
                       "buccal width",
                       "head height",
                       "head length")


these_paths <- c(paste0("output/4_empirical_haemulidae_", 1, ".log"),
     paste0("output/4_empirical_haemulidae_", 2, ".log"),
     paste0("output/4_empirical_haemulidae_", 5, ".log"),
     paste0("output/4_empirical_haemulidae_", 6, ".log"),
     paste0("output/4_empirical_haemulidae_", 7, ".log"),
     paste0("output/4_empirical_haemulidae_", 9, ".log"),
     paste0("output/4_empirical_haemulidae_", 10, ".log"),
     paste0("output/4_empirical_haemulidae_", 11, ".log"),
     paste0("output/4_empirical_haemulidae_", 12, ".log"))


trace_quant <- readTrace(path = these_paths, burnin = 1/7)
halfLife <- list()
for (i in 1:9) {
  title <- cont_trait_labels[i]
  halfLife[[i]] <- plotTrace(trace = trace_quant,
                             vars = c("halfLife[1]","halfLife[2]"),
                             color = c("halfLife[1]" = "#EE7733",
                                       "halfLife[2]" = "#0077BB"))[[i]] +
    ggtitle(title)
}


pdf("figures/haemulidae_halfLife.pdf", height = 11, width = 8.5, paper = "letter")
grid.newpage()
grid.draw(
  rbind(
    ggplotGrob(halfLife[[1]]),
    ggplotGrob(halfLife[[2]]),
    ggplotGrob(halfLife[[3]]),
    ggplotGrob(halfLife[[4]]),
    ggplotGrob(halfLife[[5]]),
    ggplotGrob(halfLife[[6]]),
    ggplotGrob(halfLife[[7]]),
    ggplotGrob(halfLife[[8]]),
    ggplotGrob(halfLife[[9]])
  )
)
dev.off()



alpha <- plotTrace(trace = trace_quant,
                      vars = c("alpha[1]","alpha[2]"),
                      color = c("alpha[1]" = "#EE7733",
                                "alpha[2]" = "#0077BB"))[[1]] +
  ggtitle("Alpha")

stV <- plotTrace(trace = trace_quant,
                      vars = c("stationaryVariance[1]","stationaryVariance[2]"),
                      color = c("stationaryVariance[1]" = "#EE7733",
                                "stationaryVariance[2]" = "#0077BB"))[[1]] +
  ggtitle("Stationary Variance")

sigma2 <- plotTrace(trace = trace_quant,
                      vars = c("sigma2[1]","sigma2[2]"),
                      color = c("sigma2[1]" = "#EE7733",
                                "sigma2[2]" = "#0077BB"))[[1]] +
  ggtitle("sigma2")

theta <- plotTrace(trace = trace_quant,
                      vars = c("theta[1]","theta[2]"),
                      color = c("theta[1]" = "#EE7733",
                                "theta[2]" = "#0077BB"))[[1]] +
  ggtitle("theta")





#cuteRidgePlot_halfLife <- looongPosteriorDf %>% 
#  filter(parameter %in% c("halfLife.1.", "halfLife.2.")) %>% 
#  group_by(trait) %>% 
#  ggplot(aes(x = value, y = trait, fill = parameter, alpha = 0.7)) +
#  geom_density_ridges() +
#  xlim(0, 3) +
#  theme_ridges() + 
#  theme(legend.position = "none",
#        axis.title.x=element_blank(),
#        axis.title.y=element_blank(),
#        plot.title = element_text(hjust = 0.5)) +
#  ggtitle("half life")


nested <- ((halfLife|alpha)/(stV|sigma2)/theta)


ggsave("figures/haemulidae_all_parameters.pdf", nested, width = 400, height = 200, units = "mm")
