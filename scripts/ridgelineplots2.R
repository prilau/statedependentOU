library(RevGadgets)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(patchwork)
library(grid)
library(tikzDevice)
library(cowplot)


cont_trait_labels <- c(
                       "Adductor mass",
                       "Ascending process",
                       "Raker length",
                       "Buccal length",
                       "Buccal width"
                       )


these_paths <- c(
     paste0("output/4_empirical_haemulidae_", 2, ".log"),
     paste0("output/4_empirical_haemulidae_", 5, ".log"),
     paste0("output/4_empirical_haemulidae_", 6, ".log"),
     paste0("output/4_empirical_haemulidae_", 9, ".log"),
     paste0("output/4_empirical_haemulidae_", 10, ".log")
     )


trace_quant <- readTrace(path = these_paths, burnin = 1/7)
halfLife <- list()
alpha <- list()
V_st <- list()
sigma2 <- list()
theta <- list()

for (i in 1:5) {
  subtitle <- cont_trait_labels[i]
  halfLife[[i]] <- plotTrace(trace = trace_quant,
                             vars = c("halfLife[1]","halfLife[2]"),
                             color = c("halfLife[1]" = "#EE7733",
                                       "halfLife[2]" = "#0077BB"))[[i]] +
    xlab("") +
    ylab("") +
    scale_fill_discrete(name = "$t_{1/2}$", labels = c("Reef-dwelling", "Non-reef-dwelling")) +
    ggtitle("", subtitle = subtitle) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 5),
      axis.text = element_text(size = 3)
    )
  if (i == 1){
    halfLife[[i]] <- halfLife[[i]] +
      theme(plot.title = element_text(size = 7)) +
      ggtitle("$t_{1/2}$")
  }
  if (i == 3) {
    halfLife[[i]] <- halfLife[[i]] +
      ylab("Density") +
      theme(axis.title = element_text(size = 5))
  }
}

for (i in 1:5) {
  subtitle <- cont_trait_labels[i]
  alpha[[i]] <- plotTrace(trace = trace_quant,
                             vars = c("alpha[1]","alpha[2]"),
                             color = c("alpha[1]" = "#EE7733",
                                       "alpha[2]" = "#0077BB"))[[i]] +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 5),
      axis.title = element_blank(),
      axis.text = element_text(size = 3)
    )
  V_st[[i]] <- plotTrace(trace = trace_quant,
                             vars = c("stationaryVariance[1]","stationaryVariance[2]"),
                             color = c("stationaryVariance[1]" = "#EE7733",
                                       "stationaryVariance[2]" = "#0077BB"))[[i]] +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 5),
      axis.title = element_blank(),
      axis.text = element_text(size = 3)
    )
  sigma2[[i]] <- plotTrace(trace = trace_quant,
                             vars = c("sigma2[1]","sigma2[2]"),
                             color = c("sigma2[1]" = "#EE7733",
                                       "sigma2[2]" = "#0077BB"))[[i]] +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 5),
      axis.title = element_blank(),
      axis.text = element_text(size = 3)
    )
  theta[[i]] <- plotTrace(trace = trace_quant,
                             vars = c("theta[1]","theta[2]"),
                             color = c("theta[1]" = "#EE7733",
                                       "theta[2]" = "#0077BB"))[[i]] +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 5),
      axis.title = element_blank(),
      axis.text = element_text(size = 3)
    )
  if (i == 1){

    alpha[[i]] <- alpha[[i]] +
      theme(plot.title = element_text(size = 7)) +
      ggtitle("$\\alpha$")
    V_st[[i]] <- V_st[[i]] +
      theme(plot.title = element_text(size = 7)) +
      ggtitle("$V_{st}$")
    sigma2[[i]] <- sigma2[[i]] +
      theme(plot.title = element_text(size = 7)) +
      ggtitle("$\\sigma^2$")
    theta[[i]] <- theta[[i]] +
      theme(plot.title = element_text(size = 7)) +
      ggtitle("$\\theta$")
  }
  if (i == 5) {
    V_st[[i]] <- V_st[[i]] +
      xlab("Value") +
      ylab("") +
      theme(axis.title = element_text(size = 5))
  }
}


halfLife_plot <-
  rbind(
    ggplotGrob(halfLife[[1]]),
    ggplotGrob(halfLife[[2]]),
    ggplotGrob(halfLife[[3]]),
    ggplotGrob(halfLife[[4]]),
    ggplotGrob(halfLife[[5]])
)

alpha_plot <-
  rbind(
    ggplotGrob(alpha[[1]]),
    ggplotGrob(alpha[[2]]),
    ggplotGrob(alpha[[3]]),
    ggplotGrob(alpha[[4]]),
    ggplotGrob(alpha[[5]])
  )

V_st_plot <-
  rbind(
    ggplotGrob(V_st[[1]]),
    ggplotGrob(V_st[[2]]),
    ggplotGrob(V_st[[3]]),
    ggplotGrob(V_st[[4]]),
    ggplotGrob(V_st[[5]])
  )

sigma2_plot <-
  rbind(
    ggplotGrob(sigma2[[1]]),
    ggplotGrob(sigma2[[2]]),
    ggplotGrob(sigma2[[3]]),
    ggplotGrob(sigma2[[4]]),
    ggplotGrob(sigma2[[5]])
  )

theta_plot <-
  rbind(
    ggplotGrob(theta[[1]]),
    ggplotGrob(theta[[2]]),
    ggplotGrob(theta[[3]]),
    ggplotGrob(theta[[4]]),
    ggplotGrob(theta[[5]])
  )


nested <- cbind(
  halfLife_plot,
  alpha_plot,
  V_st_plot,
  sigma2_plot,
  theta_plot
)





jpeg("figures/testing.jpg", width = 5.8, height = 5.8, units = "in", res = 480)
grid.newpage()
grid.draw(
  nested
)
dev.off()



tikzDevice::tikz(file = "figures/haemulidae.tex", width = 5.8, height = 5.8)
grid.newpage()
grid.draw(
  nested
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
