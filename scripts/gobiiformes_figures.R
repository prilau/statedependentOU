library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(tidyverse)
library(patchwork)
library(tikzDevice)


file <- "output/4_empirical_gobiiformes.log"

# read the trace and discard burnin
trace_quant <- readTrace(path = file, burnin = 1/7)

trace_quant[[1]][12] <- mutate(trace_quant[[1]][12], log(trace_quant[[1]][12]))
trace_quant[[1]][13] <- mutate(trace_quant[[1]][13], log(trace_quant[[1]][13]))
trace_quant[[1]][14] <- mutate(trace_quant[[1]][14], log(trace_quant[[1]][14]))
trace_quant[[1]][15] <- mutate(trace_quant[[1]][15], log(trace_quant[[1]][15]))
trace_quant[[1]][16] <- mutate(trace_quant[[1]][16], log(trace_quant[[1]][16]))
trace_quant[[1]][17] <- mutate(trace_quant[[1]][17], log(trace_quant[[1]][17]))

trace_quant[[1]][952] <- mutate(trace_quant[[1]][952], log(trace_quant[[1]][952]))
trace_quant[[1]][951] <- mutate(trace_quant[[1]][951], log(trace_quant[[1]][951]))
trace_quant[[1]][950] <- mutate(trace_quant[[1]][950], log(trace_quant[[1]][950]))
trace_quant[[1]][949] <- mutate(trace_quant[[1]][949], log(trace_quant[[1]][949]))
trace_quant[[1]][948] <- mutate(trace_quant[[1]][948], log(trace_quant[[1]][948]))
trace_quant[[1]][947] <- mutate(trace_quant[[1]][947], log(trace_quant[[1]][947]))

trace_quant[[1]][6] <- mutate(trace_quant[[1]][6], log(trace_quant[[1]][6]))
trace_quant[[1]][7] <- mutate(trace_quant[[1]][7], log(trace_quant[[1]][7]))
trace_quant[[1]][8] <- mutate(trace_quant[[1]][8], log(trace_quant[[1]][8]))
trace_quant[[1]][9] <- mutate(trace_quant[[1]][9], log(trace_quant[[1]][9]))
trace_quant[[1]][10] <- mutate(trace_quant[[1]][10], log(trace_quant[[1]][10]))
trace_quant[[1]][11] <- mutate(trace_quant[[1]][11], log(trace_quant[[1]][11]))

trace_quant[[1]][953] <- mutate(trace_quant[[1]][953], log(trace_quant[[1]][953]))
trace_quant[[1]][954] <- mutate(trace_quant[[1]][954], log(trace_quant[[1]][954]))
trace_quant[[1]][955] <- mutate(trace_quant[[1]][955], log(trace_quant[[1]][955]))
trace_quant[[1]][956] <- mutate(trace_quant[[1]][956], log(trace_quant[[1]][956]))
trace_quant[[1]][957] <- mutate(trace_quant[[1]][957], log(trace_quant[[1]][957]))
trace_quant[[1]][958] <- mutate(trace_quant[[1]][958], log(trace_quant[[1]][958]))

color_scale <- c('#6699CC', '#004488', '#EECC66', '#994455', '#997700', '#EE99AA')
names(color_scale) <- c("halfLife[1]","halfLife[2]","halfLife[3]","halfLife[4]", "halfLife[5]", "halfLife[6]")

log_halfLife <- plotTrace(trace = trace_quant,
                          vars = c("halfLife[1]", "halfLife[2]", "halfLife[3]",
                                   "halfLife[4]", "halfLife[5]", "halfLife[6]"), color = color_scale)[[1]] +
  scale_colour_manual(name = "Habitats",
                      labels = c("Brackish", "Diadromous", "Lacustrine and Riverine", "Exclusively Lacustrine", "Exclusively Marine", "Exclusively Riverine"),
                      values = color_scale) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_text(size = 4)) +
  ggtitle("$\\ln t_{1/2}$")

names(color_scale) <- c("sigma2[1]","sigma2[2]","sigma2[3]","sigma2[4]", "sigma2[5]", "sigma2[6]")
log_sigma2 <- plotTrace(trace = trace_quant,
                        vars = c("sigma2[1]","sigma2[2]","sigma2[3]","sigma2[4]", "sigma2[5]", "sigma2[6]"), color = color_scale)[[1]] +
  scale_colour_manual(name = "Habitats",
                      labels = c("Brackish", "Diadromous", "Lacustrine and Riverine", "Exclusively Lacustrine", "Exclusively Marine", "Exclusively Riverine"),
                      values = color_scale) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_text(size = 4)) +
  ggtitle("$\\ln \\sigma^2$")


names(color_scale) <- c("theta[1]","theta[2]","theta[3]",
                        "theta[4]", "theta[5]", "theta[6]")
theta <- plotTrace(trace = trace_quant,
                   vars = c("theta[1]","theta[2]","theta[3]",
                            "theta[4]", "theta[5]", "theta[6]"), color = color_scale)[[1]] +
  scale_colour_manual(name = "Habitats",
                      labels = c("Brackish", "Diadromous", "Lacustrine and Riverine", "Exclusively Lacustrine", "Exclusively Marine", "Exclusively Riverine"),
                      values = color_scale) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_text(size = 4)) +
  ggtitle("$\\theta$")


names(color_scale) <- c("alpha[1]","alpha[2]","alpha[3]",
                        "alpha[4]", "alpha[5]", "alpha[6]")
log_alpha <-  plotTrace(trace = trace_quant,
                        vars = c("alpha[1]","alpha[2]","alpha[3]",
                                 "alpha[4]", "alpha[5]", "alpha[6]"), color = color_scale)[[1]] +
  scale_color_manual(name = "Habitats",
                      labels = c("Brackish",
                                 "Diadromous",
                                 "Lacustrine and Riverine",
                                 "Exclusively Lacustrine",
                                 "Exclusively Marine",
                                 "Exclusively Riverine"),
                      values = color_scale) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_text(size = 4)) +
  ggtitle("$\\ln \\alpha$")



names(color_scale) <- c("stationaryVariance[1]",
                        "stationaryVariance[2]",
                        "stationaryVariance[3]",
                        "stationaryVariance[4]",
                        "stationaryVariance[5]",
                        "stationaryVariance[6]")
log_stV <-  plotTrace(trace = trace_quant,
                      vars = c("stationaryVariance[1]",
                               "stationaryVariance[2]",
                               "stationaryVariance[3]",
                               "stationaryVariance[4]",
                               "stationaryVariance[5]",
                               "stationaryVariance[6]"), color = color_scale)[[1]] +
  scale_colour_manual(name = "Habitats",
                      labels = c("Brackish", "Diadromous",
                                 "Lacustrine and Riverine",
                                 "Exclusively Lacustrine",
                                 "Exclusively Marine", "Exclusively Riverine"),
                      values = color_scale) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 8),
        axis.title = element_blank(),
        axis.text = element_text(size = 4)) +
  ggtitle("$V_{st}$")




nested <- ((log_halfLife|log_alpha)/(log_sigma2|log_stV)/theta)


jpeg("figures/testing.jpg", width = 5.8, height = 5, units = "in", res = 480)
nested
dev.off()

tikzDevice::tikz(file = "figures/gobiiformes.tex", width = 5.8, height = 5)
nested
dev.off()


#ggsave("figures/gobiiformes_trace.pdf", nested, width = 400, height = 80, units = "mm")
