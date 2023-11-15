library(RevGadgets)
library(coda)
library(ggplot2)
library(ggtree)
library(grid)
library(gridExtra)
library(tidyverse)
library(patchwork)

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




log_halfLife <- plotTrace(trace = trace_quant, vars = c("halfLife[1]","halfLife[2]","halfLife[3]","halfLife[4]", "halfLife[5]", "halfLife[6]"))[[1]]
log_sigma2 <- plotTrace(trace = trace_quant, vars = c("sigma2[1]","sigma2[2]","sigma2[3]","sigma2[4]", "sigma2[5]", "sigma2[6]"))[[1]]
theta <- plotTrace(trace = trace_quant, vars = c("theta[1]","theta[2]","theta[3]","theta[4]", "theta[5]", "theta[6]"))[[1]]

nested <- ((log_halfLife|log_sigma2)/theta)+
  plot_annotation(tag_levels = 'A') #add figure labels
ggsave("figures/gobiiformes_trace.pdf", nested, width = 400, height = 200, units = "mm")
