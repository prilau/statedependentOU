library(ggplot2)
library(tidyverse)
library(patchwork)
library(tikzDevice)


num_tips = c('100', '250', '500')
parameters = c("$t_{1/2}$", "$\\alpha$", "$V_{st}$", "$\\sigma^2$", "$\\theta_0$", "$\\theta_1$")
shape = c("unimodal", "multimodal")
count = 0

grid = expand.grid(shape=shape, parameters=parameters, num_tips=num_tips, count = count,
                   stringsAsFactors=FALSE)
  

grid$count[1] = 1
grid$count[2] = 15
grid$count[3] = 12
grid$count[4] = 4
grid$count[5] = 1
grid$count[6] = 15
grid$count[7] = 16
grid$count[8] = 0
grid$count[9] = 8
grid$count[10] = 8
grid$count[11] = 1
grid$count[12] = 15

grid$count[13] = 11
grid$count[14] = 5
grid$count[15] = 6
grid$count[16] = 10
grid$count[17] = 11
grid$count[18] = 5
grid$count[19] = 16
grid$count[20] = 0
grid$count[21] = 15
grid$count[22] = 1
grid$count[23] = 11
grid$count[24] = 5

grid$count[25] = 15
grid$count[26] = 1
grid$count[27] = 11
grid$count[28] = 5
grid$count[29] = 15
grid$count[30] = 1
grid$count[31] = 16
grid$count[32] = 0
grid$count[33] = 16
grid$count[34] = 0
grid$count[35] = 14
grid$count[36] = 2

grid <- grid %>% mutate(count = count/16)
grid <- grid %>% rename(frequency = count)


shape_color <- c('#004488', '#DDAA33')
names(shape_color) <- unique(grid$shape)
level_order <- c('100', '250', '500') 





sim2_shape <- ggplot(grid, aes(x = factor(num_tips, levels = level_order),
                               y = frequency, fill = shape)) +
  geom_bar(position = position_stack(reverse = FALSE), stat='identity') +
  theme(legend.position = "right") +
  scale_fill_manual(name = "Shape of distribution",values = shape_color) +
  facet_wrap(facets = vars(parameters)) +
  xlab("Number of tree tips") +
  ylab("Frequency") +
  theme(legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      #plot.title = element_text(size = 8),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 4))  +
  theme_bw()


jpeg("figures/testing.jpg", width = 5.8, height = 5, units = "in", res = 480)
sim2_shape
dev.off()

tikzDevice::tikz(file = "figures/sim2_shape.tex", width = 5.8, height = 5)
sim2_shape
dev.off()





rates = c('5', '10', '20', '50', '500')
parameters = c("$t_{1/2}$", "$\\alpha$", "$V_{st}$", "$\\sigma^2$", "$\\theta_0$", "$\\theta_1$")
shape = c("unimodal", "multimodal")
count = 0

grid = expand.grid(shape=shape, parameters=parameters, rates=rates, count = count,
                   stringsAsFactors=FALSE)
  

grid$count[1] = 15
grid$count[2] = 1
grid$count[3] = 14 
grid$count[4] = 2
grid$count[5] = 14
grid$count[6] = 2
grid$count[7] = 15
grid$count[8] = 1
grid$count[9] = 15
grid$count[10] = 1
grid$count[11] = 15
grid$count[12] = 1

grid$count[13] = 16
grid$count[14] = 0
grid$count[15] = 11
grid$count[16] = 5
grid$count[17] = 16
grid$count[18] = 0
grid$count[19] = 14
grid$count[20] = 2
grid$count[21] = 16
grid$count[22] = 0
grid$count[23] = 16
grid$count[24] = 0

grid$count[25] = 13
grid$count[26] = 3
grid$count[27] = 12
grid$count[28] = 6
grid$count[29] = 14
grid$count[30] = 2
grid$count[31] = 14
grid$count[32] = 2
grid$count[33] = 16
grid$count[34] = 0
grid$count[35] = 15
grid$count[36] = 1

grid$count[37] = 14
grid$count[38] = 2
grid$count[39] = 2
grid$count[40] = 14
grid$count[41] = 14
grid$count[42] = 2
grid$count[43] = 16
grid$count[44] = 0
grid$count[45] = 16
grid$count[46] = 0
grid$count[47] = 14
grid$count[48] = 2

grid$count[49] = 0
grid$count[50] = 16
grid$count[51] = 3
grid$count[52] = 13
grid$count[53] = 0
grid$count[54] = 16
grid$count[55] = 0
grid$count[56] = 16
grid$count[57] = 2
grid$count[58] = 14
grid$count[59] = 1
grid$count[60] = 15





grid <- grid %>% mutate(count = count/16)
grid <- grid %>% rename(frequency = count)


shape_color <- c('#004488', '#DDAA33')
names(shape_color) <- unique(grid$shape)
level_order <- c('100', '250', '500') 





sim2_shape <- ggplot(grid, aes(x = factor(num_tips, levels = level_order),
                               y = frequency, fill = shape)) +
  geom_bar(position = position_stack(reverse = FALSE), stat='identity') +
  theme(legend.position = "right") +
  scale_fill_manual(name = "Shape of distribution",values = shape_color) +
  facet_wrap(facets = vars(parameters)) +
  xlab("Number of tree tips") +
  ylab("Frequency") +
  theme(legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      #plot.title = element_text(size = 8),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 4))  +
  theme_bw()


jpeg("figures/testing.jpg", width = 5.8, height = 5, units = "in", res = 480)
sim2_shape
dev.off()

tikzDevice::tikz(file = "figures/sim2_shape.tex", width = 5.8, height = 5)
sim2_shape
dev.off()

#ggsave("figures/sim2_shape.pdf", sim2_shape,
#       width = 400, height = 300, units = "mm")

