library(ggplot2)
library(tidyverse)
library(patchwork)


num_tips = c('100', '250', '500')
parameters = c("halfLife", "alpha", "stV", "sigma2", "theta0", "theta1")
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
  scale_fill_manual(name = "shape",values = shape_color) +
  facet_wrap(facets = vars(parameters)) +
  theme_bw()
ggsave("figures/sim2_shape.pdf", sim2_shape,
       width = 400, height = 300, units = "mm")

