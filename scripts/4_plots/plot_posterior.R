library(RevGadgets)
library(tidyverse)
library(grid)

greenPrior <- readTrace(path = "output/3_empirical/3g_statedependentOU_greenPrior/3g.log", burnin = 0.1)
blackPrior <- readTrace(path = "output/3_empirical/3h_statedependentOU_blackPrior/3h.log", burnin = 0.1)

changenames <- colnames(greenPrior[[1]])
names(changenames) <- paste0("g_", colnames(greenPrior[[1]]))
greenPrior[[1]] <- greenPrior[[1]] %>% rename(all_of(changenames))

changenames <- colnames(blackPrior[[1]])
names(changenames) <- paste0("h_", colnames(blackPrior[[1]]))
blackPrior[[1]] <- blackPrior[[1]] %>% rename(all_of(changenames))

comparePriors <- list()
comparePriors[[1]] <- cbind(greenPrior[[1]], blackPrior[[1]])


plotPosteriors <- function(trace, colors, vars, plotTitle=NA){
  names(colors) <- vars
  if (is.na(plotTitle)){
    plotTitle <- strsplit(vars[1], split="[", fixed=T)[[1]][1]
  }
  plotTrace(trace, color = colors, vars = vars)[[1]] +
    ggtitle(plotTitle)
}


colors <- c("#44aa99", "#ddcc77", "#882255")
h_alpha <- plotPosteriors(blackPrior, colors, vars=c("h_alpha[1]", "h_alpha[2]", "h_alpha[3]"))
h_sigma2 <- plotPosteriors(blackPrior, colors, vars=c("h_sigma2[1]", "h_sigma2[2]", "h_sigma2[3]"))
h_theta <- plotPosteriors(blackPrior, colors, vars=c("h_theta[1]", "h_theta[2]", "h_theta[3]"))
h_halflife <- plotPosteriors(blackPrior, colors, vars=c("h_halflife[1]", "h_halflife[2]", "h_halflife[3]"))
h_stv <- plotPosteriors(blackPrior, colors, vars=c("h_stv[1]", "h_stv[2]", "h_stv[3]"))
h_rho <- plotPosteriors(blackPrior, colors, vars=c("h_rho[1]", "h_rho[2]", "h_rho[3]"))

g_alpha <- plotPosteriors(greenPrior, colors, vars=c("g_alpha[1]", "g_alpha[2]", "g_alpha[3]"))
g_sigma2 <- plotPosteriors(greenPrior, colors, vars=c("g_sigma2[1]", "g_sigma2[2]", "g_sigma2[3]"))
g_theta <- plotPosteriors(greenPrior, colors, vars=c("g_theta[1]", "g_theta[2]", "g_theta[3]"))
g_halflife <- plotPosteriors(greenPrior, colors, vars=c("g_halflife[1]", "g_halflife[2]", "g_halflife[3]"))
g_stv <- plotPosteriors(greenPrior, colors, vars=c("g_stv[1]", "g_stv[2]", "g_stv[3]"))
g_rho <- plotPosteriors(greenPrior, colors, vars=c("g_rho[1]", "g_rho[2]", "g_rho[3]"))



grid.newpage()
grid.draw( # draw the following matrix of plots
  rbind( # bind together the columns into a matrix
    cbind(ggplotGrob(h_alpha    + theme(axis.title.x    = element_blank(),
                                        legend.position = "none")),
          ggplotGrob(h_sigma2   + theme(axis.title.y    = element_blank(),
                                        axis.title.x    = element_blank(),
                                        legend.position = "none")),
          ggplotGrob(h_theta    + theme(axis.title.y    = element_blank(),
                                        axis.title.x    = element_blank())
                     + scale_color_manual(labels = c("Herbivores", "Omnivores", "Carnivores"),
                                          values = c("#44aa99", "#ddcc77", "#882255")))),
    cbind(ggplotGrob(h_halflife + theme(legend.position = "none")),
          ggplotGrob(h_stv      + theme(axis.title.y    = element_blank(),
                                        legend.position = "none")),
          ggplotGrob(h_rho      + theme(axis.title.y    = element_blank(),
                                        legend.position = "none")))
    )
)

grid.newpage()
grid.draw( # draw the following matrix of plots
  rbind( # bind together the columns into a matrix
    cbind(ggplotGrob(g_alpha    + theme(axis.title.x    = element_blank(),
                                        legend.position = "none")),
          ggplotGrob(g_sigma2   + theme(axis.title.y    = element_blank(),
                                        axis.title.x    = element_blank(),
                                        legend.position = "none")),
          ggplotGrob(g_theta    + theme(axis.title.y    = element_blank(),
                                        axis.title.x    = element_blank())
                     + scale_color_manual(labels = c("Herbivores", "Omnivores", "Carnivores"),
                                          values = c("#44aa99", "#ddcc77", "#882255")))),
    cbind(ggplotGrob(g_halflife + theme(legend.position = "none")),
          ggplotGrob(g_stv      + theme(axis.title.y    = element_blank(),
                                        legend.position = "none")),
          ggplotGrob(g_rho      + theme(axis.title.y    = element_blank(),
                                        legend.position = "none")))
  )
)








#log_norm_err <- function(mu){
#  sigma = 1.044495
#  tree_height = 200
#  quants <- c(qlnorm(0.025, meanlog=mu, sdlog=sigma), qlnorm(0.975, meanlog=mu, sdlog=sigma))
#  sq_err <- (quants - c(0.05*tree_height, 3*tree_height))^2
#  sum_sq_err <- sum(sq_err)
#  return(sum_sq_err)
#}
#
#
#log_norm_err(1, 1)
#
#
#
#optim(par=c(1), fn=log_norm_err, method = "L-BFGS-B")
#
#
#c(qlnorm(0.05, meanlog=0, sdlog=0.674), qlnorm(0.95, meanlog=0, sdlog=0.674))
#
#
#tree <- drop.tip(tree, tip = c("Galictis_cuja", "Hemibelideus_lemuroides"))
#write.tree(tree, "Desktop/statedependentOU/data/3_empirical/mammal_diet_perMY.tre")
#
#tree <- read.tree("Desktop/statedependentOU/data/3_empirical/mammal_diet.tre")


