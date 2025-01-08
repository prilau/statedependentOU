library(RevGadgets)
library(tidyverse)
library(grid)
library(latex2exp)


# run 1 was empty
# run 3 did not converge for theta
pick_tr <- readTrace(path = c("output/3_empirical/sdOU_r500_picky/sdOU_run_2.log",
                              "output/3_empirical/sdOU_r500_picky/sdOU_run_4.log",
                              "output/3_empirical/sdOU_r500_picky/sdOU_run_5.log"), burnin = 0.1)
# run 2 did not converge for theta (but the mixing is actually better hmm)
carn_tr <- readTrace(path = c("output/3_empirical/sdOU_r500_carnivory/sdOU_run_1.log",
                              "output/3_empirical/sdOU_r500_carnivory/sdOU_run_3.log",
                              "output/3_empirical/sdOU_r500_carnivory/sdOU_run_4.log",
                              "output/3_empirical/sdOU_r500_carnivory/sdOU_run_5.log"), burnin = 0.1)
# only run 4 mix fine
pica_tr <- readTrace(path = c("output/3_empirical/sdOU_r500_pica/sdOU_run_4.log"), burnin = 0.1)


plotPosteriors <- function(trace, colors, vars, plotTitle=NA){
  names(colors) <- vars
  if (is.na(plotTitle)){
    plotTitle <- strsplit(vars[1], split="[", fixed=T)[[1]][1]
  }
  plotTrace(trace, color = colors, vars = vars)[[1]] +
    ggtitle(plotTitle)
}


pick_colors <- c('#BBCCEE', '#222255')
carn_colors <- c('#FFCCCC', '#663333')
pica_colors <- c('#C4D5CC', '#A6A6BF', '#778080','#442B44')
pick_theta    <- plotPosteriors(pick_tr, pick_colors, vars=c("theta[1]", "theta[2]"))
pick_halflife <- plotPosteriors(pick_tr, pick_colors, vars=c("halflife[1]", "halflife[2]"))
pick_stv      <- plotPosteriors(pick_tr, pick_colors, vars=c("stv[1]", "stv[2]"))

carn_theta    <- plotPosteriors(carn_tr, carn_colors, vars=c("theta[1]", "theta[2]"))
carn_halflife <- plotPosteriors(carn_tr, carn_colors, vars=c("halflife[1]", "halflife[2]"))
carn_stv      <- plotPosteriors(carn_tr, carn_colors, vars=c("stv[1]", "stv[2]"))

pica_theta    <- plotPosteriors(pica_tr, pica_colors, vars=c("theta[1]", "theta[2]", "theta[3]", "theta[4]"))
pica_halflife <- plotPosteriors(pica_tr, pica_colors, vars=c("halflife[1]", "halflife[2]", "halflife[3]", "halflife[4]"))
pica_stv      <- plotPosteriors(pica_tr, pica_colors, vars=c("stv[1]", "stv[2]", "stv[3]", "stv[4]"))


grid.newpage()
grid.draw( # draw the following matrix of plots
  rbind( # bind together the columns into a matrix
    cbind(ggplotGrob(pick_theta +
                       coord_cartesian(xlim=c(-10, 15)) +
                       theme(axis.title.y    = element_blank(),
                                           axis.title.x    = element_blank(),
                                           legend.position = "none") +
                       ggtitle(TeX("$\\theta$"))),
          ggplotGrob(pick_halflife + coord_cartesian(xlim=c(0, 3000))
                     + theme(axis.title.y    = element_blank(),
                             axis.title.x    = element_blank(),
                             legend.position = "none") +
                       ggtitle(TeX("$t_{0.5}$"))),
          ggplotGrob(pick_stv + coord_cartesian(xlim=c(0, 100))
                     + theme(axis.title.y    = element_blank(),
                                        axis.title.x    = element_blank())
                     + ggtitle(TeX("$V$"))
                     + scale_color_manual(labels = c("Not carnivores", "Carnivores"),
                                          values = c('#FFCCCC', '#663333')))),
    cbind(ggplotGrob(carn_theta +
                       coord_cartesian(xlim=c(-10, 15)) +
                       theme(axis.title.x    = element_blank(),
                             legend.position = "none",
                             plot.title = element_blank())),
          ggplotGrob(carn_halflife + coord_cartesian(xlim=c(0, 3000))
                     + theme(axis.title.x    = element_blank(),
                             axis.title.y    = element_blank(), 
                             legend.position = "none", plot.title = element_blank())),
          ggplotGrob(carn_stv + coord_cartesian(xlim=c(0, 100))
                     + theme(axis.title.x    = element_blank(),
                             axis.title.y    = element_blank(),
                             plot.title = element_blank())
          + scale_color_manual(labels = c("Not picky", "Picky"),
                               values = c('#BBCCEE', '#222255')))),
    cbind(ggplotGrob(pica_theta    + coord_cartesian(xlim=c(-10, 15)) +
                       theme(axis.title.y    = element_blank(),
                             axis.title.x    = element_blank(),
                             legend.position = "none",
                             plot.title = element_blank())),
          ggplotGrob(pica_halflife + coord_cartesian(xlim=c(0, 3000))
                     + theme(axis.title.y    = element_blank(),
                             legend.position = "none", plot.title = element_blank())),
          ggplotGrob(pica_stv + coord_cartesian(xlim=c(0, 100))
                     + theme(axis.title.x    = element_blank(),
                             axis.title.y    = element_blank(),
                             plot.title = element_blank()) +
            scale_color_manual(labels = c("Not picky, not carnivores", "Not picky, carnivores", "Picky, not carnivores", "Picky, carnivores"),
                               values = c('#C4D5CC', '#A6A6BF', '#778080','#442B44'))))
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


