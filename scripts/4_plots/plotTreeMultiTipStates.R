tree <- read.tree("data/3_empirical/mammal_perMY_r500.tre")
circ <- ggtree(tree, layout = "circular")

ase <- read.simmap("output/3_empirical/ase_r500_3StateOrdered/marginal_character_run_1.tree",
                   format="phylip")
ase$maps[(ase$edge[,2] %in% 1:length(ase$tip.label))]


getlast <- function(x){
  names(tail(x, n=1))
}


missing_states <- sapply(ase$maps[(ase$edge[,2] %in% 1:length(ase$tip.label))], getlast)
two_states <- data.frame(diet=missing_states,
                 random=missing_states[sample(1:length(missing_states))])
rownames(two_states) <- tree$tip.label

p1 <- gheatmap(circ, two_states, offset=.8, width=0.05,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_manual(values=c("#44aa99", "#ddcc77", "#882255"), labels="0", "1", "2")
p1


library(ggnewscale)
p2 <- p1 + new_scale_fill()
gheatmap(p2, two_states, offset=15, width=0.05,
         colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_manual(values=c("#44aa99", "#ddcc77", "#882255"), labels="0", "1", "2")
