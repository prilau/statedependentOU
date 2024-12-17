library(phytools)
library(viridis)

# !!! In rb stochastic mapping, if discrete character index != 1, simmap output will make tip states = character[1]

t <- read.tree("testing/mammal_perMY_n2955.tre")
disc <- read.nexus.data("testing/mammal_diet_binary_Discrete.nex")
a <- read.simmap("testing/dummy2_marginal_character.tree", format="phylip")

ntip <- length(disc)

sp <- names(disc)
disc1 <- disc2 <- c()
for (i in 1:length(disc)){
  disc1[i] <- disc[[i]][1]
  disc2[i] <- disc[[i]][2]
}
disc1 <- unlist(disc1)
disc2 <- unlist(disc2)
names(disc1) <- names(disc2) <- sp
sims <- make.simmap(tree = t, x = disc2, nsim = 1)
sims_alt <- make.simmap(tree = t, x = disc1, nsim = 1)
#obj<-densityMap(sims,states=c("0","1"),plot=FALSE)
#obj<-setMap(obj,c("#0077bb", "#ee7733"))


#subset simmaps
a2 <- a
sims2 <- sims
sims_alt2 <- sims_alt
#a2$tip.label <- paste0("t", 1:ntip)
#sims2$tip.label <- paste0("t", 1:ntip)
a2 <- keep.tip(a2, tip=paste0("t", 1:100))
sims2 <- keep.tip(sims2, tip=paste0("t", 1:100))
sims_alt2 <- keep.tip(sims_alt2, tip=paste0("t", 1:100))

par(mfrow=c(2,2))
state_colors <- c("0"="#ee7733", "1"="#0077bb")
plot(a2, color=state_colors, lwd = 1)
plot(sims2, color=state_colors, lwd = 1, direction="leftwards")
plot(a2, color=state_colors, lwd = 1)
plot(sims_alt2, color=state_colors, lwd = 1, direction="leftwards")


# find tip branches
tip_branch_index <- which(a$edge[,2] <= ntip)
names(tip_branch_index) <- a$tip.label
rb_end_states <- c()
for(i in 1:length(tip_branch_index)){
  tip_seg <- length(a$maps[[tip_branch_index[i]]])
  tip_end_state <- names(a$maps[[tip_branch_index[i]]][tip_seg])
  names(tip_end_state) <- names(tip_branch_index[i])
  rb_end_states <- append(rb_end_states, tip_end_state)
}

tip_branch_index <- which(sims$edge[,2] <= ntip)
names(tip_branch_index) <- sims$tip.label
pt_end_states <- c()
for(i in 1:length(tip_branch_index)){
  tip_seg <- length(sims$maps[[tip_branch_index[i]]])
  tip_end_state <- names(sims$maps[[tip_branch_index[i]]][tip_seg])
  names(tip_end_state) <- names(tip_branch_index[i])
  pt_end_states <- append(pt_end_states, tip_end_state)
}

# tip states in simulation
disc_states <- rep(0, ntip)
names(disc_states) <- a$tip.label
for (i in names(disc_states)){
  disc_states[which(names(disc_states) == i)] <- disc[[i]]
}

compare_rb_pt <- rep(NA, ntip)
names(compare_rb_pt) <- a$tip.label
for (i in names(compare_rb_pt)){
  compare_rb_pt[which(names(compare_rb_pt) == i)] <- ifelse(rb_end_states[which(names(rb_end_states) == i)] == pt_end_states[which(names(pt_end_states) == i)],
                                                              T, F)
}

compare_rb_nex <- rep(NA, ntip)
names(compare_rb_nex) <- a$tip.label
for (i in names(compare_rb_nex)){
  compare_rb_nex[which(names(compare_rb_nex) == i)] <- ifelse(rb_end_states[which(names(rb_end_states) == i)] == disc_states[which(names(disc_states) == i)],
                                                            T, F)
}

compare_pt_nex <- rep(NA, ntip)
names(compare_pt_nex) <- a$tip.label
for (i in names(compare_pt_nex)){
  compare_pt_nex[which(names(compare_pt_nex) == i)] <- ifelse(pt_end_states[which(names(pt_end_states) == i)] == disc_states[which(names(disc_states) == i)],
                                                              T, F)
}




# simulated trees
library(TESS)
ntip = 2955
t = ladderize(tess.sim.taxa(1, ntip, 10, 1, 0.5)[[1]])
#t$edge.length = t$edge.length / max(branching.times(t))

#disc <- read.nexus.data("testing/dummy_Discrete.nex")
disc <- list()
for (i in 1:length(t$tip.label)){
  sp = t$tip.label[i]
  disc[[sp]] <- ifelse(as.numeric(gsub("t", "", t$tip.label[i])) <= ntip/2, 0, 1)
}

write.tree(t,file = "testing/dummy.tre")
write.nexus.data(disc, file = "testing/dummy_Discrete.nex", format = "standard")

par(mfrow=c(1,2))
tip_colors = ifelse(as.numeric(gsub("t", "", t$tip.label)) <= ntip/2, "#ee7733", "#0077bb")
plot(t,tip.color=tip_colors)

a <- read.simmap("testing/dummy_marginal_character.tree", format="phylip")
state_colors <- c("0"="#ee7733", "1"="#0077bb")
plot(a, color=state_colors, lwd = 1, direction="leftwards")


tip_branch_index <- which(a$edge[,2] <= ntip)
names(tip_branch_index) <- a$tip.label
rb_end_states <- c()
for(i in 1:length(tip_branch_index)){
  tip_seg <- length(a$maps[[tip_branch_index[i]]])
  tip_end_state <- names(a$maps[[tip_branch_index[i]]][tip_seg])
  names(tip_end_state) <- names(tip_branch_index[i])
  rb_end_states <- append(rb_end_states, tip_end_state)
}
disc_states <- rep(0, ntip)
names(disc_states) <- a$tip.label
for (i in names(disc_states)){
  disc_states[which(names(disc_states) == i)] <- disc[[i]]
}
compare_rb_nex <- rep(NA, ntip)
names(compare_rb_nex) <- a$tip.label
for (i in names(compare_rb_nex)){
  compare_rb_nex[which(names(compare_rb_nex) == i)] <- ifelse(rb_end_states[which(names(rb_end_states) == i)] == disc_states[which(names(disc_states) == i)],
                                                              T, F)
}


#post_colors <- magma(51, alpha = 1, begin = 0, end = 1, direction = 1)
#names(post_colors) <- as.character(100:50)
#plot(a, colors=post_colors)

# change tip labels
t <- read.tree("data/3_empirical/mammal_perMY_n2955.tre")
disc <- read.nexus.data("data/3_empirical/mammal_diet_binary_Discrete.nex")

sps <- names(disc)
names(disc) <- paste0("t", 1:ntip)
for (i in 1:ntip){
  sp <- sps[i]
  t$tip.label[which(t$tip.label == sp)] <- names(disc)[i]
}

write.tree(t,file = "testing/mammal_perMY_n2955.tre")
write.nexus.data(disc, file = "testing/mammal_diet_binary_Discrete.nex", format = "standard")
