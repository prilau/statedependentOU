library(ape)
library(tidyverse)
library(RRphylo)


tree <- read.nexus("bird_tree.nex")
Traits <- read_csv("bird_traits_joined.csv")

prunedTree <- keep.tip(tree, Traits$Species)

#check for polytomy
is.binary(prunedTree)

#resolve polytomy
prunedNoPolyTree <- fix.poly(prunedTree,type="resolve",node=NULL,tol=1e-10,random=TRUE)

write.tree(prunedNoPolyTree, "bird_tree_pruned_resolvedPoly.tre")
