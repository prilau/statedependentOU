library(tidyverse)
library(ape)
library(RRphylo)

fish <- read.csv("old/fish_DiscreteContinuous.csv")

gobii <- fish %>% 
  filter(Rab18tax.Order == "Gobiiformes") %>% 
  select(1, 9:11)

DiscTrait <- list()
ContTrait <- list()

for (i in 1:length(gobii$sciname_)){
  DiscTrait[[gobii$sciname_[i]]] <- list(gobii$Final.hab.info[i])
}

for (i in 1:length(gobii$sciname_)){
  ContTrait[[gobii$sciname_[i]]] <- list(gobii$nlog.TL[i])
}

write.nexus.data(DiscTrait, file = "data/gobiiformesDiscrete.nex", format = "Standard")
write.nexus.data(ContTrait, file = "data/gobiiformesContinuous.nex", format = "Continuous")

#gobii %>% filter(Final.hab.info == 6) %>% 
#  count()

fishTree <- read.tree(file = "old/Fish_Chang2019.tre")

prunedTree <- keep.tip(fishTree, gobii$sciname_)
is.binary(prunedTree)
write.tree(prunedTree, "data/gobiiformes.tre")


var(gobii$nlog.TL)
