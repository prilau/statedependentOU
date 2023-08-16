setwd("~/Desktop/ASSIM/lab_hoehna/statedependentOU/data/")
library(tidyverse)
library(ape)
library(RRphylo)
palms <- read.table("/Users/priscillalau/Desktop/ASSIM/lab_hoehna/test_data/palm/PalmTraits_1.0.txt", sep= "\t", header = TRUE)
palmTree <- read.nexus("palmTree.nex")


for (i in 1:length(palms$SpecName)){
  sp <- strsplit(palms$SpecName[i], split = " ")
  palms$SpecName[i] <- paste0(sp[[1]][1], "_", sp[[1]][2])
}
   
palmsDisc <- palms %>% 
  select(1, 6:8) %>% 
  drop_na() %>% 
  filter(Climbing != 2 & Acaulescent != 2 & Erect != 2)

palmsCont <- palms %>% 
  select(1, 16, 19) %>% 
  drop_na()

palms <- inner_join(palmsDisc, palmsCont)

omitSpecName <- c("Calamus_wailong",
                  "Iguanura_speciosa",
                  "Syagrus_campos-portoana",
                  "Syagrus_costae",
                  "Syagrus_matafome",
                  "Syagrus_tostana")

palms <- palms %>% 
  filter(!(SpecName %in% omitSpecName))

climb1 <- palms %>% filter(Climbing == 1) %>% 
  mutate(Form = 0)
acaul1 <- palms %>% filter(Acaulescent == 1) %>% 
  mutate(Form = 1)
erect1 <- palms %>% filter(Erect == 1) %>% 
  mutate(Form = 2)

climb1 <- sample_n(climb1, 100)
acaul1 <- sample_n(acaul1, 100)
erect1 <- sample_n(erect1, 100)

palmsJoined <- full_join(climb1, acaul1)
palmsJoined <- full_join(palmsJoined, erect1) %>% 
  select(1, 5:7)

palmsDisc <- palmsJoined %>% select(1, 4)

DiscTrait <- list()

for (i in 1:length(palmsDisc$SpecName)){
  DiscTrait[[palmsDisc$SpecName[i]]]$Form = palmsDisc$Form[i]
}


palmsCont <- palmsJoined %>% select(1:3)

ContTrait <- list()

for (i in 1:length(palmsCont$SpecName)){
  ContTrait[[palmsCont$SpecName[i]]]$maxBlade = palmsCont$Max_Blade_Length_m[i]
  ContTrait[[palmsCont$SpecName[i]]]$avgFruit = palmsCont$AverageFruitLength_cm[i]
  }



write.nexus.data(DiscTrait, file = "palmsDiscrete.nex", format = "standard")
write.nexus.data(ContTrait, file = "palmsContinuous.nex", format = "continuous")


hist(erect1$Max_Blade_Length_m, breaks = 50, ylim = c(0, 300), xlim = c(0, 20))
hist(climb1$Max_Blade_Length_m, add = TRUE, col = "red", breaks = 20)
hist(acaul1$Max_Blade_Length_m, add = TRUE, col = "yellow", breaks = 50)
var(palmsCont_omit$Max_Blade_Length_m)

hist(erect1$AverageFruitLength_cm, breaks = 100, ylim = c(0, 400), xlim = c(0, 20))
hist(climb1$AverageFruitLength_cm, add = TRUE, col = "red", breaks = 10)
hist(acaul1$AverageFruitLength_cm, add = TRUE, col = "yellow", breaks = 20)
var(palmsCont_omit$AverageFruitLength_cm)









is.binary(palmTree)

prunedTree <- keep.tip(palmTree, palmsJoined$SpecName)
# root age of tree
max(node.depth.edgelength(prunedTree))
#check for polytomy
is.binary(prunedTree)

#resolve polytomy
## prunedNoPolyTree <- fix.poly(prunedTree,type="resolve",node=NULL,tol=1e-10,random=TRUE)

write.tree(prunedTree, "palmPruned.tre")

