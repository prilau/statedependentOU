library(ape)
library(phytools)
library(RevGadgets)
library(tidyverse)
source("scripts/5_miscellaneous/functions.R")

tree <- read.tree("data/1_validation/artiodactyla/artiodactyla.tree")
dir_in = "output/2_simulation/missing_state/"
files <- list.files(dir_in, recursive = TRUE)[grepl(".trees", list.files(dir_in, recursive = TRUE))]

log_to_simmap(files)

for (i in 1:length(files)){
  file_in <- paste0(dir_in, files[i])
  file_out <- paste0(gsub(pattern="trees", replacement="", x=file_in), "log")
  simmap_to_ancStates(file_in, file_out, tree)
}




# 01.26.2024 redo the rb step! 

ase <- processAncStates("output/2_simulation/missing_state/sim_1/anc_states_nstate_2_run_1.log",
                        state_labels=c("0"="0",
                                       "1"="1",
                                       "2"="2",
                                       "3"="3"))




aug <- processAncStates("output/3_empirical/aug_r500_3StateModel/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivorous",
                                       "1"="Omnivorous",
                                       "2"="Carnivorous"))


ggplot() +
  geom_point(aes(x = ase@data$anc_state_1_pp, y = aug@data$anc_state_1_pp))

