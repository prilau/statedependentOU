library(ape)
library(phytools)
library(RevGadgets)
library(tidyverse)
source("scripts/5_miscellaneous/functions.R")

tree <- read.tree("data/3_empirical/mammal_perMY_r500.tre")
dir_in = "output/3_empirical/sdOU_r500_missingStateModel/"
files <- list.files(dir_in, recursive = TRUE)[grepl(".trees", list.files(dir_in, recursive = TRUE))]

log_to_simmap(dir_in, files)

for (i in 1:length(files)){
  file_in <- paste0(dir_in, files[i])
  file_out <- paste0(gsub(pattern="trees", replacement="", x=file_in), "log")
  simmap_to_ancStates(file_in, file_out, tree)
}
