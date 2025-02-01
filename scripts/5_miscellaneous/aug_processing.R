library(ape)
source("scripts/5_miscellaneous/functions.R")

tree <- read.tree("data/1_validation/artiodactyla/artiodactyla.tree")
dir_in = "output/2_simulation/missing_state/"
files <- list.files(dir_in, recursive = TRUE)[grepl(".trees", list.files(dir_in, recursive = TRUE))]

#log_to_simmap(files)

for (i in 1:length(files)){
  file_in <- paste0(dir_in, files[i])
  file_out <- paste0(gsub(pattern="trees", replacement="", x=file_in), "log")
  simmap_to_ancStates(file_in, file_out, tree)
}

for (i in 1:length(files)){
  file_in <- paste0(dir_in, files[i])
  file_out <- paste0(gsub(pattern="augch",
                          replacement="stoch",
                          x=gsub(pattern="trees",
                                 replacement="",
                                 x=file_in)),
                     "log")
  simmap_to_stochmap(file_in, file_out, tree)
}


# simmap object to rb format
# simulation is in state 1234, not 0123
dir_in = "data/2_simulation/missing_state/"
files <- list.files(dir_in, recursive = TRUE)[grepl(".Rda", list.files(dir_in, recursive = TRUE))]
for (i in 1:length(files)){
  file_in <- ""
  file_out <- paste0(dir_in, gsub(pattern="history.Rda", replacement="true_history.log", x=files[i]))
  load(paste0(dir_in, files[i]))
  simmaps <- list()
  simmaps[[1]] <- history2
  simmap_to_stochmap(file_in, file_out, tree, simmaps=simmaps)
}

