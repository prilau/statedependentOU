library(RevGadgets)
library(tidyverse)
tree <- read.tree("data/3_empirical/mammal_perMY_r500.tre")

index_to_rev <- RevGadgets::matchNodes(tree)

dir_in = "output/3_empirical/aug_r500_3StateOrderedModel/"
files <- list.files(dir_in)[grepl(".trees", list.files(dir_in))]
for (file in files){
  path = paste0(dir_in, file)
  log <- read_tsv(path)
  log <- as.data.frame(log$char_hist)
  write_tsv(log, file=path, col_names = FALSE)
}


simmap_to_ancStates <- function(input_path, output_path){
  #input_path <- paste0(dir_in, list.files(dir_in))
  simmaps <- read.simmap(input_path, format="phylip")
  
  df_rev <- data.frame()
  
  for (row_num in 1:length(simmaps)){
    simmap <- simmaps[[row_num]]
    
    # Iteration column
    df_rev[row_num, 1] <- row_num-1
    
    for (i in 1:(length(simmap$maps))){
      ape_node <- which(index_to_rev[,2]==i)
      ape_edge <- which(simmap$edge[,2]==ape_node)
      map <- simmap$maps[[ape_edge]]
      df_rev[row_num, i+1] <- names(map[length(map)])
    }
    df_rev[row_num, length(simmap$maps)+2] <- names(map[1])
  }
  
  header <- paste0("end_", as.character(1:(length(simmap$maps)+1)))
  colnames(df_rev) <- c("Iteration", header)
  write_tsv(df_rev, output_path)
}

for (i in 1:length(files)){
  file_in <- paste0(dir_in, "augch_run_", i, ".trees")
  file_out <- paste0(dir_in, "states_run_", i, ".log")
  simmap_to_ancStates(file_in, file_out)
}


#Do this in RevBayes!
#tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")[1]
#files_in <- ["output/3_empirical/aug_highBranchMoves/states_run_1.log",
#             "output/3_empirical/aug_highBranchMoves/states_run_2.log",
#             "output/3_empirical/aug_highNodeMoves/states_run_1.log",
#             "output/3_empirical/aug_highNodeMoves/states_run_2.log"]
#files_out <- ["output/3_empirical/aug_highBranchMoves/anc_states_run_1.tre",
#              "output/3_empirical/aug_highBranchMoves/anc_states_run_2.tre",
#              "output/3_empirical/aug_highNodeMoves/anc_states_run_1.tre",
#              "output/3_empirical/aug_highNodeMoves/anc_states_run_2.tre"]
#i=1
#for (file in files_in){
#  anc_states = readAncestralStateTrace(file)
#  anc_tree = ancestralStateTree(
#  tree=tree,
#  ancestral_state_trace_vector=anc_states,
#  include_start_states=false,
#  file=files_out[i],
#  summary_statistic="MAP",
#  reconstruction="marginal",
#  burnin=0.1,
#  nStates=3,
#  site=1)
#  i+=1
# }

ase <- processAncStates("output/3_empirical/ase_r500_3StateOrdered/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivorous",
                                       "1"="Omnivorous",
                                       "2"="Carnivorous"))
aug <- processAncStates("output/3_empirical/aug_r500_3StateModel/anc_states_run_1.tre",
                        state_labels=c("0"="Herbivorous",
                                       "1"="Omnivorous",
                                       "2"="Carnivorous"))


ggplot() +
  geom_point(aes(x = ase@data$anc_state_1_pp, y = aug@data$anc_state_1_pp))

