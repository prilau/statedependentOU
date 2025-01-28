library(ape)
library(phytools)
library(RevGadgets)
library(tidyverse)
tree <- read.tree("data/3_empirical/mammal_perMY_r500.tre")

index_to_rev <- RevGadgets::matchNodes(tree)
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

dir_in = "output/3_empirical/sdOU_r500_3StateOrderedModel/"
files <- list.files(dir_in, recursive = TRUE)[grepl(".trees", list.files(dir_in, recursive = TRUE))]
for (file in files){
  path = paste0(dir_in, file)
  log <- read_tsv(path)
  log <- as.data.frame(log$char_hist)
  write_tsv(log, file=path, col_names = FALSE)
}



for (i in 1:length(files)){
  file_in <- paste0(dir_in, files[i])
  file_out <- paste0(gsub(pattern="trees", replacement="", x=file_in), "log")
  simmap_to_ancStates(file_in, file_out)
}


#Do this in RevBayes!
#tree <- readTrees("data/3_empirical/mammal_perMY_r500.tre")[1]
#files_in <-  ["output/2_simulation/missing_state/sim_1/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_1/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_1/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_2/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_2/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_2/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_3/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_3/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_3/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_4/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_4/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_4/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_5/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_5/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_5/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_6/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_6/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_6/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_7/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_7/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_7/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_8/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_8/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_8/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_9/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_9/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_9/augch_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_10/augch_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_10/augch_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_10/augch_nstate_4_run_1.log"]
#files_out <- ["output/2_simulation/missing_state/sim_1/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_1/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_1/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_2/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_2/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_2/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_3/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_3/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_3/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_4/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_4/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_4/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_5/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_5/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_5/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_6/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_6/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_6/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_7/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_7/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_7/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_8/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_8/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_8/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_9/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_9/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_9/anc_states_nstate_4_run_1.log",
#              "output/2_simulation/missing_state/sim_10/anc_states_nstate_2_run_1.log",
#              "output/2_simulation/missing_state/sim_10/anc_states_nstate_3_run_1.log",
#              "output/2_simulation/missing_state/sim_10/anc_states_nstate_4_run_1.log"]
#i=1
# # nStates always >=3 !!!!
#for (file in files_in){
#  anc_states = readAncestralStateTrace(file)
#  anc_tree = ancestralStateTree(
#  tree=tree,
#  ancestral_state_trace_vector=anc_states,
#  include_start_states=false,
#  file=files_out[i],
#  summary_statistic="MAP",
#  reconstruction="marginal",
#  burnin=0.0,
#  nStates=4,
#  site=1)
#  i+=1
# }

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

