library(slouch)
library(tidyverse)
data("neocortex")
neocortex <- neocortex %>% select(species, diet)

# relative freq of discrete states in empirical data
length(which(neocortex$diet == "MF"))
length(which(neocortex$diet == "Br"))
length(which(neocortex$diet == "Gr"))

convert_state <- function(x){
  if (x <= 22){
    y = "2"
  } else if (x >= 23 & x <= 35){
    y = "0"
  } else {
    y = "1"
  }
  return(y)
}

for (i in 1:50){
  for (j in 1:43){
    random <- rdunif(43,1,43)
    states <- sapply(random, convert_state)
    neocortex[i+2] <- sapply(random, convert_state)
  }
}

for (i in 1:43){
  all_states <- paste(neocortex[1,3:52], collapse = '')
  neocortex$sim[i] <- all_states
}

  
lines = readLines("data/1_validation/artiodactyla/artiodactyla_Discrete.nex")
lines[3] <- "  DIMENSIONS NTAX=43 NCHAR=1;"
# make all the job scripts
bar = txtProgressBar(style=3, width=40)
for(i in 1:50) {
  these_lines = lines
  
  this_sim = neocortex[,i+2]
  this_species = neocortex[,1]
  
  for (j in 1:43){
    these_lines[5+j] <- paste0("    ", this_species[j], "        ", this_sim[j])
  }
  
  this_file = paste0('data/1_validation/artiodactyla/random_diets_sim_', i, '_Discrete.nex')
  
  cat(these_lines, file = this_file, sep="\n")
  
  setTxtProgressBar(bar, i / 50)
  
}


# causing segmentation error when clamping
#bar = txtProgressBar(style=3, width=40)
#
#lines[3] <- "  DIMENSIONS NTAX=43 NCHAR=50;"
#for (i in 1:43){
#  this_species = neocortex$species[i]
#  these_states = neocortex$sim[i]
#  
#  lines[5+i] <- paste0("    ", this_species, "        ", these_states)
#}
#
#this_file = paste0('data/1_validation/artiodactyla/sim_random_diets_Discrete.nex')
#
#cat(lines, file = this_file, sep="\n")

lines = readLines("jobs/artiodactyla_random_diets_sim_1.Rev")
# make all the job scripts
for(i in 1:50) {
  these_lines = lines
  these_lines[26] <- paste0('disc <- readDiscreteCharacterData("data/1_validation/artiodactyla/random_diets_sim_', i, '_Discrete.nex")')
  these_lines[123] <- paste0('monitors.append( mnModel(filename=\"output/1_validation/artiodactyla_random_diets/', i, '.log\", printgen=10) )')
  
  this_file = paste0('jobs/artiodactyla_random_diets_sim_', i, '.Rev')
  cat(these_lines, file = this_file, sep="\n")
}


