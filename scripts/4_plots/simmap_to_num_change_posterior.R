library(stringr)

numChangePosterior <- function(simmap_path=NULL, trace_path=NULL){
  df <- read_tsv(simmap_path) %>%
    select(!c(head(colnames(.), 1), tail(colnames(.), 1)))
  
  message("Calculating total number of changes")
  num_change <- c()
  for (i in 1:nrow(df)){
    num_change[i] <- sum(sapply(df[i,], function(item){str_count(item, ":")}))
  }
  
  trace <- read_tsv(trace_path)
  trace$total_num_change <- num_change
  write_tsv(trace, trace_path)
  
  return(num_change)
}

#dir_in <- "output/3_empirical/ase_r500_3StateOrdered/"
dir_in <- "output/3_empirical/ase_r500_4State/"

simmap_paths <- paste0(dir_in, list.files(dir_in, pattern = "simmap"))
trace_paths <- paste0(dir_in, list.files(dir_in, pattern = "trace"))

for (i in 1:length(simmap_paths)){
  numChangePosterior(simmap_paths[i], trace_paths[i])
}

