library(stringr)

numChangePosterior <- function(simmap_path, trace_path){
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


simmap_path <- "output/3_empirical/ase_carnivory/simmap.log"
trace_path <- "output/3_empirical/ase_carnivory/trace.log"

numChangePosterior(simmap_path, trace_path)

