library(ape)
library(phytools)
library(RevGadgets)
library(tidyverse)

readCharacterData = function(path, exclude) {
  
  # read the file
  lines = readLines(path)
  
  # find lines corresponding to data
  start_line = grep("MATRIX", toupper(lines), fixed=TRUE) + 1
  end_line   = grep(";", lines, fixed=TRUE)
  end_line   = end_line[which(end_line > start_line)[1]] - 1
  data_lines = lines[start_line:end_line]  
  
  # strip out white space and tabs
  data_lines = gsub("\t"," ", data_lines)
  data_lines = strsplit(data_lines, " ")
  data_lines = lapply(data_lines, function(x) x[x != ""])
  
  # get the taxa
  taxa = sapply(data_lines, function(x) x[1])
  data = do.call(rbind, lapply(data_lines, function(x) as.numeric(x[-1])))
  rownames(data) = taxa

  # excude the specified traits
  data = data[,!1:ncol(data) %in% exclude]
  
  return(data)  
  
}

writeCharacterData <- function(data, file, type) {
  
  lines <- "#NEXUS"
  lines <- c(lines, c(""))
  lines <- c(lines, "BEGIN DATA;")
  lines <- c(lines, paste0("\tDIMENSIONS NTAX=",nrow(data)," NCHAR=",ncol(data),";"))
  if( toupper(type) == "STANDARD") {
    lines <- c(lines, paste0("\tFORMAT DATATYPE=STANDARD SYMBOLS=\"",paste0(sort(unique(as.vector(data))),collapse=""),"\" MISSING=? GAP=-;"))
  } else if ( toupper(type) == "CONTINUOUS" ) {
    lines <- c(lines, paste0("\tFORMAT DATATYPE=CONTINUOUS MISSING=? GAP=-;"))
  } else {
    stop("Invalid data type")
  }
  lines <- c(lines,"MATRIX")
  
  taxa <- rownames(data)
  for(t in taxa) {
    these_traits <- as.numeric(data[t,])
    lines <- c(lines, paste(t, paste0(these_traits, collapse="\t"), sep="\t"))
  }
  
  lines <- c(lines, ";")
  lines <- c(lines, "END;")
  
  # cat(lines, sep="\n")
  
  cat(lines, file=file, sep="\n")
  
}

drawHalflife <- function(state_dependent, num_state, root_age){
  if(state_dependent == T){
    halflife <- runif(n=num_state, 0.1*root_age, root_age)
  } else {
    halflife <- rep(runif(n=1, 0.1*root_age, root_age), num_state)
  }
  names(halflife) = c(1:(num_state))
  return(halflife)
}

drawStv <- function(state_dependent, num_state, emp){
  if(isTRUE(state_dependent)){
    stv <- runif(n=num_state, 0.5*emp, 2*emp)
  } else {
    stv <- rep(runif(n=1, 0.5*emp, 2*emp), num_state)
  }
  names(stv) = c(1:(num_state))
  return(stv)
}


drawTheta <- function(linked, state_dependent, num_state){
  if (isTRUE(linked)){
    # can it be extended to more than 2 states?
    theta <- c(runif(1, -10, 10))
    if(num_state == 2){
      x <- rnorm(1, 0, 4)
      #while (abs(x+theta[1]) > 8){
      #  x <- rnorm(1, 0, 4)
      #}
      theta[2] <- theta[1] + x
    } else if (num_state == 3){
      x <- rnorm(1, 0, 4)
      #while (abs(x+theta[1]) > 8 | abs(x-theta[1]) > 8){
      #  x <- rnorm(1, 0, 4)
      #}
      theta[2] <- theta[1] + x
      theta[3] <- theta[1] - x
    } else {
      print(paste("Linked theta only supports 2 and 3-state currently."))
    }
  } else {
    if(isTRUE(state_dependent)){
      theta <- c(runif(n=num_state, -10, 10))
    } else {
      theta <- rep(runif(n=1, -10, 10), num_state)
    }
  }
  names(theta) = c(1:(num_state))
  return(theta)
}

simulateContinuous = function(tree, alpha, sigma2, theta) {
  preorder <- rev(postorder(tree))
  edges <- tree$edge
  root_node <- length(tree$tip.label) + 1
  state = tree$node.states[root_node]
  mu_at_nodes <- rep(0, length(tree$node.states))
  mu_at_nodes[root_node] <- theta[[state]]
  
  for (edge_index in preorder){
    sub_edges <- tree$maps[[edge_index]]
    parent_node <- edges[edge_index, 1]
    y <- mu_at_nodes[parent_node]
    for (j in 1:length(sub_edges)) {
      #alpha <- drawAlpha(state_dependent = stateDependencies[1])
      #sigma2 <- drawAlpha(state_dependent = stateDependencies[2])
      #theta <- drawAlpha(state_dependent = stateDependencies[3])
      
      state <- names(sub_edges[j])
      mu <- y * exp(-alpha[[state]] * sub_edges[[j]]) + theta[[state]] * (1 - exp(-alpha[[state]] * sub_edges[[j]]))
      v <- sigma2[[state]] / (2 * alpha[[state]]) * (1 - exp(-2 * alpha[[state]] * sub_edges[[j]]))
      y <- rnorm(n=1, mu, sqrt(v))    
    }
    desc_node <- edges[edge_index, 2]
    mu_at_nodes[desc_node] <- y
  }
  
  cont_list <- list()
  for (i in 1:length(tree$tip.label)){
    tip <- tree$tip.label[i]
    cont_list[[tip]] <- mu_at_nodes[i]
  }
  
  return(cont_list)
}

log_to_simmap <- function(files){
  for (file in files){
    path = paste0(dir_in, file)
    log <- read_tsv(path)
    log <- as.data.frame(log$char_hist)
    write_tsv(log, file=path, col_names = FALSE)
  }
}

simmap_to_ancStates <- function(input_path, output_path, tree){
  index_to_rev <- RevGadgets::matchNodes(tree)
  
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
