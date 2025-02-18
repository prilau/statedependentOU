library(ape)
library(phytools)
library(RevGadgets)
library(tidyverse)
library(pracma)
library(slouch)
library(tibble) 

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

log_to_simmap <- function(dir_in, files){
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


simmap_to_stochmap <- function(input_path=NULL, output_path, tree, simmaps=NULL){
  index_to_rev <- RevGadgets::matchNodes(tree)
  if (is.null(simmaps)){
    simmaps <- read.simmap(input_path, format="phylip")
    }
  df_rev <- data.frame()
  
  for (row_num in 1:length(simmaps)){
    simmap <- simmaps[[row_num]]
    
    # Iteration column
    df_rev[row_num, 1] <- ( row_num-1 ) * 10
    
    for (i in 1:(length(simmap$maps))){
      ape_node <- which(index_to_rev[,2]==i)
      ape_edge <- which(simmap$edge[,2]==ape_node)
      map <- simmap$maps[[ape_edge]]
      # {2,10.5815:1,15.7765}
      x <- ""
      for (j in 1:length(map)){
        x <- paste0(names(map[j]), ",", unname(map[j]), ":", x)
      }
      x <- gsub('.{0,1}$', '', x)
      x <- paste0("{", x, "}")
      df_rev[row_num, i+1] <- x
    }
    df_rev[row_num, length(simmap$maps)+2] <- paste0("{", names(map[length(map)]), ",0}")
  }
  
  header <- as.character(1:(length(simmap$maps)+1))
  colnames(df_rev) <- c("Iteration", header)
  write_tsv(df_rev, output_path)
}

matchNodes = function(phy) {
  
  # get some useful info
  num_tips = length(phy$tip.label)
  num_nodes = phy$Nnode
  tip_indexes = 1:num_tips
  node_indexes = num_tips + num_nodes:1
  
  node_map = data.frame(R=1:(num_tips + num_nodes), Rev=NA, visits=0)
  current_node = phy$Nnode + 2
  k = 1
  t = 1
  
  while(TRUE) {
    
    if ( current_node <= num_tips ) {
      node_map$Rev[node_map$R == current_node] = t
      current_node = phy$edge[phy$edge[,2] == current_node,1]
      t = t + 1
    } else {
      
      if ( node_map$visits[node_map$R == current_node] == 0 ) {
        node_map$Rev[node_map$R == current_node] = node_indexes[k]
        k = k + 1
      }
      node_map$visits[node_map$R == current_node] = node_map$visits[node_map$R == current_node] + 1
      
      if ( node_map$visits[node_map$R == current_node] == 1 ) {
        # go right
        current_node = phy$edge[phy$edge[,1] == current_node,2][2]
      } else if ( node_map$visits[node_map$R == current_node] == 2 ) {
        # go left
        current_node = phy$edge[phy$edge[,1] == current_node,2][1]
      } else if ( node_map$visits[node_map$R == current_node] == 3 ) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node = phy$edge[phy$edge[,2] == current_node,1]
        }
      }
    }
  }
  return(node_map[,1:2])
}

processAncStates <-
  function(path,
           state_labels = NULL,
           labels_as_numbers = FALSE,
           missing_to_NA = TRUE) {
    # read in tree
    tree <- readTrees(path)
    t <- tree[[1]][[1]]
    
    # process column names
    include_start_states <- FALSE
    if ("anc_state_1" %in% names(t@data)) {
      # do nothing
    } else if ("start_state_1" %in% names(t@data) &&
               "end_state_1" %in% names(t@data)) {
      include_start_states <- TRUE
    } else {
      stop(
        "tree file does not contain expected state labels:
                [\'anc_state\'] or [\'start_state\' and \'end_state\']"
      )
    }
    
    # add state labels
    n_states <- length(state_labels)
    t <-
      .assign_state_labels(t,
                           state_labels,
                           include_start_states,
                           labels_as_numbers,
                           missing_to_NA,
                           n_states)
    
    # add range for pp factors
    t <- .set_pp_factor_range(t, include_start_states)
    
    # return processed TreeIO object
    return(t)
    
  }

plotAncStatesMAP <- function(t,
                             # option for plotting shoulder states
                             cladogenetic = FALSE,
                             
                             # label taxa at tips
                             tip_labels = TRUE,
                             tip_labels_size = 2,
                             tip_labels_offset = 1,
                             tip_labels_italics = FALSE,
                             tip_labels_formatted = FALSE,
                             tip_labels_remove_underscore = TRUE,
                             
                             # label states at tips
                             tip_labels_states = FALSE,
                             tip_labels_states_size = 2,
                             tip_labels_states_offset = 0.1,
                             
                             # text labels at nodes
                             node_labels_as = NULL,
                             node_labels_size = 2,
                             node_labels_offset = 0.1,
                             node_labels_centered = FALSE,
                             
                             # what to plot at nodes
                             node_size_as = "state_posterior",
                             node_color_as = "state",
                             node_shape_as = NULL,
                             
                             # aesthetics for plotting at nodes
                             node_shape = 19,
                             node_color = "default",
                             node_size = c(2, 6),
                             
                             # aesthetics for tip states (inherents additional
                             # aesthetics from nodes)
                             tip_states = TRUE,
                             tip_states_size = node_size,
                             tip_states_shape = node_shape,
                             
                             state_transparency = 0.75,
                             tree_layout = "rectangular",
                             
                             timeline = FALSE,
                             geo = timeline,
                             geo_units = list("epochs", "periods"),
                             time_bars = timeline,
                             
                             tree_color = "black",
                             tree_linewidth = 1,
                             ...) {
  ##### parameter compatability checks! #####
  if (!methods::is(t, "treedata"))
    stop("t should be a treedata objects")
  if (is.logical(cladogenetic) == FALSE)
    stop("cladogenetic should be TRUE or FALSE")
  if (is.logical(tip_labels) == FALSE)
    stop("tip_labels should be TRUE or FALSE")
  if (is.numeric(tip_labels_size) == FALSE)
    stop("tip_labels_size should be a number")
  if (is.numeric(tip_labels_offset) == FALSE)
    stop("tip_labels_offset should be a number")
  if (is.logical(tip_labels_italics) == FALSE)
    stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_formatted) == FALSE)
    stop("tip_labels_formatted should be TRUE or FALSE")
  if (tip_labels_italics == TRUE & tip_labels_formatted == TRUE) 
    stop("tip_labels_italics and tip_labels_formatted may not both be TRUE")
  if (is.logical(tip_labels_remove_underscore) == FALSE)
    stop("tip_labels_remove_underscore should be TRUE or FALSE")
  if (is.logical(tip_labels_states) == FALSE)
    stop("tip_labels_states should be TRUE or FALSE")
  if (is.numeric(tip_labels_states_size) == FALSE)
    stop("tip_labels_states_size should be a number")
  if (is.numeric(tip_labels_states_offset) == FALSE)
    stop("tip_labels_states_offsetshould be a number")
  if (is.null(node_labels_as) == FALSE) {
    node_labels_as <-
      match.arg(node_labels_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.numeric(node_labels_size) == FALSE)
    stop("node_labels_size should be a number")
  if (is.numeric(node_labels_offset) == FALSE)
    stop("node_labels_offset should be a number")
  if (is.logical(node_labels_centered) == FALSE)
    stop("node_labels_centered should be TRUE or FALSE")
  if (is.null(node_size_as) == FALSE) {
    node_size_as <-
      match.arg(node_size_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_color_as) == FALSE) {
    node_color_as <-
      match.arg(node_color_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as != "state")
      stop("node_shape_as should be NULL or 'state'")
  }
  if (is.numeric(node_shape) == FALSE)
    stop("node_shape should be a number indicating symbol type")
  if (is.character(node_color) == FALSE)
    stop ("node_color should be 'default' or valid color(s)")
  if (node_color[1] != "default" &
      any(.isColor(node_color) == FALSE))
    stop("node_color should be valid color(s)")
  if (any(is.numeric(node_size) == FALSE))
    stop("node_size should be a single number or a vector of two numbers")
  if (length(node_size) > 2)
    stop("node_size should be a single number or a vector of two numbers")
  if (is.logical(tip_states) == FALSE)
    stop("tip_states should be TRUE or FALSE")
  if (is.numeric(tip_states_size) == FALSE)
    stop("tip_states_size should be a number")
  if (is.numeric(tip_states_shape) == FALSE)
    stop("tip_states_shape should be a number indicating symbol type")
  if (is.numeric(state_transparency) == FALSE)
    stop("state_transparency should be a number between 0 - 1")
  if (state_transparency > 1 |
      state_transparency < 0)
    stop("state_transparency should be a number between 0 - 1")
  if (cladogenetic == FALSE) {
    tree_layout <-
      match.arg(
        tree_layout,
        choices = c(
          'rectangular',
          'slanted',
          'ellipse',
          'roundrect',
          'fan',
          'circular',
          'inward_circular',
          'radial',
          'equal_angle',
          'daylight',
          'ape'
        )
      )
  } else if (cladogenetic == TRUE) {
    tree_layout <-
      match.arg(tree_layout, choices = c('rectangular', 'circular'))
  }
  if (is.logical(timeline) == FALSE)
    stop("timeline should be TRUE or FALSE")
  if (tree_layout != "rectangular") {
    if (timeline == TRUE) {
      stop("timeline is only compatible with
                                 tree_layout = 'rectangular'")
    }
  }
  if (is.list(geo_units)) {
    if (length(geo_units) != 2)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[1]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[2]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  } else {
    if (geo_units %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras') == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  }
  ##### calculate helper variables #####
  tree <- attributes(t)$phylo
  
  ##### create basic tree plot #####
  p <- ggtree::ggtree(t, layout = tree_layout, color = tree_color, linewidth = tree_linewidth, ...)
  
  # get dimensions
  n_node <- ape::Nnode(tree, internal.only = FALSE)
  tree_height <- max(phytools::nodeHeights(t@phylo))
  ntips <- sum(p$data$isTip)
  
  ##### process column names #####
  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE &
             "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE &
             "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }
  
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      if (!is.factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))) {
        p$data$node_color_as <-
          factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))
      } else {
        p$data$node_color_as <-
          dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))
      }
      # double check levels are alphabetical (if not numeric)
      if (suppressWarnings(any(is.na(as.integer(levels(p$data$node_color_as)))))) {
        levels(p$data$node_color_as) <-
          sort(levels(p$data$node_color_as))
      }
      
    }
    if (node_color_as == "node_posterior") {
      p$data$node_color_as <- as.numeric(p$data$posterior)
    }
    if (node_color_as == "state_posterior") {
      p$data$node_color_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }
  
  if (is.null(node_size_as) == FALSE) {
    if (node_size_as == "state") {
      size_tmp <- dplyr::pull(p$data, paste0(state_pos_str_base[1], "1"))
      if (is.factor(size_tmp)) {
        p$data$node_size_as <- as.integer(levels(size_tmp))[size_tmp]
      } else {
        p$data$node_size_as <- as.integer(size_tmp)
      }
    }
    if (node_size_as == "node_posterior") {
      p$data$node_size_as <- as.numeric(p$data$posterior)
    }
    if (node_size_as == "state_posterior") {
      p$data$node_size_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }
  
  if (is.null(node_shape_as) == FALSE) {
    if (node_shape_as == "state") {
      
      p$data$node_shape_as <-
        factor(dplyr::pull(p$data, paste0(state_pos_str_base[1], "1")))
    }
    if (node_shape_as == "node_posterior") {
      p$data$node_shape_as <- as.numeric(p$data$posterior)
    }
    if (node_shape_as == "state_posterior") {
      p$data$node_shape_as <-
        as.numeric(dplyr::pull(p$data,
                               paste0(state_pos_str_base[1], "1", "_pp")))
    }
  }
  
  if (cladogenetic == TRUE) {
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state") {
        if (!is.factor(dplyr::pull(p$data,
                                   paste0(state_pos_str_base[2], "1")))) {
          p$data$clado_node_color_as <-
            factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))
        } else {
          p$data$clado_node_color_as <-
            dplyr::pull(p$data, paste0(state_pos_str_base[2], "1"))
        }
        
        if (suppressWarnings(any(is.na(as.integer(levels(p$data$clado_node_color_as)))))) {
          levels(p$data$clado_node_color_as) <-
            sort(levels(p$data$clado_node_color_as))
        }
      }
      if (node_color_as == "node_posterior") {
        p$data$clado_node_color_as <- 1
      }
      if (node_color_as == "state_posterior") {
        p$data$clado_node_color_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
    
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state") {
        clado_size_tmp <- dplyr::pull(p$data, paste0(state_pos_str_base[2], "1"))
        if (is.factor(clado_size_tmp)) {
          p$data$clado_node_size_as <- as.integer(levels(clado_size_tmp))[clado_size_tmp]
        } else {
          p$data$clado_node_size_as <- as.integer(clado_size_tmp)
        }
      }
      if (node_size_as == "node_posterior") {
        p$data$clado_node_size_as <- 1
      }
      if (node_size_as == "state_posterior") {
        p$data$clado_node_size_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
    
    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state") {
        p$data$clado_node_shape_as <-
          factor(dplyr::pull(p$data, paste0(state_pos_str_base[2], "1")))
      }
      if (node_shape_as == "node_posterior") {
        p$data$clado_node_shape_as <- 1
      }
      if (node_shape_as == "state_posterior") {
        p$data$clado_node_shape_as <-
          as.numeric(dplyr::pull(p$data,
                                 paste0(state_pos_str_base[2], "1", "_pp")))
      }
    }
  }
  
  # gather list of all character states from data
  if (cladogenetic == TRUE) {
    all_states <- unique(c(p$data$start_state_1, p$data$end_state_1))
  } else {
    all_states <-
      na.omit(unique(factor(dplyr::pull(
        p$data, paste0(state_pos_str_base[1], "1")
      ))))
  }
  all_states <- sort(all_states)
  
  ##### color processing and checks #####
  # check if number of states exceeds default color palette options
  if (!is.null(node_color_as) && node_color_as == "states") {
    if (node_color[1] == "default") {
      nstates <- length(all_states)
      if (nstates <= 12) {
        node_color <- colFun(nstates)
      } else {
        node_color <- grDevices::colorRampPalette(colFun(12))(nstates)
      }
    }
    
    # check if number of states not equal to provided colors
    if (node_color[1] != "default" &
        length(node_color) < length(all_states)) {
      stop(
        paste0(
          "You provided fewer colors in node_color than states
          in your dataset. There are ",
          length(all_states),
          " states and you provide ",
          length(node_color),
          " colors."
        )
      )
    }
    if (node_color[1] != "default" &
        length(node_color) > length(all_states)) {
      stop(
        paste0(
          "You provided more colors in node_color than states
          in your dataset. There are ",
          length(all_states),
          " states and you provide ",
          length(node_color),
          " colors."
        )
      )
    }
  }
  
  # set default colors
  if (any(node_color == "default")) {
    if (is.null(node_color_as) == TRUE) {
      colors <- colFun(1)
    } else if (node_color_as == "state") {
      nstates <- length(all_states)
      colors <- colFun(nstates)
      # name colors if unnamed
      names(colors) <- sort(all_states)
    } else if (node_color_as == "node_posterior" |
               node_color_as == "state_posterior") {
      colors <- colFun(2)
    }
  } else {
    colors <- node_color
  }
  
  
  ##### adjust aesthetics lengths if needed #####
  # shape
  if (is.null(node_shape_as) == TRUE) {
    if (length(node_shape) > 1) {
      node_shape <- node_shape[1]
    }
  }
  # color
  if (is.null(node_color_as) == TRUE) {
    if (length(colors) > 1) {
      colors <- colors[1]
    }
  }
  # size
  if (is.null(node_size_as) == TRUE) {
    if (length(node_size) > 1) {
      node_size <- node_size[1]
    }
  }
  ##### reformat labels if necessary #####
  if (tip_labels_remove_underscore) {
    p$data$label <- gsub("_", " ", p$data$label)
  }
  
  ##### get hjust values #####
  if (node_labels_centered) {
    hjust_node <- 0.5
    hjust_shoulder <- 0.5
  }
  if (!node_labels_centered) {
    hjust_node <- 0
    hjust_shoulder <- 1
  }
  ##### calculate cladogenetic plotting data #####
  if (cladogenetic == TRUE) {
    x <- ggtree::fortify(tree)$x
    y <- ggtree::fortify(tree)$y
    x_anc <- numeric(n_node)
    node_index <- numeric(n_node)
    for (i in 1:n_node) {
      if (.getParent(tree, i) != 0) {
        # if not the root, get the x coordinate for the parent node
        x_anc[i] <- x[.getParent(tree, i)]
        node_index[i] <- i
      }
    }
    shoulder_data <-
      data.frame(node = node_index,
                 x_anc = x_anc,
                 y = y)
    if (timeline == TRUE) {
      shoulder_data$x_anc <- shoulder_data$x_anc - tree_height
    }
    `%<+%` <- ggtree::`%<+%`
    p <- p %<+% shoulder_data
  }
  
  ##### start plotting #####
  
  # add timeline
  if (timeline == TRUE) {
    max_age <- tree_height
    
    if (max_age > 100) {
      interval <- 50
    } else {
      interval <- 10
    }
    dx <- max_age %% interval
    # set coordinates
    ### Fix the xlim and ylims - if no error bars, should be a function of
    ### max age and n nodes, respectively.
    ### If error bars, -x lim should be as old as the max of the error bar
    tick_height <- ntips / 100
    if (geo == TRUE) {
      #determine whether to include quaternary
      if (tree_height > 50) {
        skipit <- c("Quaternary", "Holocene", "Late Pleistocene")
      } else {
        skipit <- c("Holocene", "Late Pleistocene")
      }
      # add deep timescale
      if (length(geo_units) == 1) {
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height * 1.1, tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          height = grid::unit(4, "line"),
          skip = skipit,
          abbrv = FALSE,
          rot = 90,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      } else if (length(geo_units) == 2) {
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height * 1.05, tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          skip = skipit,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
        
      }
    }
    #add axis title
    p <- p + ggplot2::scale_x_continuous(name = "Age (Ma)",
                                         limits = c(-tree_height, tree_height /
                                                      2))
    p <- ggtree::revts(p)
    # add ma ticks and labels
    xline <- pretty(c(0, max_age))[pretty(c(0, max_age)) < max_age]
    df <-
      data.frame(
        x = -xline,
        y = rep(-tick_height * 5, length(xline)),
        vx = -xline,
        vy = rep(-tick_height * 5 + tick_height, length(xline))
      )
    
    p <-
      p + ggplot2::geom_segment(ggplot2::aes(
        x = 0,
        y = -tick_height * 5,
        xend = -max_age,
        yend = -tick_height * 5
      )) +
      ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      )) +
      ggplot2::annotate(
        "text",
        x = -rev(xline),
        y = -tick_height * 5 + tick_height * 2,
        label = rev(xline),
        size = tip_labels_size
      )
    
    # add vertical gray bars
    if (time_bars) {
      if (geo) {
        if ("epochs" %in% geo_units) {
          x_pos <- -rev(c(0, deeptime::get_scale_data("epochs")$max_age))
        } else {
          x_pos <-  -rev(c(0, deeptime::get_scale_data("periods")$max_age))
        }
      } else if (!geo) {
        x_pos <- -rev(xline)
      }
      for (k in 2:(length(x_pos))) {
        box_col <- "gray92"
        if (k %% 2 == 1)
          box_col <- "white"
        box <-
          ggplot2::geom_rect(
            xmin = x_pos[k - 1],
            xmax = x_pos[k],
            ymin = -tick_height * 5,
            ymax = ntips,
            fill = box_col
          )
        p <- gginnards::append_layers(p, box, position = "bottom")
      }
    }
    if (tip_labels) {
      # recenter legend
      tot <- max_age + tree_height / 2
      p <-
        p + ggplot2::theme(axis.title.x =
                             ggplot2::element_text(hjust =
                                                     max_age /  (2 * tot)))
    }
  }
  
  # add tip labels
  if (tip_labels == TRUE) {
    if (tip_labels_italics == TRUE) {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = paste0('italic(`', label, '`)')),
          size = tip_labels_size,
          offset = tip_labels_offset,
          parse = TRUE
        )
    }
    else if (tip_labels_formatted == TRUE ) {
      p <- p + ggtree::geom_tiplab(
        ggplot2::aes(label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        parse = TRUE
      )
    } else {
      p <-
        p + ggtree::geom_tiplab(size = tip_labels_size,
                                offset = tip_labels_offset)
    }
  }
  
  # add the tip states
  if (tip_states == TRUE) {
    # unless node size should vary by state, don't allow tip sizes to vary
    if (is.null(node_size_as) == TRUE || node_size_as != "state") {
      tip_states_size <- tip_states_size[1]
    }
    
    # vary tip symbols by color only
    # when shape is null and size is not state
    if (is.null(node_color_as) == FALSE) {
      if (node_color_as == "state" &
          is.null(node_shape_as) == TRUE &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state"))  {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(colour = node_color_as),
          size = tip_states_size,
          alpha = state_transparency,
          shape = tip_states_shape
        )
      }
    }
    
    # vary tip symbols by shape only
    # when shape is state, color is not state, size is not state
    if (is.null(node_shape_as) == FALSE) {
      if (node_shape_as == "state" &
          (is.null(node_color_as) == TRUE ||
           node_color_as != "state") &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state")) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(shape = node_shape_as),
          size = tip_states_size,
          alpha = state_transparency,
          color = colors
        )
      }
    }
    
    # vary tip symbol by shape and color
    # when shape is state, color is state, and size is anything but state
    if (is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      if (node_color_as == "state" &
          node_shape_as == "state" &
          (is.null(node_size_as) == TRUE ||
           node_size_as != "state")) {
        p <-  p + ggtree::geom_tippoint(
          ggplot2::aes(shape = node_shape_as,
                       color = node_color_as),
          size = tip_states_size,
          alpha = state_transparency
        )
      }
    }
    
    # vary tip symbol by size only
    # when size is state, color is not state, and shape is null
    if (is.null(node_size_as) == FALSE) {
      if (node_size_as == "state" &
          (is.null(node_color_as) == TRUE ||
           node_color_as != "state") &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(size = node_size_as),
          shape = tip_states_shape,
          alpha = state_transparency,
          color = "grey"
        )
      }
    }
    
    # vary tip symbol by size and color
    # when size is state, color is state or PP, and shape is null
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE) {
      if (node_size_as == "state" &
          is.null(node_shape_as) == TRUE) {
        p <- p + ggtree::geom_tippoint(
          ggplot2::aes(size = node_size_as,
                       color = node_color_as),
          shape = tip_states_shape,
          alpha = state_transparency
        )
      }
    }
  }
  
  # plot symbols at nodes and shoulders
  blank_nodes <-
    is.null(node_color_as) == TRUE &
    is.null(node_size_as) == TRUE & is.null(node_shape_as) == TRUE
  if (blank_nodes == FALSE) {
    # plot if color, size, and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <- p + ggtree::geom_nodepoint(
        ggplot2::aes(
          colour = node_color_as,
          size = node_size_as,
          shape = node_shape_as
        ),
        alpha = state_transparency
      )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    
    # plot if color and size vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as, size = node_size_as),
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if color and shape vary
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as, shape = node_shape_as),
          size = node_size,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              colour = clado_node_color_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if size and shape vary
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(shape = node_shape_as, size = node_size_as),
          color = colors,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              size = clado_node_size_as,
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if just color varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == FALSE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(colour = node_color_as),
          size = node_size,
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            size = node_size,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if just size varies
    if (is.null(node_size_as) == FALSE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == TRUE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(size = node_size_as),
          color = colors,
          shape = node_shape,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              size = clado_node_size_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              color = clado_node_color_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            shape = node_shape,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
    
    #plot if just shape varies
    if (is.null(node_size_as) == TRUE &
        is.null(node_color_as) == TRUE &
        is.null(node_shape_as) == FALSE) {
      p <-
        p + ggtree::geom_nodepoint(
          ggplot2::aes(shape = node_shape_as),
          color = colors,
          size = node_size,
          alpha = state_transparency
        )
      #add start (shoulder) states
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_nodepoint(
            ggplot2::aes(
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          ) +
          ggtree::geom_tippoint(
            ggplot2::aes(
              shape = clado_node_shape_as,
              x = x_anc,
              y = y
            ),
            color = colors,
            size = node_size,
            na.rm = TRUE,
            alpha = state_transparency
          )
      }
    }
  }
  
  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    if (node_labels_as == "state") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = start_state_1,
              x = x_anc,
              y = y
            ),
            hjust = hjust_shoulder,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = anc_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "state_posterior") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = .convertAndRound(start_state_1_pp),
              x = x_anc,
              y = y
            ),
            hjust = hjust_shoulder,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(anc_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = hjust_node,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "node_posterior") {
      p <-
        p + ggtree::geom_nodelab(
          ggplot2::aes(label = .convertAndRound(posterior)),
          hjust = hjust_node,
          nudge_x = node_labels_offset,
          size = node_labels_size
        )
    }
  }
  
  # add tip states labels (text)
  if (tip_labels_states == TRUE) {
    if (state_pos_str_base[1] == "anc_state_") {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = anc_state_1),
          hjust = hjust_node,
          offset = tip_labels_states_offset,
          size = tip_labels_states_size
        )
    } else {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = end_state_1),
          hjust = hjust_node,
          offset = tip_labels_states_offset,
          size = tip_labels_states_size
        )
    }
  }
  
  # add custom colors, shapes, and sizes
  if (is.null(node_size_as) == FALSE) {
    p <-
      p + ggplot2::scale_size(range = node_size,
                              name = .titleFormat(node_size_as))
  }
  if (is.null(node_color_as) == FALSE) {
    if (node_color_as == "state") {
      if (is.null(names(colors))) {
        p <- p + ggplot2::scale_color_manual(
          values = colors,
          na.translate = FALSE,
          name = .titleFormat(node_color_as)
        )
      } else {
        p <- p + ggplot2::scale_color_manual(
          values = colors,
          na.translate = FALSE,
          name = .titleFormat(node_color_as),
          breaks = names(colors)
        )
      }
      
    } else if (node_color_as == "state_posterior" |
               node_color_as == "node_posterior") {
      if (cladogenetic) {
        prettify <- c(p$data$node_color_as, p$data$clado_node_color_as)
      } else {
        prettify <- p$data$node_color_as
      }
      p <- p + ggplot2::scale_color_gradient(
        low = colors[1],
        high = colors[2],
        breaks = pretty(prettify),
        name = .titleFormat(node_color_as)
      )
    }
  }
  if (is.null(node_shape_as) == FALSE) {
    p <-
      p + ggplot2::scale_shape_manual(values = node_shape,
                                      name = .titleFormat(node_shape_as))
  }
  
  # add space on x axis for tip labels
  if (tip_labels == TRUE) {
    if (timeline == FALSE) {
      p <- p + ggtree::xlim(0, tree_height + tree_height / 2)
    }
  }
  
  return(p)
}

plotAncStatesPie <- function(t,
                             # option for plotting shoulder states
                             cladogenetic = FALSE,
                             
                             # label taxa at tips
                             tip_labels = TRUE,
                             tip_labels_size = 2,
                             tip_labels_offset = 1,
                             tip_labels_italics = FALSE,
                             tip_labels_formatted = FALSE,
                             tip_labels_remove_underscore = TRUE,
                             
                             # label states at tips
                             tip_labels_states = FALSE,
                             tip_labels_states_size = 2,
                             tip_labels_states_offset = 0.1,
                             
                             # text labels at nodes
                             node_labels_as = NULL,
                             node_labels_size = 2,
                             node_labels_offset = 0.1,
                             
                             # pies aesthetics
                             pie_colors = "default",
                             node_pie_size = 1,
                             shoulder_pie_size = node_pie_size,
                             tip_pies = TRUE,
                             tip_pie_size = 0.5,
                             
                             # nudges to center pies
                             node_pie_nudge_x = 0,
                             node_pie_nudge_y = 0,
                             tip_pie_nudge_x = node_pie_nudge_x,
                             tip_pie_nudge_y = node_pie_nudge_y,
                             shoulder_pie_nudge_x = node_pie_nudge_x,
                             shoulder_pie_nudge_y = node_pie_nudge_y,
                             
                             state_transparency = 0.75,
                             tree_layout = "rectangular",
                             
                             timeline = FALSE,
                             geo = timeline,
                             geo_units = list("epochs", "periods"),
                             time_bars = timeline,
                             
                             tree_color = "black",
                             tree_linewidth = 1,
                             ...) {
  ##### parameter compatibility checks #####
  if (!methods::is(t, "treedata"))
    stop("t should be a treedata object")
  if (is.logical(cladogenetic) == FALSE)
    stop("cladogenetic should be TRUE or FALSE")
  if (is.logical(tip_labels) == FALSE)
    stop("tip_labels should be TRUE or FALSE")
  if (is.numeric(tip_labels_size) == FALSE)
    stop("tip_labels_size should be a number")
  if (is.numeric(tip_labels_offset) == FALSE)
    stop("tip_labels_offset should be a number")
  if (is.logical(tip_labels_italics) == FALSE)
    stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_formatted) == FALSE)
    stop("tip_labels_formatted should be TRUE or FALSE")
  if (tip_labels_italics == TRUE & tip_labels_formatted == TRUE) 
    stop("tip_labels_italics and tip_labels_formatted may not both be TRUE")
  if (is.logical(tip_labels_remove_underscore) == FALSE)
    stop("tip_labels_remove_underscore should be TRUE or FALSE")
  if (is.logical(tip_labels_states) == FALSE)
    stop("tip_labels_states should be TRUE or FALSE")
  if (is.numeric(tip_labels_states_size) == FALSE)
    stop("tip_labels_states_size should be a number")
  if (is.numeric(tip_labels_states_offset) == FALSE)
    stop("tip_labels_states_offsetshould be a number")
  if (is.null(node_labels_as) == FALSE) {
    node_labels_as <-
      match.arg(node_labels_as,
                choices = c("state", "state_posterior", "node_posterior"))
  }
  if (is.numeric(node_labels_size) == FALSE)
    stop("node_labels_size should be a number")
  if (is.numeric(node_labels_offset) == FALSE)
    stop("node_labels_offset should be a number")
  if (is.character(pie_colors) == FALSE)
    stop ("pie_colors should be 'default' or valid color(s)")
  if (pie_colors[1] != "default" &
      any(.isColor(pie_colors) == FALSE))
    stop("pie_colors should be valid color(s)")
  if (any(is.numeric(node_pie_size) == FALSE))
    stop("node_pie_size should be a single number")
  if (length(node_pie_size) > 1)
    stop("node_pie_size should be a single number")
  if (any(is.numeric(shoulder_pie_size) == FALSE))
    stop("shoulder_pie_size should be a single number")
  if (length(shoulder_pie_size) > 1)
    stop("shoulder_pie_size should be a single number")
  if (any(is.numeric(tip_pie_size) == FALSE))
    stop("tip_pie_size should be a single number")
  if (length(tip_pie_size) > 1)
    stop("tip_pie_size should be a single number")
  if (is.numeric(node_pie_nudge_x) == FALSE)
    stop("node_pie_nudge_x should be a single number")
  if (is.numeric(node_pie_nudge_y) == FALSE)
    stop("node_pie_nudge_y should be a single number")
  if (is.numeric(tip_pie_nudge_x) == FALSE)
    stop("tip_pie_nudge_x should be a single number")
  if (is.numeric(tip_pie_nudge_y) == FALSE)
    stop("tip_pie_nudge_y should be a single number")
  if (is.numeric(shoulder_pie_nudge_x) == FALSE)
    stop("shoulder_pie_nudge_x should be a single number")
  if (is.numeric(shoulder_pie_nudge_y) == FALSE)
    stop("shoulder_pie_nudge_y should be a single number")
  if (length(node_pie_nudge_x) != 1)
    stop("node_pie_nudge_x should be a single number")
  if (length(node_pie_nudge_y) != 1)
    stop("node_pie_nudge_y should be a single number")
  if (length(tip_pie_nudge_x) != 1)
    stop("tip_pie_nudge_x should be a single number")
  if (length(tip_pie_nudge_y) != 1)
    stop("tip_pie_nudge_y should be a single number")
  if (length(shoulder_pie_nudge_x) != 1)
    stop("shoulder_pie_nudge_x should be a single number")
  if (length(shoulder_pie_nudge_y) != 1)
    stop("shoulder_pie_nudge_y should be a single number")
  if (is.numeric(state_transparency) == FALSE)
    stop("state_transparency should be a number between 0 - 1")
  if (state_transparency < 0 ||
      state_transparency > 1)
    stop("state_transparency should be a number between 0 - 1")
  if (cladogenetic == FALSE) {
    tree_layout <-
      match.arg(
        tree_layout,
        choices = c(
          'rectangular',
          'slanted',
          'ellipse',
          'roundrect',
          'fan',
          'circular',
          'inward_circular',
          'radial',
          'equal_angle',
          'daylight',
          'ape'
        )
      )
  } else if (cladogenetic == TRUE) {
    tree_layout <-
      match.arg(tree_layout, choices = c('rectangular', 'circular'))
  }
  if (is.logical(timeline) == FALSE)
    stop("timeline should be TRUE or FALSE")
  if (is.list(geo_units)) {
    if (length(geo_units) != 2)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[1]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[2]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  } else {
    if (geo_units %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras') == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  }
  
  ##### create basic tree plot #####
  p <- ggtree::ggtree(t, layout = tree_layout,
                      color = tree_color, linewidth = tree_linewidth, ...)
  
  #p <- ggtree::ggtree(t)
  
  ##### specify temp directory for intermediary files #####
  tmp <- tempdir()
  
  ##### calculate helper variables #####
  tree <- attributes(t)$phylo
  tree_height <- max(phytools::nodeHeights(t@phylo))
  n_node <- ape::Nnode(tree, internal.only = FALSE)
  ntips <- length(tree$tip.label)
  node_idx <- (ntips + 1):n_node
  tip_idx <- 1:ntips
  all_idx <- 1:n_node
  
  ##### transform nudge parameter #####
  tip_pie_nudge_x <- -tip_pie_nudge_x
  node_pie_nudge_x <- -node_pie_nudge_x
  shoulder_pie_nudge_x <- -shoulder_pie_nudge_x
  
  ##### reorder labels #####
  state_labels <- as.factor(attributes(t)$state_labels)
  
  ##### calculate pie sizes #####
  node_pie_size <-  node_pie_size / 30
  shoulder_pie_size <- shoulder_pie_size / 30
  tip_pie_size <- tip_pie_size / 30
  
  if (cladogenetic == TRUE) {
    state_pos_str_base <- c("end_state_", "start_state_")
  } else if (cladogenetic == FALSE &
             "start_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "end_state_"
  } else if (cladogenetic == FALSE &
             "anc_state_1" %in% colnames(p$data)) {
    state_pos_str_base <- "anc_state_"
  }
  
  ##### color and label processing #####
  
  # check if number of states exceeds default color palette options
  if (pie_colors[1] == "default") {
    nstates <- length(state_labels)
    if (nstates <= 12) {
      pie_colors <- colFun(nstates)
    } else {
      pie_colors <- grDevices::colorRampPalette(colFun(12))(nstates)
    }
  }
  
  # check if number of states not equal to provided colors
  if (pie_colors[1] != "default" &
      length(pie_colors) < length(state_labels)) {
    stop(
      paste0(
        "You provided fewer colors in node_color
                than states in your dataset. There are ",
        length(state_labels),
        " states and you provide ",
        length(pie_colors),
        " colors."
      )
    )
  }
  
  # add names to colors if none present
  if (is.null(names(pie_colors))) {
    names(pie_colors) <- state_labels
  }
  
  # set colors, add "other" if necessary
  otherpp <- as.numeric(dplyr::pull(p$data,
                                    var = paste0(state_pos_str_base[1],
                                                 "other_pp")))
  if (sum(otherpp, na.rm = TRUE) == 0) {
    # set default colors
    if (any(pie_colors == "default")) {
      nstates <- length(state_labels)
      colors <- c(colFun(nstates))
      names(colors) <- state_labels
    } else {
      colors <- pie_colors
    }
    
  } else if (sum(otherpp, na.rm = TRUE) != 0) {
    
    state_labels <- as.factor(c(as.character(t@state_labels), "other"))
    
    if ("anc_state_" %in% state_pos_str_base) {
      p$data$anc_state_other <- "other"
    }
    if ("end_state_" %in% state_pos_str_base) {
      p$data$end_state_other <- "other"
    }
    
    # add other to user-set colors
    if (pie_colors[1] != "default") {
      pc_names <- names(pie_colors)
      pie_colors <- c(pie_colors, "grey50")
      names(pie_colors) <- c(pc_names, "other")
    }
    
    # set default colors
    if (any(pie_colors == "default")) {
      nstates <- length(state_labels) - 1
      colors <- c(colFun(nstates), "grey50")
      names(colors) <- state_labels
    } else {
      colors <- pie_colors
    }
    
  }
  
  ##### reformat labels if necessary #####
  if (tip_labels_remove_underscore) {
    p$data$label <- gsub("_", " ", p$data$label)
  }
  
  ##### start plotting #####
  
  # add timeline
  if (timeline == TRUE) {
    
    max_age <- tree_height
    
    if (max_age > 100) {
      interval <- 50
    } else {
      interval <- 10
    }
    dx <- max_age %% interval
    # set coordinates
    ### Fix the xlim and ylims - if no error bars, should be a function of
    ### max age and n nodes, respectively.
    ### If error bars, -x lim should be as old as the max of the error bar.
    tick_height <- ntips / 100
    if (geo == TRUE) {
      #determine whether to include quaternary
      if (tree_height > 50) {
        skipit <- c("Quaternary", "Holocene", "Late Pleistocene")
      } else {
        skipit <- c("Holocene", "Late Pleistocene")
      }
      # add deep timescale
      if (length(geo_units) == 1) {
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height, tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          height = grid::unit(4, "line"),
          skip = skipit,
          abbrv = FALSE,
          rot = 90,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      } else if (length(geo_units) == 2) {
        p <- p + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-tree_height, tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          skip = skipit,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
        
      }
    }
    #add axis title
    p <- p + ggplot2::scale_x_continuous(name = "Age (Ma)",
                                         limits = c(-tree_height, tree_height /
                                                      2))
    p <- ggtree::revts(p)
    # add ma ticks and labels
    xline <- pretty(c(0, max_age))[pretty(c(0, max_age)) < max_age]
    df <-
      data.frame(
        x = -xline,
        y = rep(-tick_height * 5, length(xline)),
        vx = -xline,
        vy = rep(-tick_height * 5 + tick_height, length(xline))
      )
    
    p <-
      p + ggplot2::geom_segment(ggplot2::aes(
        x = 0,
        y = -tick_height * 5,
        xend = -max_age,
        yend = -tick_height * 5
      )) +
      ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      )) +
      ggplot2::annotate(
        "text",
        x = -rev(xline),
        y = -tick_height * 5 + tick_height * 2,
        label = rev(xline),
        size = tip_labels_size
      )
    
    # add vertical gray bars
    if (time_bars) {
      if (geo) {
        if ("epochs" %in% geo_units) {
          x_pos <- -rev(c(0, deeptime::get_scale_data("epochs")$max_age))
        } else {
          x_pos <-  -rev(c(0, deeptime::get_scale_data("periods")$max_age))
        }
      } else if (!geo) {
        x_pos <- -rev(xline)
      }
      for (k in 2:(length(x_pos))) {
        box_col <- "gray92"
        if (k %% 2 == 1)
          box_col <- "white"
        box <-
          ggplot2::geom_rect(
            xmin = x_pos[k - 1],
            xmax = x_pos[k],
            ymin = -tick_height * 5,
            ymax = ntips,
            fill = box_col
          )
        p <- gginnards::append_layers(p, box, position = "bottom")
      }
    }
    if (tip_labels) {
      # recenter legend
      tot <- max_age + tree_height / 2
      p <-
        p + ggplot2::theme(axis.title.x =
                             ggplot2::element_text(hjust = max_age /
                                                     (2 * tot)))
    }
  }
  
  # add tip labels
  if (tip_labels == TRUE) {
    if (tip_labels_italics == TRUE) {
      p <-
        p + ggtree::geom_tiplab(
          ggplot2::aes(label = paste0('italic(`', label, '`)')),
          size = tip_labels_size,
          offset = tip_labels_offset,
          parse = TRUE
        )
    } 
    else if (tip_labels_formatted == TRUE ) {
      p <- p + ggtree::geom_tiplab(
        ggplot2::aes(label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        parse = TRUE
      )
    } else {
      p <-
        p + ggtree::geom_tiplab(size = tip_labels_size,
                                offset = tip_labels_offset)
    }
  }
  
  # set up guides
  if (cladogenetic == TRUE) {
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_1),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_2),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_3),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(end_state_1),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(end_state_2),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(end_state_3),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    
    if ("other" %in% state_labels) {
      p <-
        p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_other),
                                               size = 0),
                                   na.rm = TRUE,
                                   alpha = 0.0)
    }
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(start_state_1),
                                             size =  0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(start_state_2),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(start_state_3),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    
  } else if (cladogenetic == FALSE &
             "anc_state_1" %in% colnames(t@data)) {
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(anc_state_1),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(anc_state_2),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    if ("anc_state_3" %in% colnames(t@data)) {
      p <-
        p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(anc_state_3),
                                               size = 0),
                                   na.rm = TRUE,
                                   alpha = 0.0)
    }
    
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(anc_state_1),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(anc_state_2),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    if ("anc_state_3" %in% colnames(t@data)){
      p <-
        p + ggtree::geom_tippoint(ggtree::aes(colour = factor(anc_state_3),
                                              size = 0),
                                  na.rm = TRUE,
                                  alpha = 0.0)
    }
    
    
    if ("other" %in% state_labels) {
      p <-
        p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(anc_state_other),
                                               size = 0),
                                   na.rm = TRUE,
                                   alpha = 0.0)
    }
    
  } else {
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_1),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    p <-
      p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_2),
                                             size = 0),
                                 na.rm = TRUE,
                                 alpha = 0.0)
    
    if ("end_state_3" %in% colnames(t@data)){
      p <-
        p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_3),
                                               size = 0),
                                   na.rm = TRUE,
                                   alpha = 0.0)
    }
    
    
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(end_state_1),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    p <-
      p + ggtree::geom_tippoint(ggtree::aes(colour = factor(end_state_2),
                                            size = 0),
                                na.rm = TRUE,
                                alpha = 0.0)
    if ("end_state_3" %in% colnames(t@data)){
      p <-
        p + ggtree::geom_tippoint(ggtree::aes(colour = factor(end_state_3),
                                              size = 0),
                                  na.rm = TRUE,
                                  alpha = 0.0)
    }
    
    if ("other" %in% state_labels) {
      p <-
        p + ggtree::geom_nodepoint(ggtree::aes(colour = factor(end_state_other),
                                               size = 0),
                                   na.rm = TRUE,
                                   alpha = 0.0)
    }
  }
  
  if (is.null(names(colors))) {
    breaks <- levels(state_labels)
  } else {breaks <- names(colors)}
  
  p <-
    p + ggplot2::scale_color_manual(values = colors, breaks = breaks)
  p <-
    p + ggplot2::guides(colour =
                          ggplot2::guide_legend("State",
                                                override.aes =
                                                  list(size = 4, alpha = 1.0)),
                        order = 1)
  p <- p + ggplot2::guides(size = "none")
  
  # import theme
  theme_transparent <- ggimage::theme_transparent()
  
  # plot pies at nodes (and shoulders)
  if (cladogenetic == TRUE) {
    
    # create state matrices (matrix of nodes (rows) and all
    # possible states (columns), values are pp. )
    state_probs <-
      .build_state_probs(t, state_labels, include_start_states = TRUE)
    dat_state_end <- state_probs$end
    dat_state_start <- state_probs$start
    
    # make pie plots
    pies_start <-
      .nodepie(
        dat_state_start,
        cols = 1:(ncol(dat_state_start) - 1),
        color = colors,
        alpha = state_transparency
      ) 
    pies_end <-
      .nodepie(
        dat_state_end,
        cols = 1:(ncol(dat_state_end) - 1),
        color = colors,
        alpha = state_transparency
      )
    
    # change 0s to avoid dividing by zero when calculating coordinates
    zeros <- which(dplyr::pull(p$data, "x") == 0)
    p$data[zeros, "x"] <- 0.0001
    
    # convert pie plots to lists
    
    # NODE PIES
    # save pies as images and plot as raster grobs
    pies_end_to_plot <- pies_end[node_idx]
    results_end <- list()
    for (i in seq_len(length(pies_end_to_plot))) {
      ggplot2::ggsave(
        paste0(tmp,"/.temp.png"),
        plot = pies_end_to_plot[[i]],
        bg = "transparent",
        width = 3,
        height = 3,
        units = "cm",
        dpi = 200
      )
      pie <- png::readPNG(paste0(tmp,"/.temp.png"))
      results_end[[i]] <-
        ggplotify::as.ggplot(grid::rasterGrob(pie, interpolate = TRUE))
    }
    # plotting data
    df_pies_end <- p$data[p$data$isTip == FALSE, ]
    # adjust nudges
    df_pies_end$x <- df_pies_end$x - node_pie_nudge_x
    df_pies_end$y <- df_pies_end$y - node_pie_nudge_y
    
    # SHOULDER PIES
    # save pies as images and plot as raster grobs
    results_start <- list()
    for (i in seq_len(length(pies_start))) {
      ggplot2::ggsave(
        paste0(tmp,"/.temp.png"),
        plot = pies_start[[i]],
        bg = "transparent",
        width = 3,
        height = 3,
        units = "cm",
        dpi = 200
      )
      pie <- png::readPNG(paste0(tmp,"/.temp.png"))
      results_start[[i]] <-
        ggplotify::as.ggplot(grid::rasterGrob(pie, interpolate = TRUE))
    }
    # plotting data
    df_pies_start <- p$data
    df_pies_start$x <- df_pies_start$x[match(df_pies_start$parent,
                                             df_pies_start$node)]
    # adjust nudges
    df_pies_start$x <- df_pies_start$x - shoulder_pie_nudge_x
    df_pies_start$y <- df_pies_start$y - shoulder_pie_nudge_y
    
    # save pies as images and plot as raster grobs
    # TIP PIES
    if (tip_pies == TRUE) {
      pies_tip <- pies_end[tip_idx]
      results_tip <- list()
      for (i in seq_len(length(pies_tip))) {
        ggplot2::ggsave(
          paste0(tmp,"/.temp.png"),
          plot = pies_end[[i]],
          bg = "transparent",
          width = 3,
          height = 3,
          units = "cm",
          dpi = 200
        )
        pie <- png::readPNG(paste0(tmp,"/.temp.png"))
        results_tip[[i]] <-
          ggplotify::as.ggplot(grid::rasterGrob(pie, interpolate = TRUE))
      }
      # plotting data
      df_pies_tip <- p$data[p$data$isTip == TRUE, ]
      # adjust nudges
      df_pies_tip$x <- df_pies_tip$x - tip_pie_nudge_x
      df_pies_tip$y <- df_pies_tip$y - tip_pie_nudge_y
    }
    
    if (tip_pies == TRUE) {
      df_pies <- rbind(df_pies_end, df_pies_start, df_pies_tip)
      results <- c(results_end, results_start, results_tip)
      sizes <- c(rep(node_pie_size, nrow(df_pies_end)),
                 rep(shoulder_pie_size, nrow(df_pies_start)),
                 rep(tip_pie_size, nrow(df_pies_tip)))
    } else {
      df_pies <- rbind(df_pies_end, df_pies_start)
      results <- c(results_end, results_start)
      sizes <- c(rep(node_pie_size, nrow(df_pies_end)),
                 rep(shoulder_pie_size, nrow(df_pies_start)))
    }
    
    p <-
      p + ggpp::geom_plot(data = df_pies,
                          mapping = ggplot2::aes(
                            x = x,
                            y = y,
                            label = results
                          ),
                          vp.width = sizes,
                          vp.height = sizes,
                          hjust = 0.5,
                          vjust = 0.5
      )
    
  } else {
    # create state matrices (matrix of nodes (rows) and all
    # possible states (columns), values are pp. )
    dat_state_anc <-
      .build_state_probs(t, state_labels, include_start_states = FALSE)[[1]]
    if (sum(otherpp, na.rm = TRUE) == 0) {
      dat_state_anc$other <- NULL
    }
    pies_anc <-
      .nodepie(
        dat_state_anc,
        cols = 1:(ncol(dat_state_anc) - 1),
        color = colors,
        alpha = state_transparency
      )
    zeros <- which(dplyr::pull(p$data, "x") == 0)
    p$data[zeros, "x"] <- 0.0001
    
    # convert pie plots to lists
    # NODE PIES
    # save pies as images and plot as raster grobs
    pies_anc_to_plot <- pies_anc[node_idx]
    results_anc <- list()
    for (i in seq_len(length(pies_anc_to_plot))) {
      ggplot2::ggsave(
        paste0(tmp,"/.temp.png"),
        plot = pies_anc_to_plot[[i]],
        bg = "transparent",
        width = 3,
        height = 3,
        units = "cm",
        dpi = 200
      )
      pie <- png::readPNG(paste0(tmp,"/.temp.png"))
      results_anc[[i]] <-
        ggplotify::as.ggplot(grid::rasterGrob(pie, interpolate = TRUE))
    }
    # plotting data
    df_pies_anc <- p$data[p$data$isTip == FALSE, ]
    # adjust nudges
    df_pies_anc$x <- df_pies_anc$x - node_pie_nudge_x
    df_pies_anc$y <- df_pies_anc$y - node_pie_nudge_y
    
    
    # TIP PIES
    if (tip_pies == TRUE) {
      pies_tip <- pies_anc[tip_idx]
      results_tip <- list()
      for (i in seq_len(length(pies_tip))) {
        ggplot2::ggsave(
          paste0(tmp,"/.temp.png"),
          plot = pies_tip[[i]],
          bg = "transparent",
          width = 3,
          height = 3,
          units = "cm",
          dpi = 200
        )
        pie <- png::readPNG(paste0(tmp,"/.temp.png"))
        results_tip[[i]] <-
          ggplotify::as.ggplot(grid::rasterGrob(pie, interpolate = TRUE))
      }
      # plotting data
      df_pies_tip <- p$data[p$data$isTip == TRUE, ]
      # adjust nudges
      df_pies_tip$x <- df_pies_tip$x - tip_pie_nudge_x
      df_pies_tip$y <- df_pies_tip$y - tip_pie_nudge_y
    }
    
    if (tip_pies == TRUE) {
      df_pies <- rbind(df_pies_anc, df_pies_tip)
      results <- c(results_anc, results_tip)
      sizes <- c(rep(node_pie_size, nrow(df_pies_anc)),
                 rep(tip_pie_size, nrow(df_pies_tip)))
    } else {
      df_pies <- df_pies_anc
      results <- results_anc
      sizes <- rep(node_pie_size, nrow(df_pies_anc))
    }
    
    
    # save pies as images and plot as raster grobs
    p <-
      p + ggpp::geom_plot(data = df_pies,
                          mapping = ggplot2::aes(
                            x = x,
                            y = y,
                            label = results
                          ),
                          vp.width = sizes,
                          vp.height = sizes,
                          hjust = 0.5,
                          vjust = 0.5
      )
  }
  
  # add node labels (text)
  if (is.null(node_labels_as) == FALSE) {
    # add clado plotting data for node labels
    if (cladogenetic == TRUE) {
      x <- ggtree::fortify(tree)$x
      y <- ggtree::fortify(tree)$y
      x_anc <- numeric(n_node)
      node_index <- numeric(n_node)
      for (i in 1:n_node) {
        if (.getParent(tree, i) != 0) {
          # if not the root, get the x coordinate for the parent node
          x_anc[i] <- x[.getParent(tree, i)]
          node_index[i] <- i
        }
      }
      shoulder_data <-
        data.frame(node = node_index,
                   x_anc = x_anc,
                   y = y)
      if (timeline == TRUE) {
        shoulder_data$x_anc <- shoulder_data$x_anc - tree_height
      }
      `%<+%` <- ggtree::`%<+%`
      p <- p %<+% shoulder_data
    }
    
    if (node_labels_as == "state") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = 1,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = start_state_1,
              x = x_anc,
              y = y
            ),
            hjust = 0,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = anc_state_1, subset = !isTip),
            hjust = 1,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(label = end_state_1, subset = !isTip),
            hjust = 1,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "state_posterior") {
      if (cladogenetic == TRUE) {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = 1,
            nudge_x = node_labels_offset,
            size = node_labels_size
          ) +
          ggtree::geom_text(
            ggplot2::aes(
              label = .convertAndRound(start_state_1_pp),
              x = x_anc,
              y = y
            ),
            hjust = 0,
            nudge_x = node_labels_offset,
            size = node_labels_size,
            na.rm = TRUE
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base[1] == "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(anc_state_1_pp),
              subset = !isTip
            ),
            hjust = 1,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      } else if (cladogenetic == FALSE &
                 state_pos_str_base != "anc_state_") {
        p <-
          p + ggtree::geom_text2(
            ggplot2::aes(
              label = .convertAndRound(end_state_1_pp),
              subset = !isTip
            ),
            hjust = 1,
            nudge_x = node_labels_offset,
            size = node_labels_size
          )
      }
    } else if (node_labels_as == "node_posterior") {
      p <-
        p + ggtree::geom_nodelab(
          ggplot2::aes(label = .convertAndRound(posterior)),
          hjust = 1,
          nudge_x = node_labels_offset,
          size = node_labels_size
        )
    }
  }
  
  # add tip states labels (text)
  if (tip_labels_states == TRUE) {
    if (state_pos_str_base[1] == "anc_state_") {
      p <- p + ggtree::geom_tiplab(
        ggplot2::aes(label = anc_state_1),
        hjust = "center",
        offset = tip_labels_states_offset,
        size = tip_labels_states_size
      )
    } else {
      p <- p + ggtree::geom_tiplab(
        ggplot2::aes(label = end_state_1),
        hjust = "center",
        offset = tip_labels_states_offset,
        size = tip_labels_states_size
      )
    }
  }
  
  # add space on x axis for tip labels
  if (tip_labels == TRUE & timeline == FALSE) {
    p <- p + ggtree::xlim(0, tree_height + tree_height / 2)
  }
  
  # clean up pngs
  unlink(paste0(tmp,"/.temp.png"))
  
  return(p)
}

plotStochMaps <- function(tree,
                          maps,
                          colors = "default",
                          color_by = "prob",
                          tree_layout = "rectangular",
                          line_width = 1,
                          tip_labels = TRUE,
                          tip_labels_italics = FALSE,
                          tip_labels_formatted = FALSE,
                          tip_labels_remove_underscore = TRUE,
                          tip_labels_color = "black",
                          tip_labels_size = 3,
                          tip_labels_offset = 0,
                          timeline = FALSE,
                          geo_units = list("epochs", "periods"),
                          geo = timeline,
                          time_bars = timeline,
                          label_sampled_ancs = FALSE,
                          ...) {
  
  # pull tree from list object if necessary
  if (inherits(tree,"list")) {
    if (length(tree) == 1){
      tree <- tree[[1]]
    } else {stop("tree should contain only one tree object")}
  }
  
  if (inherits(tree,"list")) {
    if (length(tree) == 1){
      tree <- tree[[1]]
    } else {stop("tree should contain only one tree object")}
  } 
  
  p <-  plotTreeFull(
    tree = list(list(tree)),
    tree_layout = tree_layout,
    line_width = line_width,
    
    tip_labels = tip_labels,
    tip_labels_italics = tip_labels_italics,
    tip_labels_formatted = tip_labels_formatted,
    tip_labels_remove_underscore = tip_labels_remove_underscore,
    tip_labels_color = tip_labels_color,
    tip_labels_size = tip_labels_size,
    tip_labels_offset = tip_labels_offset,
    
    timeline = timeline,
    geo_units = geo_units,
    geo = timeline,
    time_bars = timeline,
    
    label_sampled_ancs = label_sampled_ancs,
    
    node_age_bars = FALSE,
    age_bars_color = "blue",
    age_bars_colored_by = NULL,
    
    node_labels = NULL,
    node_labels_color = "black",
    node_labels_size = 3,
    node_labels_offset = 0,
    
    node_pp = FALSE,
    node_pp_shape = 16,
    node_pp_color = "black",
    node_pp_size = "variable",
    
    branch_color = "black",
    color_branch_by = NULL,
    
    tip_age_bars = FALSE,
    lineend = "square",
    ...
  )
  
  if (colors[1] != "default") {
    # error checking
    if (is.null(names(colors))) 
    {stop("colors must be a NAMED vector of colors where names correspond to the character states")}
    states <- names(colors)
  } else {
    states <- colnames(maps)[-c(1:5)]
    colors <- colFun(length(states))
    names(colors) <- states
  }
  
  dat <- dplyr::left_join(maps, p$data, by = "node")
  
  #set up colors 
  if (color_by == "MAP") {
    max <- apply(dat[, states], MARGIN = 1, which.max)
    seg_col <- colors[unlist(max)]
    dat$seg_col <- seg_col
    names(seg_col) <- seg_col
  } else if (color_by == "prob") {
    rgbcols <- col2rgb(colors)
    rgb_values_per_seg <- t(rgbcols %*% t(dat[,states]))
    seg_col <- tolower(grDevices::rgb(red   = rgb_values_per_seg[ ,1],
                                      green = rgb_values_per_seg[ ,2],
                                      blue  = rgb_values_per_seg[ ,3],
                                      maxColorValue = 255))
    dat$seg_col <- seg_col
    names(seg_col) <- seg_col
  }
  
  # horizontal segments
  dat_horiz <- dat[dat$vert == FALSE,]
  
  seg_horiz <- data.frame(
    x    = dat_horiz$x - dat_horiz$x0,
    xend = dat_horiz$x - dat_horiz$x1,
    y    = dat_horiz$y,
    yend = dat_horiz$y,
    col  = dat_horiz$seg_col
  )
  
  #vertical segments
  dat_vert <- dat[dat$vert == TRUE,]
  
  m <- match(x = dat_vert$parent, dat_vert$node)
  dat_vert$y_parent <- dat_vert[m, "y"]
  dat_vert$x_parent <- dat_vert[m, "x"]
  
  seg_vert <- data.frame(
    x = dat_vert$x_parent,
    xend = dat_vert$x_parent,
    y = dat_vert$y,
    yend = dat_vert$y_parent,
    col = dat_vert$seg_col
  )
  
  # plot! 
  
  p + ggplot2::geom_segment(
    data = seg_horiz,
    ggplot2::aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      color = col
    ),
    lineend = "square",
    size = line_width,
  ) +
    ggplot2::geom_segment(
      data = seg_vert,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        color = col
      ),
      lineend = "square",
      size = line_width, 
      
    ) +
    ggplot2::scale_color_manual(values = seg_col, 
                                breaks = colors,
                                name = "State",
                                labels = names(colors))
  
}


processStochMaps <- function(tree,
                             paths = NULL,
                             simmap = NULL,
                             states,
                             num_intervals = 1000,
                             verbose = TRUE,
                             ...) {
  
  # pull tree from list object if necessary
  if (inherits(tree,"list")) {
    if (length(tree) == 1){
      tree <- tree[[1]]
    } else {stop("tree should contain only one tree object")}
  }
  
  if (inherits(tree,"list")) {
    if (length(tree) == 1){
      tree <- tree[[1]]
    } else {stop("tree should contain only one tree object")}
  } 
  
  # compute the number of states
  nstates <- length(states)
  
  # create the index map
  map <- matchNodes(tree@phylo)
  
  # either paths or simmap must be provided
  if ( !is.null(paths) ) { # samples in files
    
    # read traces
    samples <- readTrace(paths, verbose = verbose, ...)
    
    # combine multiple samples together
    if ( length(samples) > 1 ) {
      samples <- combineTraces(samples)
      samples <- samples[[1]]
    } else {
      samples <- samples[[1]]
    }
    
    # compute the number of samples
    nsamples <- nrow(samples)
    
  } else if ( !is.null(simmap) ) { # samples in phytools format
    
    message("Reformatting simmap objects")
    
    # make the samples
    samples <- as.data.frame(do.call(rbind, lapply(simmap, function(simmap) {
      sapply(simmap$maps, function(edge) {
        edge <- rev(edge)
        return(paste0("{", paste0(paste0(names(edge),",", edge), collapse = ":"),"}"))
      })
    })))
    
    # add a root edge
    root_edge_samples <- sapply(simmap, function(map) {
      return(paste0("{", tail(names(map$maps[[1]]), n = 1), ",0}"))
    })
    samples <- cbind(samples, root_edge_samples)
    
    # get the nodes
    nodes <- c(tree@phylo$edge[,2], ape::Ntip(tree@phylo) + 1)
    colnames(samples) <- map$Rev[nodes]
    
    # compute the number of samples
    nsamples <- length(simmap)
    
  } else {
    stop("Please provide either a paths or simmap argument.")
  }
  
  message("Processing stochastic maps")
  
  # get the number of branches
  # including the root branch
  num_branches <- length(tree@phylo$edge.length) + 1
  root_index   <- ape::Ntip(tree@phylo) + 1
  
  # get the dt
  root_age <- max(ape::branching.times(tree@phylo))
  if (!is.null(tree@phylo$root.edge)) {
    root_age <- root_age + tree@phylo$root.edge
  } else {
    tree@phylo$root.edge <- 0
  }
  dt <- root_age / num_intervals
  
  # loop over branches
  dfs <- vector("list", num_branches)
  
  if (verbose) { pb <- txtProgressBar(min = 0, max = num_branches, initial = 0) }
  
  for(i in 1:num_branches) {
    
    # get the branch indexes
    R_index   <- map$R[i]
    Rev_index <- as.character(map[R_index,2])
    
    # get the time points
    if ( R_index == root_index ) {
      this_edge_length <- tree@phylo$root.edge
    } else {
      this_edge_length <- tree@phylo$edge.length[tree@phylo$edge[,2] == R_index]
    }
    these_pts <- seq(0, this_edge_length, by = dt)
    
    # get the samples
    branch_samples <- samples[,Rev_index]
    branch_samples <- gsub("{", "", branch_samples, fixed = TRUE)
    branch_samples <- gsub("}", "", branch_samples, fixed = TRUE)
    
    # split the per event
    branch_samples <- strsplit(branch_samples, ":")
    
    # get the times per state
    branch_samples <- lapply(branch_samples, function(sample) {
      sample <- do.call(rbind, strsplit(sample, ","))
      sample_states <- sample[,1]
      sample_times  <- as.numeric(sample[,2])
      names(sample_times) <- sample_states
      # sample_times <- rev(sample_times) # turn this on for plotting the output from (old) tensorphylo runs
      return(cumsum(sample_times))
    })
    
    # get the state per interval
    if ( this_edge_length == 0 ) {
      branch_states_per_interval <- t(t(match(names(unlist(branch_samples)), states)))
    } else {
      branch_states_per_interval <- do.call(rbind, lapply(branch_samples, function(sample) {
        match(names(sample)[findInterval(these_pts, sample) + 1], states)
      }))
    }
    
    # compute probability of each state per interval
    branch_prob_per_state <- apply(branch_states_per_interval, 2, tabulate, nbins = nstates) / nsamples
    rownames(branch_prob_per_state) <- states
    
    # now do the vertical segments
    vert_prob_per_state <- t(branch_prob_per_state[,ncol(branch_prob_per_state), drop = FALSE])
    
    # make the df
    this_df <- data.frame(index = Rev_index, bl = this_edge_length, x0 = these_pts, x1 = c(these_pts[-1], this_edge_length), vert = FALSE)
    this_df <- cbind(this_df, t(branch_prob_per_state))
    vert_df <- cbind(data.frame(index = Rev_index, bl = this_edge_length, x0 = this_edge_length, x1 = this_edge_length, vert = TRUE), vert_prob_per_state)
    this_df <- rbind(this_df, vert_df)        
    
    # store
    dfs[[i]] <- this_df
    
    if (verbose) { setTxtProgressBar(pb,i) }
    
  }
  if (verbose) { close(pb) }
  
  # combine the branches
  dfs <- do.call(rbind, dfs)
  
  # get node instead of index 
  # node is R's standard numbering for nodes
  # index is RevBayes specific 
  colnames(map) <- c("node", "index")
  map$index  <- as.character(map$index)
  dfs <- dplyr::full_join(map,dfs, by = "index")
  dfs$index <- NULL
  
  return(dfs)
  
}

plotTreeFull <- function(tree,
                         
                         timeline,
                         geo,
                         geo_units,
                         time_bars,
                         
                         node_age_bars,
                         tip_age_bars,
                         age_bars_color,
                         age_bars_colored_by,
                         
                         node_labels,
                         node_labels_color,
                         node_labels_size,
                         node_labels_offset,
                         
                         tip_labels,
                         tip_labels_italics,
                         tip_labels_formatted,
                         tip_labels_remove_underscore,
                         tip_labels_color,
                         tip_labels_size,
                         tip_labels_offset,
                         
                         label_sampled_ancs,
                         
                         node_pp,
                         node_pp_shape,
                         node_pp_color,
                         node_pp_size,
                         
                         branch_color,
                         color_branch_by,
                         line_width,
                         
                         tree_layout,
                         ...) {
  # enforce argument matching
  if (!is.list(tree))
    stop("tree should be a list of lists of treedata objects")
  if (!methods::is(tree[[1]][[1]], "treedata"))
    stop("tree should be a list of lists of treedata objects")
  vars <- colnames(tree[[1]][[1]]@data)
  if (is.logical(timeline) == FALSE)
    stop("timeline should be TRUE or FALSE")
  if (is.logical(node_age_bars) == FALSE)
    stop("node_age_bars should be TRUE or FALSE")
  if (any(.isColor(age_bars_color) == FALSE))
    stop("age_bars_color should be valid color(s)")
  if (is.null(age_bars_colored_by) == FALSE &
      any(vars %in% age_bars_colored_by) == FALSE)
    stop("age_bars_colored_by should be a column in your tidytree object")
  if (is.null(node_labels) == FALSE &
      any(vars %in% node_labels) == FALSE)
    stop("node_labels should be NULL or a column in your tidytree object")
  if (is.null(node_labels_color) == FALSE &
      .isColor(node_labels_color) == FALSE)
    stop("node_labels_color should be NULL or a recognized color")
  if (is.logical(tip_labels) == FALSE)
    stop("tip_labels should be TRUE or FALSE")
  if (is.logical(tip_labels_italics) == FALSE)
    stop("tip_labels_italics should be TRUE or FALSE")
  if (is.logical(tip_labels_formatted) == FALSE)
    stop("tip_labels_formatted should be TRUE or FALSE")
  if (tip_labels_italics == TRUE & tip_labels_formatted == TRUE) 
    stop("tip_labels_italics and tip_labels_formatted may not both be TRUE")
  if (.isColor(tip_labels_color) == FALSE)
    stop("tip_labels_color should be a recognized color")
  if (!methods::is(node_pp,"logical"))
    stop("node_pp should be TRUE or FALSE")
  if (node_pp) {
    if (length(node_pp_color) > 2)
      stop("node_pp_color should be of length 1 or 2")
    if (.isColor(node_pp_color) == FALSE)
      stop("node_pp_color should be a recognized color")
    if (node_pp_shape %in% 0:25 == FALSE)
      stop("node_pp_shape should be a recognized shape
           (value between 0 and 25)")
    if (is.numeric(node_pp_size) == FALSE &
        node_pp_size != "variable")
      stop("node_pp_size should be numeric or 'variable'")
  }
  if (is.logical(tip_age_bars) == FALSE)
    stop("tip_age_bars should be TRUE or FALSE")
  if (length(branch_color) == 1 &
      !.isColor(branch_color))
    stop("branch_color should be a recognized color")
  if (length(branch_color) == 2) {
    if (.isColor(branch_color[1] == FALSE) &
        .isColor(branch_color[2]) == FALSE)
      stop("Neither values of branch_color are a recognized color")
    if (.isColor(branch_color[1] == FALSE) &
        .isColor(branch_color[2]))
      stop("branch_color[1] is not a recognized color")
    if (.isColor(branch_color[1]) &
        .isColor(branch_color[2]) == FALSE)
      stop("branch_color[2] is not a recognized color")
  } else if (length(branch_color) > 2)
    stop("only 2 colors may be specified in branch_color")
  if (is.null(color_branch_by) == FALSE &
      any(vars %in% color_branch_by) == FALSE)
    stop("color_branch_by should be NULL or a column in your tidytree object")
  if (is.numeric(line_width) == FALSE)
    stop ("line_width should be numeric")
  if (is.logical(label_sampled_ancs) == FALSE)
    stop("label_sampled_ancs should be TRUE or FALSE")
  if (is.list(geo_units)) {
    if (length(geo_units) != 2)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[1]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
    if (geo_units[[2]] %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras')  == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  } else {
    if (geo_units %in% 
        c('periods', 'epochs', 'stages', 'eons', 'eras') == FALSE)
      stop(
        "geo_units should be 'periods', 'epochs', 'stages', 'eons', 
         'eras', or a list of two of those units, such as:
        list('epochs','periods')"
      )
  }
  if (is.numeric(tip_labels_offset) == FALSE)
    stop ("tip_labels_offset should be a number")
  if (is.numeric(node_labels_offset) == FALSE)
    stop ("node_labels_offset should be a number")
  tree_layout <-
    match.arg(
      tree_layout,
      choices = c(
        'rectangular',
        'slanted',
        'ellipse',
        'cladogram',
        'roundrect',
        'fan',
        'circular',
        'inward_circular',
        'radial',
        'equal_angle',
        'daylight',
        'ape'
      )
    )
  if (tree_layout != "rectangular") {
    if (timeline == TRUE) {
      stop("timeline is only compatible with
                                 tree_layout = 'rectangular'")
    }
  }
  # grab single tree from input
  phy <- tree[[1]][[1]]
  
  ### fix for trees with sampled ancestors ###
  phylo    <- phy@phylo
  node_map <- .matchNodesTreeData(phy, phylo)
  if ("sampled_ancestor" %in% colnames(phy@data)) {
    phy@data$node <-
      as.character(node_map[match(as.numeric(phy@data$index),
                                  node_map$Rev), ]$R)
  }
  
  ### set up tree layout ###
  
  if (tree_layout == "cladogram") {
    tree_layout <- "rectangular"
    BL <- "none"
  } else {
    BL <- "branch.length"
  }
  
  # initiate plot
  if (is.null(color_branch_by)) {
    pp <- ggtree::ggtree(
      phy,
      right = FALSE,
      size = line_width,
      color = branch_color,
      branch.length = BL,
      layout = tree_layout,
      ...
    )
  } else if (!is.null(color_branch_by)) {
    pp <- ggtree::ggtree(
      phy,
      right = FALSE,
      size = line_width,
      branch.length = BL,
      layout = tree_layout,
      ...
    )
  }
  
  #### parameter compatibility checks ###
  if (length(node_pp_color) == 2 &
      length(branch_color) == 2)
    stop(
      "You may only include variable colors for either node_pp_label or
      branch_color, not for both"
    )
  
  #check that if user wants node_age_bars, there are dated intervals in  file
  if (node_age_bars == TRUE) {
    if (!"age_0.95_HPD" %in% colnames(phy@data))
      stop(
        "You specified node_age_bars, but there is no age_0.95_HPD column
        in the treedata object."
      )
  }
  
  # get dimensions
  n_node <- treeio::Nnode(phy)
  tree_height <- max(phytools::nodeHeights(phy@phylo))
  ntips <- sum(pp$data$isTip)
  
  # reformat labels if necessary
  if (tip_labels_remove_underscore) {
    pp$data$label <- gsub("_", " ", pp$data$label)
  }
  
  # check that if user wants to label sampled ancs,
  # there are sampled ancs in the files
  if (label_sampled_ancs == TRUE &
      "sampled_ancestor" %in% colnames(pp$data)) {
    sampled_ancs <-
      pp$data[!pp$data$isTip & !is.na(pp$data$sampled_ancestor),]
    if (nrow(sampled_ancs) < 1) {
      label_sampled_acs <- FALSE
    }
  }
  
  # add timeline
  if (timeline == TRUE) {
    warning("Plotting with default axis label (Age (Ma))")
    if (node_age_bars == FALSE) {
      minmax <- phytools::nodeHeights(phy@phylo)
      max_age <- tree_height
    } else {
      pp$data$age_0.95_HPD <- lapply(pp$data$age_0.95_HPD, function(z) {
        if (any(is.null(z)) ||
            any(is.na(z))) {
          return(c(NA, NA))
        } else {
          return(as.numeric(z))
        }
      })
      minmax <- t(matrix(unlist(pp$data$age_0.95_HPD), nrow = 2))
      max_age <- max(minmax, na.rm = TRUE)
    }
    
    interval <- max_age / 5
    dx <- max_age %% interval
    
    # add geo
    tick_height <- ntips / 100
    if (geo == TRUE) {
      #determine whether to include quaternary
      if (tree_height > 50) {
        skipit <- c("Quaternary", "Holocene", "Late Pleistocene")
      } else {
        skipit <- c("Holocene", "Late Pleistocene")
      }
      # add deep timescale
      if (length(geo_units) == 1) {
        pp <- pp + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-max(minmax, na.rm = TRUE), tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          height = grid::unit(4, "line"),
          skip = skipit,
          abbrv = FALSE,
          rot = 90,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      } else if (length(geo_units) == 2) {
        pp <- pp + deeptime::coord_geo(
          dat  = geo_units,
          pos  = lapply(seq_len(length(geo_units)), function(x)
            "bottom"),
          size = lapply(seq_len(length(geo_units)), function(x)
            tip_labels_size),
          xlim = c(-max(minmax, na.rm = TRUE), tree_height /
                     2),
          ylim = c(-tick_height * 5, ntips *
                     1.1),
          skip = skipit,
          center_end_labels = TRUE,
          bord = c("right", "top", "bottom"),
          neg  = TRUE
        )
      }
    }
    #add axis title
    pp <- pp + ggplot2::scale_x_continuous(
      name = "Age (Ma)",
      expand = c(0, 0),
      limits = c(-max_age, tree_height /
                   2),
      breaks = -rev(seq(0, max_age +
                          dx, interval)),
      labels = rev(seq(0, max_age +
                         dx, interval))
    )
    pp <- ggtree::revts(pp)
    
    # add ma ticks and labels
    xline <- pretty(c(0, max_age))[pretty(c(0, max_age)) < max_age]
    df <-
      data.frame(
        x = -xline,
        y = rep(-tick_height * 5, length(xline)),
        vx = -xline,
        vy = rep(-tick_height * 5 + tick_height, length(xline))
      )
    
    pp <-
      pp + ggplot2::geom_segment(ggplot2::aes(
        x = 0,
        y = -tick_height * 5,
        xend = -max_age,
        yend = -tick_height * 5
      )) +
      ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      )) +
      ggplot2::annotate(
        "text",
        x = -rev(xline),
        y = -tick_height * 5 + tick_height * 2,
        label = rev(xline),
        size = tip_labels_size
      )
    
    # add vertical gray bars
    if (time_bars) {
      if (geo) {
        if ("epochs" %in% geo_units) {
          x_pos <- -rev(c(0, deeptime::getScaleData("epochs")$max_age))
        } else {
          x_pos <-  -rev(c(0, deeptime::getScaleData("periods")$max_age))
        }
      } else if (!geo) {
        x_pos <- -rev(xline)
      }
      for (k in 2:(length(x_pos))) {
        box_col <- "gray92"
        if (k %% 2 == 1)
          box_col <- "white"
        box <-
          ggplot2::geom_rect(
            xmin = x_pos[k - 1],
            xmax = x_pos[k],
            ymin = -tick_height * 5,
            ymax = ntips,
            fill = box_col
          )
        pp <- gginnards::append_layers(pp, box, position = "bottom")
      }
    }
    if (tip_labels) {
      # recenter legend
      tot <- max_age + tree_height / 2
      pp <-
        pp +
        ggplot2::theme(axis.title.x =
                         ggplot2::element_text(hjust = max_age / (2 * tot)))
    }
  }
  
  # processing for node_age_bars and tip_age_bars
  if (node_age_bars == TRUE) {
    # Encountered problems with using geom_range to plot age HPDs in ggtree. It
    # appears that geom_range incorrectly rotates the HPD relative to the height
    # of the node unnecessarily. My guess for this would be because older
    # version of ggtree primarily supported length measurements, and not
    # height measurements so the new capability to handle height might contain
    # a "reflection" bug.
    # For example, suppose a node has height 3 with HPD [2, 7]. You can think of
    # this offset as h + [2-h, 7-h]. ggtree seems to "rotate" this
    # causing the HPD to appear as [-1, 4]. Figtree displays this correctly.
    #
    # See this excellent trick by Tauana:
    # https://groups.google.com/forum/#!msg/bioc-ggtree/wuAlY9phL9Q/L7efezPgDAAJ
    # Adapted this code to also plot fossil tip uncertainty in red
    pp$data$age_0.95_HPD <-
      lapply(pp$data$age_0.95_HPD, function(z) {
        if (any(is.null(z)) ||
            any(is.na(z))) {
          return(c(NA, NA))
        } else {
          return(as.numeric(z))
        }
      })
    
    minmax <- t(matrix(unlist(pp$data$age_0.95_HPD), nrow = 2))
    bar_df <-
      data.frame(
        node_id = as.integer(pp$data$node),
        isTip = pp$data$isTip,
        as.data.frame(minmax)
      )
    names(bar_df) <- c("node_id", "isTip", "min", "max")
    if (tip_age_bars == TRUE) {
      tip_df <- dplyr::filter(bar_df, isTip == TRUE & !is.na(min))
      tip_df <-
        dplyr::left_join(tip_df, pp$data, by = c("node_id" = "node"))
      tip_df <- dplyr::select(tip_df, node_id, min, max, y)
    }
    if (is.null(age_bars_colored_by) == TRUE) {
      # plot age densities
      bar_df <-
        dplyr::left_join(bar_df, pp$data, by = c("node_id" = "node"))
      bar_df <- dplyr::select(bar_df,  node_id, min, max, y)
      pp <-
        pp + ggplot2::geom_segment(
          ggplot2::aes(
            x = -min,
            y = y,
            xend = -max,
            yend = y
          ),
          data = bar_df,
          color = age_bars_color,
          size = 1.5,
          alpha = 0.8
        )
    } else if (is.null(age_bars_colored_by) == FALSE) {
      if (length(age_bars_color) == 1) {
        age_bars_color <- colFun(2)[2:1]
      }
      
      if ("sampled_ancestor" %in% colnames(pp$data) == TRUE) {
        sampled_tip_probs <-
          1 - as.numeric(pp$data$sampled_ancestor[pp$data$isTip == TRUE])
        sampled_tip_probs[is.na(sampled_tip_probs)] <- 0
      } else {
        sampled_tip_probs <- rep(1, sum(pp$data$isTip))
      }
      
      pp$data$olena <-
        c(sampled_tip_probs,
          as.numeric(.convertAndRound(L =
                                        unlist(pp$data[pp$data$isTip == FALSE,
                                                       age_bars_colored_by]))))
      
      bar_df <-
        dplyr::left_join(bar_df, pp$data, by = c("node_id" = "node"))
      bar_df <-
        dplyr::select(bar_df,  node_id, min, max, y, olena, isTip = isTip.x)
      if (tip_age_bars == FALSE) {
        bar_df <- dplyr::filter(bar_df, isTip == FALSE)
      }
      pp <-
        pp + ggplot2::geom_segment(
          ggplot2::aes(
            x = -min,
            y = y,
            xend = -max,
            yend = y,
            color = olena
          ),
          data = bar_df,
          size = 1.5,
          alpha = 0.8
        ) +
        ggplot2::scale_color_gradient(
          low = age_bars_color[1],
          high = age_bars_color[2],
          name = .titleFormat(age_bars_colored_by),
          breaks = pretty(pp$data$olena)
        )
    }
  }
  
  # label sampled ancestors
  if (label_sampled_ancs == TRUE &
      "sampled_ancestor" %in% colnames(pp$data)) {
    sampled_ancs <-
      pp$data[!pp$data$isTip & !is.na(pp$data$sampled_ancestor),]
    space_labels <- ntips / 30
    if (tip_labels_italics == TRUE) {
      pp <-
        pp + ggplot2::annotate(
          "text",
          x = sampled_ancs$x,
          y = sampled_ancs$y,
          label = seq_len(nrow(sampled_ancs)),
          vjust = -.5,
          size = tip_labels_size,
          color = tip_labels_color
        ) +
        ggplot2::annotate(
          "text",
          x = rep(-max(
            unlist(pp$data$age_0.95_HPD), na.rm = TRUE
          ),
          times = nrow(sampled_ancs)),
          y = seq(
            from = ntips - space_labels,
            by = -space_labels,
            length.out = nrow(sampled_ancs)
          ),
          label = paste0(seq_len(nrow(sampled_ancs)),
                         ": italic(`",
                         sampled_ancs$label,
                         "`)"),
          size = tip_labels_size,
          color = tip_labels_color,
          hjust = 0,
          parse = TRUE
        )
    } else if (tip_labels_italics == TRUE) {
      pp <-
        pp + ggplot2::annotate(
          "text",
          x = sampled_ancs$x,
          y = sampled_ancs$y,
          label = seq_len(nrow(sampled_ancs)),
          vjust = -.5,
          size = tip_labels_size,
          color = tip_labels_color
        ) +
        ggplot2::annotate(
          "text",
          x = rep(-max(
            unlist(pp$data$age_0.95_HPD), na.rm = TRUE
          ),
          times = nrow(sampled_ancs)),
          y = seq(
            from = ntips - space_labels,
            by = -space_labels,
            length.out = nrow(sampled_ancs)
          ),
          label = paste0(seq_len(nrow(sampled_ancs)),
                         ": ",
                         sampled_ancs$label),
          size = tip_labels_size,
          color = tip_labels_color,
          hjust = 0,
          parse = TRUE
        )
    } else {
      pp <-
        pp + ggplot2::annotate(
          "text",
          x = sampled_ancs$x,
          y = sampled_ancs$y,
          label = seq_len(nrow(sampled_ancs)),
          vjust = -.5,
          size = tip_labels_size,
          color = tip_labels_color
        ) +
        ggplot2::annotate(
          "text",
          x = rep(-max(
            unlist(pp$data$age_0.95_HPD), na.rm = TRUE
          ),
          times = nrow(sampled_ancs)),
          y = seq(
            from = ntips - space_labels,
            by = -space_labels,
            length.out = nrow(sampled_ancs)
          ),
          label = paste0(seq_len(nrow(sampled_ancs)),
                         ": ",
                         sampled_ancs$label),
          size = tip_labels_size,
          color = tip_labels_color,
          hjust = 0
        )
    }
    
    t_height <- ntips / 200
    df <- data.frame(
      x = sampled_ancs$x,
      vx = sampled_ancs$x,
      y = sampled_ancs$y + t_height,
      vy = sampled_ancs$y - t_height
    )
    pp <-
      pp + ggplot2::geom_segment(data = df, ggplot2::aes(
        x = x,
        y = y,
        xend = vx,
        yend = vy
      ))
    
  }
  
  # add node labels (text)
  if (is.null(node_labels) == FALSE) {
    # catch some funkiness from importing an unrooted tree
    if (node_labels == "posterior") {
      pp$data[grep("]:", unlist(pp$data[, node_labels])), node_labels] <-
        NA
    }
    pp$data$kula <-
      c(rep(NA, times = ntips),
        .convertAndRound(L =
                           unlist(pp$data[pp$data$isTip == FALSE,
                                          node_labels])))
    
    # change any NAs that got converted to characters back to NA
    pp$data$kula[pp$data$kula == "NA"] <- NA
    pp <- pp + ggtree::geom_text(
      ggplot2::aes(label = kula),
      color = node_labels_color,
      nudge_x = node_labels_offset,
      hjust = 0,
      size = node_labels_size
    )
  }
  
  # add tip labels (text)
  if (tip_labels == TRUE) {
    if (tip_age_bars == TRUE) {
      pp$data$extant <- !pp$data$node %in% tip_df$node_id
    } else {
      pp$data$extant <- TRUE
    }
    if (tip_labels_italics == TRUE) {
      pp <- pp + ggtree::geom_tiplab(
        ggplot2::aes(
          subset = extant & isTip,
          label = paste0('italic(`', label, '`)')
        ),
        size = tip_labels_size,
        offset = tip_labels_offset,
        color = tip_labels_color,
        parse = TRUE
      )
      if (tip_age_bars == TRUE) {
        new_tip_df <-
          dplyr::left_join(tip_df,
                           pp$data[, c("label", "node")],
                           by = c("node_id" = "node"))
        pp <-
          pp + ggplot2::annotate(
            "text",
            x = -new_tip_df$min + tip_labels_offset,
            y = new_tip_df$y,
            label = paste0('italic(`', new_tip_df$label, '`)'),
            hjust = 0,
            color = tip_labels_color,
            size = tip_labels_size,
            parse = TRUE
          )
      }
    } else if (tip_labels_formatted == TRUE ) {
      pp <- pp + ggtree::geom_tiplab(
        ggplot2::aes(subset = extant & isTip,
                     label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        color = tip_labels_color,
        parse = TRUE
      )
      if (tip_age_bars == TRUE) {
        new_tip_df <-
          dplyr::left_join(tip_df,
                           pp$data[, c("label", "node")],
                           by = c("node_id" = "node"))
        pp <-
          pp + ggplot2::annotate(
            "text",
            x = -new_tip_df$min + tip_labels_offset,
            y = new_tip_df$y,
            label = new_tip_df$label,
            hjust = 0,
            color = tip_labels_color,
            size = tip_labels_size,
            parse = TRUE 
          )
      }
    } else {
      pp <- pp + ggtree::geom_tiplab(
        ggplot2::aes(subset = extant & isTip,
                     label = label),
        size = tip_labels_size,
        offset = tip_labels_offset,
        color = tip_labels_color
      )
      if (tip_age_bars == TRUE) {
        new_tip_df <-
          dplyr::left_join(tip_df,
                           pp$data[, c("label", "node")],
                           by = c("node_id" = "node"))
        pp <-
          pp + ggplot2::annotate(
            "text",
            x = -new_tip_df$min + tip_labels_offset,
            y = new_tip_df$y,
            label = new_tip_df$label,
            hjust = 0,
            color = tip_labels_color,
            size = tip_labels_size
          )
      }
    }
  }
  
  # add node PP (symbols)
  if (node_pp == TRUE) {
    # reformat posterior
    pp$data$posterior <- as.numeric(pp$data$posterior)
    
    if (length(node_pp_color) == 1 & node_pp_size == "variable") {
      pp <- pp + ggtree::geom_nodepoint(color = node_pp_color,
                                        ggplot2::aes(size = posterior),
                                        shape = node_pp_shape) +
        ggplot2::scale_size_continuous(name = "Posterior")
    } else if (length(node_pp_color) == 2 &
               node_pp_size != "variable") {
      pp <- pp + ggtree::geom_nodepoint(size = node_pp_size,
                                        ggplot2::aes(color = posterior),
                                        shape = node_pp_shape) +
        ggplot2::scale_color_gradient(
          name = "Posterior",
          low = node_pp_color[1],
          high = node_pp_color[2],
          breaks = pretty(pp$data$posterior)
        )
    }
  }
  
  # add branch coloration by variable
  if (is.null(color_branch_by) == FALSE) {
    #set default colors if none provided
    if (length(branch_color) != 2) {
      branch_color <- c("#005ac8", "#fa7850")
    }
    col_num <- which(colnames(pp$data) == color_branch_by)
    pp$data[, col_num] <-
      as.numeric(as.data.frame(pp$data)[, col_num])
    name <- .titleFormat(color_branch_by)
    pp <- pp + 
      ggplot2::aes(color = as.data.frame(pp$data)[, col_num]) +
      ggplot2::scale_color_gradient(
        low = branch_color[1],
        high = branch_color[2],
        breaks = pretty(as.data.frame(pp$data)[, col_num]),
        name = name
      )
  }
  
  # readjust axis for non-timeline plots
  if (timeline == FALSE & BL != "none") {
    if (node_age_bars == FALSE) {
      xlim_min <- -tree_height
    } else {
      xlim_min <- -max(t(matrix(
        unlist(pp$data$age_0.95_HPD),
        nrow = 2
      )), na.rm = TRUE)
    }
    
    if (tip_labels == TRUE) {
      xlim_max <- tree_height / 2
    } else {
      xlim_max <- 0
    }
    
    pp <- pp + ggtree::xlim(xlim_min, xlim_max)
    pp <- ggtree::revts(pp)
    
  }
  
  # readjust axis for cladograms
  if (timeline == FALSE & BL == "none") {
    xlim_min <- range(pp$data$x)[1]
    
    if (tip_labels == TRUE) {
      xlim_max <- range(pp$data$x)[2] * 1.5
    } else {
      xlim_max <- range(pp$data$x)[2]
    }
    
    pp <- pp + ggtree::xlim(xlim_min, xlim_max)
    
  }
  
  return(pp)
}


# Non-exported utility functions for RevGadgets

# set custom state labels
.assign_state_labels <-
  function(t,
           state_labels,
           include_start_states,
           labels_as_numbers,
           missing_to_NA,
           n_states = 2) {
    # what is the ancestral state name tag?
    if (include_start_states) {
      state_pos_str_base <- c("start_state_", "end_state_")
    } else {
      state_pos_str_base <- c("anc_state_")
    }
    
    # send error if state_labels are provided without names
    if (!is.null(state_labels) && is.null(names(state_labels))) {
      stop(
        "names(state_labels) must identify all unlabeled state names in
        attributes(t)$data"
      )
    }
    
    # make matrix of all anc state values
    col_num <- grep(state_pos_str_base[1], colnames(t@data))
    if (length(state_pos_str_base) > 1) {
      col_num2 <- grep(state_pos_str_base[2], colnames(t@data))
      col_num <- c(col_num, col_num2)
    }
    pps <- grep("_pp", colnames(t@data))
    columns <- col_num[!col_num %in% pps]
    
    # change ? to NA
    if (missing_to_NA == TRUE) {
      for (c in columns) {
        x_state <- attributes(t)$data[[c]]
        x_state <- as.vector(x_state)
        x_state[x_state == "?"] <- "NA"
        attributes(t)$data[[c]] <- x_state
      }
    }
    
    all_anc_states <- unique(c(as.matrix(t@data[, columns])))
    
    # send error if state labels are provided but there are any
    # states without a corresponding state label
    if (!is.null(state_labels) &&
        any(all_anc_states %in% c("NA", names(state_labels)) == FALSE)) {
      stop(paste0(
        "names(state_labels): ",
        paste0(names(state_labels), collapse = ", "),
        " do not match data in tree file: ",
        paste0(sort(all_anc_states[all_anc_states != "NA"]), collapse = ", ")
      ))
    }
    
    # generate state labels if none provided and not a chromosome analysis
    if (is.null(state_labels) == TRUE & labels_as_numbers == FALSE) {
      warning("State labels not provided by user.
              Will be generated automatically.")
      states <-
        unique(unlist(attributes(t)$data[grepl(paste0("state_", "[0-9]$"),
                                               names(attributes(t)$data))]))
      states <- states[!states == "NA"]
      states <- states[order(states)]
      state_labels <- list()
      for (i in seq_len(length(states))) {
        state_labels[as.character(states[i])] <- LETTERS[i]
      }
      state_labels["other"] <- "other"
    }
    
    # for chromosome analyses, just keep the names as is (numbers of chromos)
    if (is.null(state_labels) == TRUE & labels_as_numbers == TRUE) {
      state_labels <-
        unique(unlist(attributes(t)$data[grepl(paste0("state_", "[0-9]$"),
                                               names(attributes(t)$data))]))
      state_labels <- state_labels[-which(state_labels == "NA")]
      names(state_labels) <- state_labels
    }
    
    # create list of ancestral state name tags
    state_pos_str_to_update <-
      c(unlist(lapply(1:n_states, function(x) {
        paste(state_pos_str_base, x, sep = "")
      })))
    
    # overwrite state labels
    for (m in state_pos_str_to_update) {
      # get the states
      x_state <- attributes(t)$data[[m]]
      x_state <- as.vector(x_state)
      x_state_valid <- which(x_state != "NA")
      x_state_invalid <- which(x_state == "NA")
      x_state_tmp <-
        unlist(lapply(x_state, function(z) {
          state_labels[names(state_labels) == z]
        }))
      x_state[x_state_valid] <- x_state_tmp
      x_state[x_state_invalid] <- NA
      if (labels_as_numbers) {
        x_state <-
          factor(x_state, levels = as.character(sort(as.integer(
            unique(state_labels)
          ))))
      }
      attributes(t)$data[[m]] <- x_state
    }
    
    # Just add the USED state_labels here
    used_state_labels <-
      na.omit(unique(c(as.matrix(t@data[, columns]))))
    if (labels_as_numbers) {
      attributes(t)$state_labels <-
        factor(used_state_labels, levels = as.character(sort(as.integer(
          unique(state_labels)
        ))))
    } else {
      attributes(t)$state_labels <- sort(as.character(used_state_labels))
    }
    
    return(t)
  }


.build_state_probs <-
  function(t,
           state_labels,
           include_start_states,
           p_threshold = 0) {
    n_states <- length(state_labels)
    n_tips <- length(attributes(t)$phylo$tip.label)
    n_node <- 2 * n_tips - 1
    
    dat <- list()
    
    if (include_start_states == TRUE) {
      state_tags <- c("end", "start")
    } else if (include_start_states == FALSE &
               "anc_state_1" %in% colnames(t@data)) {
      state_tags <- c("anc")
    } else if (include_start_states == FALSE &
               "end_state_1" %in% colnames(t@data)) {
      state_tags <- c("end")
    }
    
    for (s in state_tags) {
      dat[[s]] <- data.frame(matrix(0, nrow = n_node, ncol = n_states))
      #dat[[s]] = cbind(node=1:n_node, dat[[s]])
      
      for (i in 1:3)
      {
        m <- paste(s, "_state_", i, sep = "")
        pp_str <- paste(m, "_pp", sep = "")
        n_tmp <-
          as.numeric(as.vector(attributes(t)$data$node)) # node index
        x_tmp <- as.vector(attributes(t)$data[[m]])
        pp_tmp <- as.numeric(as.vector(attributes(t)$data[[pp_str]]))
        
        for (j in seq_len(length(x_tmp)))
        {
          if (!is.na(x_tmp[j])) {
            if (pp_tmp[j] > p_threshold) {
              k <- which(x_tmp[j] == state_labels)
              dat[[s]][n_tmp[j], k] <- pp_tmp[j]
            }
          }
        }
      }
      
      # format column names
      colnames(dat[[s]]) <- as.vector(unlist(state_labels))
      
      # add probs for >3rd state under "other" label
      rem_prob <- c()
      for (i in seq_len(nrow(dat[[s]]))) {
        rem_prob[i] <- 1
        for (j in seq_len(length(dat[[s]][i, ]))) {
          rem_prob[i] <- rem_prob[i] - dat[[s]][i, j]
        }
      }
      dat[[s]]$"other" <- rem_prob
      dat[[s]]$node <- 1:n_node
    }
    return(dat)
  }

.buildTranslateDictionary <- function(lines) {
  start_tree_block <- grep("begin trees;", lines, ignore.case = TRUE)
  end_tree_block <-
    grep("end;", lines[start_tree_block:length(lines)],
         ignore.case = TRUE)[1] + start_tree_block - 1
  tree_block <- lines[start_tree_block:end_tree_block]
  
  # look for translate block
  start_translations <-
    grep("translate", tree_block, ignore.case = TRUE)
  end_translations <-
    grep(";",
         tree_block[start_translations:length(tree_block)])[1] +
    start_translations -  1
  translations <- tree_block[start_translations:end_translations]
  
  # grab only the numbers and taxa names
  translations <- translations[grep("[1-9]", translations)]
  
  # remove commas
  translations <- gsub(",", "", translations)
  
  # replace tabs with space
  translations <- gsub("\\\t", " ", translations)
  
  # split at white space
  translations_split <- strsplit(translations, " ")
  
  # strip out empty elements
  translation_table <-
    do.call(rbind, lapply(translations_split, function(x)
      x[x != ""]))
  
  # create the data frame
  dictionary <-
    as.data.frame(translation_table, stringsAsFactors = FALSE)
  colnames(dictionary) <- c("number", "taxon")
  
  return(dictionary)
}

.collect_probable_states <- function(p, p_threshold = 0.005) {
  labels <- c("end_state", "start_state")
  index <- c(1, 2, 3)
  
  codes <- c()
  labels_pp <- c()
  for (l in labels) {
    for (i in index) {
      label_index <- paste(l, "_", i, sep = "")
      label_index_pp <- paste(l, "_", i, "_pp", sep = "")
      index_threshold <- p$data[[label_index_pp]] > p_threshold
      codes <-
        c(codes, unique(p$data[[label_index]][index_threshold]))
    }
  }
  codes <- unique(codes)
  codes <- c(codes, "other")
  return(codes)
}

.computeInterval <- function(item, rates, probs, summary = "mean") {
  interval_times <-
    unlist(rates[["speciation time"]][1, grepl("interval_times",
                                               names(rates$`speciation time`))])
  # For some reason these are ordered differently than rate vectors
  interval_times <-
    sort(interval_times)
  
  rate <- rates[[item]]
  rate <- rate[, grep("[0-9]", colnames(rate))]
  
  #mean_rate <- colMeans(rate)
  summary_rate <- apply(rate, 2, summary)
  quantiles <- apply(rate, 2,
                     quantile,
                     probs = probs)
  
  df <- dplyr::tibble(.rows = length(summary_rate))
  df["value"] <- summary_rate
  df["lower"] <- quantiles[1, ]
  df["upper"] <- quantiles[2, ]
  df$time <- interval_times
  df$item <- item
  
  return(df)
}

.convertAndRound <- function(L) {
  #sometimes there will be NAs before forcing to convert
  # got to remove nas before doing this test!
  k <- L[!is.na(L)]
  if (any(is.na(as.numeric(k))) == FALSE) {
    # if integer or numeric
    if (sum(as.numeric(L) %% 1, na.rm = TRUE) == 0) {
      # if integer
      labs <- L
      labs[labs == "1.000000"] <-
        "1" # catch case of all posterios of 1
    } else {
      # if numeric
      labs <- sprintf("%.3f", as.numeric(L)) # round nicely
      labs[labs == "1.000"] <- "1"
    }
  } else {
    # if character
    labs <- L
  }
  return(labs)
}

.findTreeLines <- function(lines) {
  # pull out tree block only
  start_tree_block <-
    grep("begin trees;", lines, ignore.case = TRUE)
  end_tree_block <-
    grep("end;",
         lines[start_tree_block:length(lines)],
         ignore.case = TRUE)[1] + start_tree_block - 1
  tree_block <- lines[start_tree_block:end_tree_block]
  
  # pull out trees
  
  # find all starting lines by searching for "tree"
  trees_start <- grep("tree ", tree_block, ignore.case = TRUE)
  # find all ending lines by searching for ";"
  # except for the last line of the tree block
  semicols <- grep("\\;", tree_block)
  semicols <- semicols[semicols >= trees_start[1]]
  trees_end <- semicols[1:(length(semicols) - 1)]
  # if tree are each on one line, return tree strings,
  # else concatenate multiple lines
  if (all(trees_start  == trees_end)) {
    tree_strings <-
      tree_block[grep("tree ", tree_block, ignore.case = TRUE)]
  } else {
    stop("RevGadgets currently doesn't support line breaks
         in trees in nexus files")
  }
  #  search  for  semicolon to signal end of line
  
  # return tree strings
  return(tree_strings)
  
}

# identify parent of a node
.getParent <- function(phylo, node) {
  if (.getRoot(phylo) == node) {
    return(0)
  } else {
    parent <- phylo$edge[phylo$edge[,2] == node,1]
    return(parent)
  }
}

# from ggtree::ggpie (called by .nodepie())
.ggpie <- function(data, y, fill, color, alpha=1, outline.color="transparent", outline.size=0) {
  p <- ggplot2::ggplot(data, ggplot2::aes_(x=1, y=y, fill=fill)) +
    ggplot2::geom_bar(stat='identity', alpha=alpha, color=outline.color, linewidth=outline.size, show.legend = F) +
    ggplot2::coord_polar(theta='y') + ggtree::theme_inset()
  
  if (methods::missingArg(color) == TRUE || is.null(color) == TRUE || any(is.na(color)) == TRUE) {
    ## do nothing
  } else {
    p <- p + ggplot2::scale_fill_manual(values=color)
  }
  return(p)
}

.isColor <- function(var) {
  if (is.null(var)) {
    return(FALSE)
  } else {
    t <- try(col2rgb(var), silent = TRUE)
    if (length(t) == 1 && methods::is(t, "try-error")) {
      return(FALSE)
    }
    else
      return(TRUE)
  }
}

.isNexusFile <- function(file)
  readLines(file, n = 1) == "#NEXUS"

.isSingleNewick <-
  function(file)
    strsplit(readLines(file, n = 1), split = "")[[1]][1] == "("

.makeNodeNames <- function(tree) {
  pr <- ape::prop.part(tree)
  labels <- attributes(pr)$labels
  names(labels) <- seq_len(length(labels))
  nodes <- lapply(pr[seq_len(length(pr))], dplyr::recode,!!!labels)
  nodes <- append(attributes(pr)$labels, nodes)
  
  node_names <- numeric()
  node_names_op <- numeric()
  for (i in seq_len(length(nodes))) {
    node_names[i] <-
      paste(as.numeric(sort(tree$tip.label) %in% nodes[[i]]),
            sep = "",
            collapse = "")
    node_names_op[i] <-
      paste(as.numeric(!sort(tree$tip.label) %in% nodes[[i]]),
            sep = "",
            collapse = "")
  }
  return(data.frame(node_names = node_names,
                    node_names_op = node_names_op))
}

.makePlotData <- function(rates, probs, summary) {
  rates <- .removeNull(rates)
  res <-
    lapply(names(rates), function(e)
      .computeInterval(
        e,
        rates = rates,
        probs = probs,
        summary = summary
      ))
  plotdata <- do.call(rbind, res)
  plotdata$item <- factor(
    plotdata$item,
    levels = c(
      "speciation rate",
      "extinction rate",
      "speciation time",
      "extinction time",
      "net-diversification rate",
      "relative-extinction rate"
    )
  )
  return(plotdata)
}

.makeStates <- function(label_fn, color_fn) {
  # generate colors for ranges
  range_color_list <-
    read.csv(color_fn,
             header = TRUE,
             sep = ",",
             colClasses = "character")
  
  # get area names
  area_names <-
    unlist(lapply(range_color_list$range, function(y) {
      if (nchar(y) == 1) {
        return(y)
      }
    }))
  
  # get state labels
  state_descriptions <-
    read.csv(label_fn,
             header = TRUE,
             sep = ",",
             colClasses = "character")
  
  # map presence-absence ranges to area names
  range_labels <-
    unlist(lapply(state_descriptions$range[2:nrow(state_descriptions)],
                  function(x) {
                    present <- as.vector(gregexpr(pattern = "1", x)[[1]])
                    paste(area_names[present], collapse = "")
                  }))
  
  # map labels to colors
  range_colors <-
    range_color_list$color[match(range_labels, range_color_list$range)]
  
  # generate state/color labels
  idx <- 1
  st_lbl <- list()
  st_colors <- c()
  for (j in 1:(nrow(state_descriptions) - 1)) {
    st_lbl[[as.character(j)]] <- range_labels[j]
    st_colors[j] <- range_colors[j]
  }
  st_colors[length(st_colors) + 1] <- "lightgray"
  st_lbl[["other"]] <- "other"
  
  return(list(state_labels = st_lbl, state_color = st_colors))
}

# from ggtree::nodepie, calls .ggpie
.nodepie <- function(data, cols, color, alpha=1, outline.color="transparent", outline.size=0) {
  if (! "node" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL
  if (methods::missingArg(color)) {
    color <- NA
  }
  `%>%` <- dplyr::`%>%`
  ldf <- tidyr::gather(data, type, value, !! cols) %>% split(., .$node)
  lapply(ldf, function(df) .ggpie(df, y=~value, fill=~type, color, alpha, outline.color, outline.size))
}


.parseTreeString <- function(string) {
  # recover()
  text <- sub("[^(]*", "", string)
  # stats <- treeio:::read.stats_beast_internal( "", text )
  
  # stats <- .read.stats_revbayes_internal("", text)
  # tree <- ape::read.tree(text = text)
  # obj <- .beast("", text, stats, tree)
  
  obj <- treeio::read.beast.newick(textConnection(text))
  
  if ("index" %in% colnames(obj@data)) {
    obj@data$index <- as.character(obj@data$index)
  } else {
    warning("No index found in tree file.
            This file may not work with downstream plotting functions.")
  }
  
  return(obj)
}

# Right tail probability of the horseshoe: integrates the density function
# via grid
.pRightTailHorseshoeGrid <- function(x,
                                     gamma = 1,
                                     grid_size = 5000) {
  quants <- seq(1e-10, 1 - 1e-10, length.out = grid_size)
  # Transform so we can look up quantiles under regular cauchy distribution
  quants <- 1.0 - (1.0 - quants) / 2.0
  probs <-
    1 / length(quants) # we're using quantiles, each gamma is equally likely
  sigmas <- qcauchy(quants, 0, gamma)
  sum(pnorm(x, 0, sigmas, lower.tail = FALSE) * probs)
}

# for reading in RevBayes output, especially if the number
# of elements per line varies
.readOutputFile <- function(path, burnin = 0.25) {
  
  `%>%` <- dplyr::`%>%`
  
  res <- path %>% 
    readLines() %>%
    utils::tail(n = -1)
  
  names <- path %>% 
    readLines() %>%
    utils::head(n = 1) %>%
    strsplit("\t")
  
  if (burnin >= length(res))
    stop("Burnin larger than provided trace file")
  
  if (burnin >= 1) {
    res <- res[(burnin + 1):length(res)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin * length(res))
    res <- res[(discard + 1):length(res)]
  } else if (burnin == 0) {
    res <- res
  } else {
    stop("What have you done?")
  }
  
  names_to_exclude = c("Iteration|Replicate_ID|Posterior|Likelihood|Prior")
  cols_to_exclude = length(grep(pattern = names_to_exclude, names[[1]]))
  
  res <- res %>%
    strsplit("\t")
  
  if (cols_to_exclude > 0) {
    res <- res %>%
      lapply(function(x) utils::tail(x, n = -cols_to_exclude))
  }
  
  res <- res %>%
    lapply(as.numeric)
  
  return(res)
}

.readNexusTrees <- function(path, burnin, verbose) {
  # read the lines
  lines <- readLines(path)
  
  # the line with a tree
  tree_strings <- .findTreeLines(lines)
  
  # discard burnin (if provided)
  if (burnin >= 1) {
    tree_strings <- tree_strings[(burnin + 1):length(tree_strings)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin * length(tree_strings))
    tree_strings <- tree_strings[(discard + 1):length(tree_strings)]
  } else if (burnin == 0) {
    tree_strings <- tree_strings
  } else {
    stop("What have you done?")
  }
  
  # get the trees
  n_trees <- length(tree_strings)
  if (verbose == TRUE) {
    bar <- txtProgressBar(style = 3, width = 40)
  }
  trees <- vector("list", n_trees)
  for (i in 1:n_trees) {
    trees[[i]] <- .parseTreeString(tree_strings[i])
    if (verbose == TRUE) {
      setTxtProgressBar(bar, i / n_trees)
    }
  }
  if (verbose == TRUE) {
    close(bar)
  }
  
  # translate using dictionary if translate block present in file
  if (length(grep("translate", lines, ignore.case = TRUE)) >= 1) {
    dictionary <- .buildTranslateDictionary(lines = lines)
    for (i in 1:n_trees) {
      n_tips <- length(trees[[i]]@phylo$tip.label)
      for (j in 1:n_tips) {
        ind <- which(trees[[i]]@phylo$tip.label[j] == dictionary[, 1])
        trees[[i]]@phylo$tip.label[j] <- dictionary[ind, 2]
      }
    }
  }
  
  # return the trees
  return(trees)
  
}

.readTreeLogs <- function(path, tree_name, burnin, verbose) {
  # read the samples
  samples <-
    utils::read.table(
      path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  
  # check that there is a column with the given name
  if (tree_name %in% colnames(samples) == FALSE) {
    stop(paste0("No column named ", tree_name, " found."))
  }
  
  # get the tree strings
  tree_strings <- samples[, tree_name]
  
  # discard burnin (if provided)
  if (burnin >= 1) {
    tree_strings <- tree_strings[(burnin + 1):length(tree_strings)]
  } else if (burnin < 1 & burnin > 0) {
    discard <- ceiling(burnin * length(tree_strings))
    tree_strings <- tree_strings[(discard + 1):length(tree_strings)]
  } else if (burnin == 0) {
    tree_strings <- tree_strings
  } else {
    stop("What have you done?")
  }
  
  # get the trees
  n_trees <- length(tree_strings)
  if (verbose == TRUE) {
    bar <- txtProgressBar(style = 3, width = 40)
  }
  trees <- vector("list", n_trees)
  for (i in 1:n_trees) {
    trees[[i]] <- .parseTreeString(tree_strings[i])
    if (verbose == TRUE) {
      setTxtProgressBar(bar, i / n_trees)
    }
  }
  if (verbose == TRUE) {
    close(bar)
  }
  
  # return the trees
  return(trees)
  
}

# functionally the same as .rootNode, but our own function
# (not modified from treeio)
# currently called by getParent()
.getRoot <- function(phylo) {
  edge1 <- phylo$edge[,1]
  edge2 <- phylo$edge[,2]
  root <- unique(edge1)[!unique(edge1) %in% unique(edge2)]
  return(root)
}

# Calculates global scale parameter for a Gaussian Markov random fielf from
# the prior mean number of "effective shifts" in the rate.
.setHSMRFGlobalScaleExpectedNumberOfJumps <-
  function(n_episodes,
           prior_n_shifts = log(2),
           shift_size = 2) {
    # We treat the change between each grid cell as a Bernoulli RV, so the
    # collection of changes becomes binomial
    # From this we can calculate the expected number of cells where a
    #shift occurs
    
    # Move to log-scale
    shift <- log(shift_size)
    
    # Probability of a shift for a value of zeta
    # We average the conditional p(shift | gamma) over p(gamma)
    quants <- seq(0.0001, 0.9999, length.out = 2000)
    
    # Transform so we can look up quantiles under regular cauchy distribution
    quants <- 1.0 - (1.0 - quants) / 2.0
    probs <-
      1 / length(quants) # we're using quantiles, each gamma is equally likely
    
    # Function to optimize
    fn <- function(zeta) {
      # Grid of gammas
      gammas <- qcauchy(quants, 0, zeta)
      # Number of expected shifts for each value of sigma
      num_expected_shifts <- unlist(lapply(gammas, function(x) {
        p_shift_one_cell_this_gamma <-
          .pRightTailHorseshoeGrid(shift, x, grid_size = 2000) / 0.5
        return(p_shift_one_cell_this_gamma * (n_episodes - 1))
      }))
      # Average the per-sigma E(n_shifts) over p(sigma) to get overall
      # expectation given zeta
      this_expected_num_shifts <- sum(probs * num_expected_shifts)
      # Distance to target
      return((log(this_expected_num_shifts) - log(prior_n_shifts)) ^ 2)
    }
    
    # Find best value of zeta
    opts <- optimize(fn, c(0, 1))
    zeta <- opts$minimum
    
    # Compute the prior on number of shifts for this zeta (to show user
    # how well we approximated the target)
    gammas <- qcauchy(quants, 0, zeta)
    num_expected_shifts <- unlist(lapply(gammas, function(x) {
      p_shift_one_cell_this_gamma <-
        .pRightTailHorseshoeGrid(shift, x, grid_size = 2000) / 0.5
      return(p_shift_one_cell_this_gamma * (n_episodes - 1))
    }))
    
    # Estimate the error of our chosen global scale hyperprior
    computed_num_expected_shifts <- sum(probs * num_expected_shifts)
    return(list(hyperprior = zeta, E.n = computed_num_expected_shifts))
  }

# Calculates global scale parameter for a Gaussian Markov random fielf
# from the prior mean number of "effective shifts" in the rate.
.setGMRFGlobalScaleExpectedNumberOfJumps <-
  function(n_episodes,
           prior_n_shifts = log(2),
           shift_size = 2) {
    # We treat the change between each grid cell as a Bernoulli RV, so the
    # collection of changes becomes binomial
    # From this we can calculate the expected number of cells where a shift
    # occurs
    
    # Move to log-scale
    shift <- log(shift_size)
    
    # Probability of a shift for a value of zeta
    # We average the conditional p(shift | sigma) over p(sigma)
    quants <- seq(0.0001, 0.9999, length.out = 2000)
    
    # Transform so we can look up quantiles under regular cauchy distribution
    quants <- 1.0 - (1.0 - quants) / 2.0
    probs <-
      1 / length(quants) # we're using quantiles, each gamma is equally likely
    
    # Function to optimize
    fn <- function(zeta) {
      # Grid of sigmas
      sigmas <- qcauchy(quants, 0, zeta)
      # Number of expected shifts for each value of sigma
      num_expected_shifts <- unlist(lapply(sigmas, function(x) {
        p_shift_one_cell_this_sigma <- pnorm(shift, 0, x, lower.tail = FALSE) /
          0.5
        return(p_shift_one_cell_this_sigma * (n_episodes - 1))
      }))
      # Average the per-sigma E(n_shifts) over p(sigma) to get overall
      # expectation given zeta
      this_expected_num_shifts <- sum(probs * num_expected_shifts)
      # Distance to target
      return((log(this_expected_num_shifts) - log(prior_n_shifts)) ^ 2)
    }
    
    # Find best value of zeta
    opts <- optimize(fn, c(0, 1))
    zeta <- opts$minimum
    
    # Compute the prior on number of shifts for this zeta (to show user how
    # well we approximated the target)
    sigmas <- qcauchy(quants, 0, zeta)
    num_expected_shifts <- unlist(lapply(sigmas, function(x) {
      p_shift_one_cell_this_sigma <- pnorm(shift, 0, x, lower.tail = FALSE) /
        0.5
      return(p_shift_one_cell_this_sigma * (n_episodes - 1))
    }))
    
    # Estimate the error of our chosen global scale hyperprior
    computed_num_expected_shifts <- sum(probs * num_expected_shifts)
    return(list(hyperprior = zeta, E.n = computed_num_expected_shifts))
  }

.titleFormatLabeller <- function(string) {
  lapply(string, .titleFormat)
}

# capitalize and remove hyphens
.titleFormat <- function(string) {
  string <- gsub("-", " ", string)
  string <- gsub("_", " ", string)
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  return(string)
}

### Functions required by densiTreeWithBranchData
# attribute colors to a vector based the value in a range
color_gradient <-
  function(x,
           intervals = seq(0, 11, 0.1),
           colors = c("red", "yellow", "green"),
           bias = 1) {
    colfun <- grDevices::colorRampPalette(colors, bias = bias)
    return(colfun(length(intervals)) [findInterval(x,
                                                   intervals,
                                                   all.inside = TRUE)])
  }

# function to sort a treedata
sort_tips <- function(x) {
  x <- reorder_treedata(x)
  nTip <- as.integer(length(x@phylo$tip.label))
  e2 <- x@phylo$edge[, 2]
  x@data <-
    x@data[c(e2[e2 <= nTip], (nTip + 1):(nTip + x@phylo$Nnode)), ]
  x@phylo$tip.label <- x@phylo$tip.label[e2[e2 <= nTip]]
  x@phylo$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
  x
}

# idem but with phylo
sort_tips_phylo <- function(x) {
  x <- ape::reorder.phylo(x)
  nTip <- as.integer(length(x$tip.label))
  e2 <- x$edge[, 2]
  x$tip.label <- x$tip.label[e2[e2 <= nTip]]
  x$edge[e2 <= nTip, 2] <- as.integer(1L:nTip)
  x
}

# get MRCA height from tree(s)
get_MRCA_heights <- function(x) {
  fun <- function(t)
    max(ape::node.depth.edgelength(t))
  height <- NULL
  if (inherits(x, "phylo"))
    height <- fun(x)
  if (inherits(x, "multiPhylo")) {
    if (!is.null(attr(x, "TipLabel"))) {
      x <- ape::.uncompressTipLabel(x)
      x <- unclass(x)
      height <- vapply(x, fun, 0)
    }
    else {
      x <- unclass(x)
      height <- vapply(x, fun, 0)
    }
  }
  else {
    height <- vapply(x, fun, 0)
  }
  height
}

# add tip labels to a tree plot - copied from phangorn
add_tiplabels <-
  function(xy,
           tip.label,
           direction,
           adj,
           font,
           srt = 0,
           cex = 1,
           col = 1,
           label_offset = 0) {
    direction <-
      match.arg(direction,
                c("rightwards", "leftwards",  "upwards",
                  "downwards"))
    horizontal <- direction %in% c("rightwards", "leftwards")
    nTips <- length(tip.label)
    xx <- rep(1, nrow(xy))
    yy <- xy[, 2]
    if (direction == "leftwards" |
        direction == "downwards")
      xx <- xx * 0
    if (!horizontal) {
      #    tmp <- yy
      yy <- xx
      xx <- xy[, 1]
    }
    MAXSTRING <- max(strwidth(tip.label, cex = cex))
    loy <- 0
    if (direction == "rightwards")
      lox <- label_offset + MAXSTRING * 1.05 * adj
    if (direction == "leftwards")
      lox <- -label_offset - MAXSTRING * 1.05 * (1 - adj)
    if (!horizontal) {
      psr <- par("usr")
      MAXSTRING <-
        MAXSTRING * 1.09 * (psr[4] - psr[3]) / (psr[2] - psr[1])
      loy <- label_offset + MAXSTRING * 1.05 * adj
      lox <- 0
      srt <- 90 + srt
      if (direction == "downwards") {
        loy <- -loy
        srt <- 180 + srt
      }
    }
    text(
      xx[1:nTips] + lox,
      yy[1:nTips] + loy,
      tip.label,
      adj = adj,
      font = font,
      srt = srt,
      cex = cex,
      col = col
    )
  }

# adapted from treeplyr
# https://github.com/uyedaj/treeplyr/blob/master/R/treeplyr_functions.R
# treeplyr::reorder (but not equivalent)
reorder_treedata <- function(tdObject, order = "postorder") {
  dat.attr <- attributes(tdObject@data)
  phy <- tdObject@phylo
  ntips <- length(phy$tip.label)
  phy$node.label <- (ntips + 1):(ntips + phy$Nnode)
  phy <- ape::reorder.phylo(phy, order)
  index <- match(tdObject@phylo$tip.label, phy$tip.label)
  index.node <- match((ntips + 1):(ntips + phy$Nnode), phy$node.label)
  
  tdObject@data <- tdObject@data[c(index, index.node), ]
  attributes(tdObject@data) <- dat.attr
  attributes(tdObject)$tip.label <- phy$tip.label
  tdObject@phylo <- phy
  
  tdObject
}

## End functions required by densiTreeWithBranchData

.removeNull <- function(x) {
  res <- x[which(!unlist(lapply(x, is.null)))]
}

# set prob factors
.set_pp_factor_range <-
  function(t, include_start_states, n_states = 1) {
    # what is the ancestral state name tag?
    if (include_start_states) {
      state_pos_str_base <- c("start_state_", "end_state_")
    } else {
      state_pos_str_base <- c("anc_state_")
    }
    
    # create list of ancestral state name tags
    state_pos_str_to_update <-
      c(unlist(lapply(1:n_states, function(x) {
        paste(state_pos_str_base, x, "_pp", sep = "")
      })))
    
    # overwrite state labels
    for (m in state_pos_str_to_update)
    {
      x_state <- attributes(t)$data[[m]]
      #levels(x_state) = c(levels(x_state))
      attributes(t)$data[[m]] <- x_state
    }
    return(t)
  }

.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)),
        substring(s, 2),
        sep = "",
        collapse = " ")
}

.matchNodesTreeData <- function(treedata, phy) {
  # get some useful info
  num_sampled_anc <- sum(phy$node.label != "")
  num_tips        <- length(phy$tip.label)
  num_nodes       <- phy$Nnode
  sampled_ancs    <- which(tabulate(phy$edge[, 1]) == 1)
  tip_indexes     <- 1:(num_tips + num_sampled_anc)
  node_indexes    <- (num_tips + num_sampled_anc) + num_nodes:1
  
  node_map     <-
    data.frame(R = 1:(num_tips + num_nodes),
               Rev = NA,
               visits = 0)
  current_node <- num_tips + 1
  k <- 1
  t <- 1
  
  while (TRUE) {
    # compute the number of descendants of this tip
    current_num_descendants <- sum(phy$edge[, 1] == current_node)
    
    if (current_node <= num_tips) {
      treedata_node <-
        which(as.character(treedata@data$node) == current_node)
      node_map$Rev[node_map$R == current_node] <-
        as.numeric(treedata@data[treedata_node, ]$index)
      current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
      t <- t + 1
      
    } else if (current_node %in% sampled_ancs) {
      if (node_map$visits[node_map$R == current_node] == 0) {
        node_map$Rev[node_map$R == current_node] <- node_indexes[k]
        k <- k + 1
      }
      node_map$visits[node_map$R == current_node] <-
        node_map$visits[node_map$R == current_node] + 1
      
      if (node_map$visits[node_map$R == current_node] == 1) {
        # go left
        current_node <- phy$edge[phy$edge[, 1] == current_node, 2][1]
      } else if (node_map$visits[node_map$R == current_node] == 2) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
        }
      }
      
    } else {
      if (node_map$visits[node_map$R == current_node] == 0) {
        node_map$Rev[node_map$R == current_node] <- node_indexes[k]
        k <- k + 1
      }
      node_map$visits[node_map$R == current_node] <-
        node_map$visits[node_map$R == current_node] + 1
      
      num_visits <- node_map$visits[node_map$R == current_node]
      
      if (num_visits <= current_num_descendants) {
        # go to next descendant
        current_node <-
          phy$edge[phy$edge[, 1] == current_node, 2][current_num_descendants -
                                                       num_visits + 1]
      } else if (num_visits > current_num_descendants) {
        # go down
        if (current_node == num_tips + 1) {
          break
        } else {
          current_node <- phy$edge[phy$edge[, 2] == current_node, 1]
        }
      }
      
    }
    
  }
  
  return(node_map[, 1:2])
  
}

funcShaded <- function(x, qLower, qUpper, mean, sd) {
  y = dnorm(x, mean = mean, sd = sd)
  lower_bound <- qnorm(qLower, mean = mean, sd = sd)
  upper_bound <- qnorm(qUpper, mean = mean, sd = sd)
  y[x < lower_bound] <- NA
  y[x > upper_bound] <- NA
  return(y)
}

get_legend2 <- function(plot, legend = NULL) {
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}


###################################################
#                                                 #
#                 Stateless pruning               #
#                                                 #
###################################################
postorder <- function(node_index, edge, tree, continuousChar,
                      μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta){
  ntip = length(tree$tip.label)
  #root_node = ntip + 1
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- postorder(left, edge, tree, continuousChar,
                             μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
    μ <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- postorder(right, edge, tree, continuousChar,
                              μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
    μ <- output_right[[1]]
    V <- output_right[[2]]
    log_norm_factor <- output_right[[3]]
    
    
    bl_left = branch_lengths[left_edge] # all branch of left child edge
    bl_right = branch_lengths[right_edge] # all branch of right child edge
    
    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
    
    v_left = sigma2/(2*alpha) *expm1(2.0*alpha*bl_left)
    var_left = v_left + V[left] * exp(2.0 * alpha * bl_left)
    
    v_right = sigma2/(2*alpha) *expm1(2.0*alpha*bl_right)
    var_right = v_right + V[right] * exp(2.0 * alpha * bl_right)
    
    # 2) mean of the normal variable
    mean_left = exp(alpha*bl_left)*(μ[left] - theta) + theta
    mean_right = exp(alpha*bl_right)*(μ[right] - theta) + theta
    
    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    μ[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node
    
    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    log_nf_left = bl_left * alpha
    log_nf_right = bl_right * alpha
    
    contrast = mean_left - mean_right
    a = -(contrast*contrast / (2*(var_left+var_right)))
    b = log(2*pi*(var_left+var_right))/2.0
    #b = log(2*pi)/2.0 + log(var_left+var_right)/2.0
    log_nf = log_nf_left + log_nf_right + a - b
    log_norm_factor[node_index] = log_nf
    
    
    return(list(μ, V, log_norm_factor))
  }
  
  
  # if is tip
  else{
    #edge_index = which(edge[,2] == node_index) # find edge index by tip node index
    #subedge = tree$maps[[edge_index]]
    species = tree$tip.label[node_index]
    
    μ[node_index] = as.numeric(continuousChar[[which(names(continuousChar) == species)]])
    V[node_index] = 0.0 ## if there is no observation error
    
    return(list(μ, V, log_norm_factor))
  }
}


## Pruning method (stateless)
logL_pruning <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes
  
  V = numeric(max_node_index)
  μ = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)
  
  branch_lengths = tree$edge.length
  
  root_index = ntip + 1
  
  
  
  
  output <- postorder(root_index, edge, tree, continuousChar,
                      μ, V, log_norm_factor, branch_lengths, alpha, sigma2, theta)
  μ <- output[[1]]
  V <- output[[2]]
  log_norm_factor <- output[[3]]
  
  ## assume root value equal to theta
  μ_root = μ[root_index]
  v_root = V[root_index]
  lnl = dnorm(theta, mean = μ_root, sd = sqrt(v_root), log = TRUE) # are \theta and \mu in correct positions?
  
  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}

###################################################
#                                                 #
#            State-dependent pruning              #
#                                                 #
###################################################
sd_postorder <- function(node_index, edge, tree, continuousChar,
                         μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta){
  ntip = length(tree$tip.label)
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    output_left <- sd_postorder(left, edge, tree, continuousChar,
                                μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta)
    μ <- output_left[[1]]
    V <- output_left[[2]]
    log_norm_factor <- output_left[[3]]
    
    output_right <- sd_postorder(right, edge, tree, continuousChar,
                                 μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta)
    μ <- output_right[[1]]
    V <- output_right[[2]]
    log_norm_factor <- output_right[[3]]
    
    
    sub_bl_left = subedges_lengths[[left_edge]] # all subedges of left child edge
    sub_bl_right = subedges_lengths[[right_edge]] # all subedges of right child edge
    
    # for the sake of readability, computation of variance, mean, and log_nf are done in separate loops
    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
    ## Is 'delta_left* exp(2.0 * alpha * bl_left)' added in each sub-edge?
    
    delta_left = V[left]
    v_left = 0 # initialise v_left
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left[i])
      v_left = sigma2[[state]]/(2*alpha[[state]]) *expm1(2.0*alpha[[state]]
                                                         *sub_bl_left[[i]])
      delta_left = v_left + delta_left * exp(2.0 * alpha[[state]] * sub_bl_left[[i]])
    }
    
    delta_right = V[right]
    v_right = 0 # initialise v_right
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right[i])
      v_right = sigma2[[state]]/(2*alpha[[state]]) *expm1(2.0*alpha[[state]]*sub_bl_right[[i]])
      delta_right = v_right + delta_right * exp(2.0 * alpha[[state]] * sub_bl_right[[i]])
    }
    
    var_left = delta_left
    var_right = delta_right
    
    # 2) mean of the normal variable
    mean_left = μ[left]
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left[i])
      mean_left = exp(alpha[[state]]*sub_bl_left[[i]])*(mean_left - theta[[state]]) + theta[[state]]
    }
    
    mean_right = μ[right]
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right[i])
      mean_right = exp(alpha[[state]]*sub_bl_right[[i]])*(mean_right - theta[[state]]) + theta[[state]]
    }
    
    
    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    μ[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node
    
    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    ## this is the problem. I think in RevBayes we compute log_nf with the oldest sub-edge only
    log_nf_left = 0
    for (i in rev(1:length(sub_bl_left))){
      state <- names(sub_bl_left[i])
      log_nf_left = log_nf_left + sub_bl_left[[i]] * alpha[[state]]
    }
    
    log_nf_right = 0
    for (i in rev(1:length(sub_bl_right))){
      state <- names(sub_bl_right[i])
      log_nf_right = log_nf_right + sub_bl_right[[i]] * alpha[[state]]
    }
    
    contrast = mean_left - mean_right
    a = -(contrast*contrast / (2*(var_left+var_right)))
    b = log(2*pi*(var_left+var_right))/2.0
    #b = log(2*pi)/2.0 + log(var_left+var_right)/2.0
    log_nf = log_nf_left + log_nf_right + a - b
    log_norm_factor[node_index] = log_nf
    
    return(list(μ, V, log_norm_factor))
  }
  
  
  # if is tip
  else{
    species = tree$tip.label[node_index]
    
    μ[node_index] = as.numeric(continuousChar[[which(names(continuousChar) == species)]])
    V[node_index] = 0.0 ## if there is no observation error
    
    return(list(μ, V, log_norm_factor))
  }
}

sd_logL_pruning <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes
  
  V = numeric(max_node_index)
  μ = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)
  
  subedges_lengths = tree$maps
  
  root_index = ntip + 1
  
  output <- sd_postorder(root_index, edge, tree, continuousChar,
                         μ, V, log_norm_factor, subedges_lengths, alpha, sigma2, theta)
  μ <- output[[1]]
  V <- output[[2]]
  log_norm_factor <- output[[3]]
  
  ## assume root value equal to theta
  μ_root = μ[root_index]
  v_root = V[root_index]
  left_edge_from_root <- which(edge[,1] == root_index)[1] # obtain left child edge index of root node
  left_subedges_from_root <- subedges_lengths[[left_edge_from_root]] # obtain sub-edge lengths
  ### note here
  root_state = names(head(left_subedges_from_root, 1)) # obtain root state, assuming it equals last state at left child edge
  lnl = dnorm(theta[[root_state]], mean = μ_root, sd = sqrt(v_root), log = TRUE)
  
  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}


###################################################
#                                                 #
#                  Stateless vcv                  #
#                                                 #
logL_vcv <- function(tree, continuousChar, alpha, sigma2, theta){
  ntip <- length(tree$tip.label)
  mrca1 <- ape::mrca(tree) # get the ancestral node label for each pair of tips
  times <- ape::node.depth.edgelength(tree) # get time at each node from root
  ta <- matrix(times[mrca1], nrow=ntip, dimnames = list(tree$tip.label, tree$tip.label)) # get time of divergence for each pair of tips
  T.term <- times[1:ntip] # get time at tips
  tia <- times[1:ntip] - ta
  tja <- t(tia)
  tij <- tja + tia # distance in time unit between two tips
  
  vy = sigma2 / (2*alpha)
  
  #V = vy * (1 - exp(-2 * alpha * ta)) * exp(-alpha * tij)
  V = vy * -1 * expm1(-2 * alpha * ta) * exp(-alpha * tij) ### ta = time tgt; tij = time not tgt (sum of two branches)
  
  X = matrix(1, ntip)
  
  C = chol(V) # upper triangular matrix
  L = t(C) # lower triangular matrix
  log_det_V = 0
  for (i in 1:ntip){
    log_det_V = log_det_V + log(L[i,i])
  }
  log_det_V = log_det_V *2.0 # equals to julia implementation to 12 sig. fig.
  
  y = NULL
  for (species in tree$tip.label){
    y[species] = as.numeric(continuousChar[species])
  }
  
  r = solve(L) %*% y - solve(L) %*% X * theta # what does de-correlated residuals mean?
  
  # res = - (n/2) * log(2*pi) - 0.5 * log_det_V - 0.5 * dot(r, r)
  #     = exp(-n/2)^(2*pi) * exp(-0.5)^det_V * exp(-0.5)^dot(r, r) ?
  res = 0.0
  res = res - 0.5 * ntip * log(2*pi)
  res = res - 0.5 * log_det_V
  res = res - 0.5 * dot(r, r) # what is r and what is  dot product of r?
  
  return(res)
}

###################################################
#                                                 #
#               State-dependent vcv               #
#                                                 #
###################################################

# find parent node by providing child node
parentNode <- function(tree, x){
  m <- which(tree$edge[, 2] == x)
  return(tree$edge[m, 1])
}

# find nodes along a lineage towards root node by providing initial child (presumably tip) node
nodesAlongLineage <- function(tree, old_node, young_node){
  k <- young_node
  while(young_node != old_node){
    k <- c(k, parentNode(tree, young_node))
    young_node <- tail(k, n = 1)
  }
  return(k)
}

# find subedges of a lineage
lineage.constructor <- function(tree, root_node, e){
  nodes <- nodesAlongLineage(tree, root_node, e)
  edges <- which(tree$edge[,2] %in% nodes) # from root to tip
  subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # tip to root
  
  state_changes <- names(subedge_lengths) # from tip to root
  #state_changes <- c(state_changes[1], state_changes) # add root state, assuming root state equals the state of the closest subedge
  
  #lineage$state_indicator <- lapply(all_states, function(x) {res <- match(lineage$state_changes, x); res[is.na(res)] <- 0; return(res)})
  #names(lineage$state_indicator) <- all_states
  
  # recording time-related numbers of each subedge (root is a subedge with length = 0)
  #times <- cumsum(unname(subedge_lengths))
  #time_tip <- tail(times, n = 1)
  
  #time_begin <- time_tip - c(0, head(times, n = -1))
  #time_end <- time_tip - times
  #time_span <- time_begin - time_end
  
  return(tibble(state = state_changes,
                #time_begin = time_begin,
                #time_end = time_end,
                time_span = subedge_lengths))
}

# not yet finished - weight matrix function
## need updates
weights.lineage <- function(tree, alpha, e){
  root_node = length(tree$tip.label) + 1
  lineage <- lineage.constructor(tree, root_node, e)
  lineage[["alpha"]] = alpha[lineage[["state"]]]
  
  W = matrix(0, ncol = length(alpha), nrow = 1)
  colnames(W) = sort(names(alpha))
  
  if (length(lineage[[1]]) > 1){
    lineage <- lineage %>%
      mutate(
        exp1 = -1 * expm1(-1 * alpha * time_span),
        sum2_temp = -1 * alpha * time_span)
    lineage$exp1[length(lineage$exp1)] = 1
    lineage$sum2 = 0
    
    for (i in 2:length(lineage[[1]])){
      lineage$sum2[i] = lineage$sum2_temp[i-1]
      lineage$sum2_temp[i] = lineage$sum2[i] + lineage$sum2_temp[i]
    }
    
    all_weights = lineage %>% mutate(exp_final = exp1 * exp(sum2)) %>% 
      group_by(state) %>% 
      summarise(weight = sum(exp_final))
    
    for (i in 1:nrow(all_weights)){
      W[, all_weights$state[i]] = all_weights$weight[i]
    }
  } else {
    W[, lineage$state[1]] = 1
  }
  
  return(W)
}

# combine to form weight matrix
weight.matrix <- function(tree, alpha){
  ntip = length(tree$tip.label)
  weight_matrix = matrix(0, nrow = ntip, ncol = length(alpha))
  rownames(weight_matrix) <- tree$tip.label
  colnames(weight_matrix) <- c(sort(names(alpha)))
  for (i in 1:ntip){
    weight_matrix[i, ] <- weights.lineage(tree, alpha, i)
  }
  return(weight_matrix)
}


cov.accum <- function(tree, mrca_node, alpha, sigma2){
  root_node = length(tree$tip.label) + 1
  if (mrca_node == root_node){
    cov_accum = 0.0
  } else {
    nodes <- nodesAlongLineage(tree, root_node, mrca_node)
    edges <- which(tree$edge[,2] %in% nodes) # from root to mcra_node
    subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from mcra_node to root
    
    subedge_lengths <- tibble(state = names(subedge_lengths),
                              time_span = subedge_lengths,
                              alpha = alpha[names(subedge_lengths)],
                              sigma2 = sigma2[names(subedge_lengths)]) %>% 
      mutate(exp1 = -1 * expm1(-2 * alpha * time_span),
             sum2_temp = -2 * alpha * time_span)
    subedge_lengths$sum2= 0
    
    if (length(subedge_lengths[[1]]) == 1){
      subedge_lengths = subedge_lengths %>% 
        mutate(cov = sigma2 / (2 * alpha) * exp1)
      cov_accum = subedge_lengths$cov[[1]]
    } else {
      for (i in 2:length(subedge_lengths[[1]])){
        subedge_lengths$sum2[i] = subedge_lengths$sum2_temp[i-1]
        subedge_lengths$sum2_temp[i] = subedge_lengths$sum2[i] + subedge_lengths$sum2_temp[i]
      }
      cov_accum = subedge_lengths %>% mutate(exp3 = exp1 * exp(sum2)) %>% 
        group_by(state) %>% 
        summarise(sum4 = sum(sigma2 / (2 * alpha) * exp3)) %>% 
        reframe(sum_final = sum(sum4)) %>% 
        unlist() %>% 
        unname()
    }
  }
  return(cov_accum)
}

cov.loss <- function(tree, mrca_node, alpha, tip){
  if (mrca_node == tip){
    cov_loss_rate = 0
  } else {
    nodes <- nodesAlongLineage(tree, mrca_node, tip)
    nodes <- head(nodes, n = -1)
    edges <- which(tree$edge[,2] %in% nodes) # from root to mcra_node
    subedge_lengths <- rev(unlist(lapply(edges, function(i) tree$maps[[i]]))) # from mcra_node to root
    subedge_lengths <- tibble(time_span = subedge_lengths,
                              alpha = alpha[names(subedge_lengths)])
    cov_loss_rate = subedge_lengths %>% 
      mutate(sum1 = -1 * alpha * time_span) %>% 
      reframe(sum_final = sum(sum1))
  }
  return(cov_loss_rate)
}

vcv.pairwise <- function(tree, alpha, sigma2, tip1, tip2){
  mrca_node <- ape::mrca(tree)[tip1, tip2]
  cov_accum = cov.accum(tree, mrca_node, alpha, sigma2)
  cov_loss1 = cov.loss(tree, mrca_node, alpha, tip1)
  cov_loss2 = cov.loss(tree, mrca_node, alpha, tip2)
  cov = cov_accum * exp(cov_loss1 + cov_loss2)
  return(unlist(unname(cov)))
}


vcv.matrix <- function(tree, alpha, sigma2){
  ntip <- length(tree$tip.label)
  V <- matrix(nrow = ntip, ncol = ntip)
  j = ntip
  while (j != 0){
    for (i in 1:ntip){
      V[i,j] <- vcv.pairwise(tree, alpha, sigma2, i, j)
      V[j,i] <- V[i,j]
    }
    j = j-1
  }
  colnames(V) <- tree$tip.label
  rownames(V) <- tree$tip.label
  return(V)
}


sd_logL_vcv <- function(tree, continuousChar, alpha, sigma2, theta){
  alpha = alpha[sort(names(alpha))]
  sigma2 = sigma2[sort(names(sigma2))]
  theta = theta[sort(names(theta))]
  theta = as.matrix(theta, nrow = 3)
  
  ntip <- length(tree$tip.label)
  V = vcv.matrix(tree, alpha, sigma2)
  
  W = weight.matrix(tree, alpha)
  
  C = chol(V) # upper triangular matrix
  L = t(C) # lower triangular matrix
  log_det_V = 0
  for (i in 1:ntip){
    log_det_V = log_det_V + log(L[i,i])
  }
  log_det_V = log_det_V * 2.0 # equals to julia implementation to 12 sig. fig.
  
  y = NULL
  for (species in tree$tip.label){
    y[species] = as.numeric(continuousChar[species])
  }
  
  # inverse of L
  r = solve(L) %*% y - solve(L) %*% W %*% theta # what does de-correlated residuals mean?
  
  # res = - (n/2) * log(2*pi) - 0.5 * log_det_V - 0.5 * dot(r, r)
  #     = exp(-n/2)^(2*pi) * exp(-0.5)^det_V * exp(-0.5)^dot(r, r) ?
  res = 0.0
  res = res - (ntip/2) * log(2*pi)
  res = res - 0.5 * log_det_V
  res = res - 0.5 * dot(r, r) # is it dot product? what is dot product of r?
  
  return(res)
}
