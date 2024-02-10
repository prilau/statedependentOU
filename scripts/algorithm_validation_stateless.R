library(ape)
library(phytools)
library(knitr)

## Postorder function
postorder <- function(node_index, edge, tree, continuousChar,
                      μ, V, log_norm_factor, branch_lengths, σ2, α, θ){
  ntip = length(tree$tip.label)
  #root_node = ntip + 1
  
  # if is internal node
  if (node_index > ntip){
    
    left_edge  = which(edge[,1] == node_index)[1] # index of left child edge
    right_edge = which(edge[,1] == node_index)[2] # index of right child edge
    left = edge[left_edge,2] # index of left child node
    right = edge[right_edge,2] # index of right child node
    
    postorder(left, edge, tree, continuousChar,
              μ, V, log_norm_factor, branch_lengths, σ2, α, θ)
    postorder(right, edge, tree, continuousChar,
              μ, V, log_norm_factor, branch_lengths, σ2, α, θ)
    
    bl_left = branch_lengths[left_edge] # all branch of left child edge
    bl_right = branch_lengths[right_edge] # all branch of right child edge
    
    # 1) variance of the normal variable: this branch (v_left) and the subtree (V[left])
    
    v_left = σ2/(2*α) *expm1(2.0*α*bl_left)
    var_left = v_left + V[left] * exp(2.0 * α * bl_left)
    
    v_right = σ2/(2*α) *expm1(2.0*α*bl_right)
    var_right = v_right + V[right] * exp(2.0 * α * bl_right)
    
    # 2) mean of the normal variable
    mean_left = exp(α*bl_left)*(μ[left] - θ) + θ
    mean_right = exp(α*bl_right)*(μ[right] - θ) + θ
    
    ## compute the mean and variance of the node
    mean_ancestor = (mean_left * var_right + mean_right * var_left) / (var_left + var_right)
    μ[node_index] = mean_ancestor
    var_node = (var_left * var_right) / (var_left + var_right)
    V[node_index] = var_node
    
    ## compute the normalizing factor, the left-hand side of the pdf of the normal variable
    ## this is the problem. I think in RevBayes we compute log_nf with the oldest sub-edge only
    log_nf_left = bl_left * α
    log_nf_right = bl_right * α

    contrast = mean_left - mean_right
    a = -(contrast*contrast / (2*(var_left+var_right)))
    b = log(2*pi*(var_left+var_right))/2.0
    #b = log(2*pi)/2.0 + log(var_left+var_right)/2.0
    log_nf = log_nf_left + log_nf_right + a - b
    log_norm_factor[node_index] = log_nf
  }
  
  
  # if is tip
  else{
    #edge_index = which(edge[,2] == node_index) # find edge index by tip node index
    #subedge = tree$maps[[edge_index]]
    species = tree$tip.label[node_index]
    
    μ[node_index] = as.numeric(continuousChar[which(names(continuousChar) == species)][[1]])
    V[node_index] = 0.0 ## if there is no observation error
  }
}

## Pruning method (stateless)
logL_pruning <- function(tree, continuousChar, σ2, α, θ){
  ntip = length(tree$tip.label) # number of tips
  edge = tree$edge # equals tree[:edge] in Julia
  n_edges = length(edge[,1]) # number of edges
  max_node_index = max(tree$edge) # total number of nodes
  
  V = numeric(max_node_index)
  μ = numeric(max_node_index)
  log_norm_factor = numeric(max_node_index)
  
  branch_lengths = tree$edge.length
  
  root_index = ntip + 1
  
  postorder(root_index, edge, tree, continuousChar,
            μ, V, log_norm_factor, branch_lengths, σ2, α, θ)
  
  ## assume root value equal to theta
  μ = μ[root_index]
  v = V[root_index]
  lnl = dnorm(θ, mean = μ, sd = sqrt(v), log = TRUE) # are \theta and \mu in correct positions?
  
  ## add norm factor
  for (log_nf in log_norm_factor){
    lnl = lnl + log_nf
  }
  return(lnl)
}




# reading in tree in simmap format
tree <- read.simmap("data/1_validation/dummy_threetaxon_simmap.tre",format="phylip",version=1)

#reading in continuous character (somehow it recognise 0.1234 as 6 characters so I read in manually)
continuousChar <- read.nexus.data("data/1_validation/dummy_threetaxon_Continuous.nex")

logL_pruning(tree, continuousChar,
             σ2 = 2,
             α = 1.5,
             θ = 4)


