HierName <- function(Smat) {
  h.depth <- max(colSums(Smat))
  names.list <- list(length = h.depth)
  all_names = rownames(Smat)
  names.list[[1L]] <- all_names[which.max(rowSums(Smat))]
  remaining.nodes = setdiff(all_names, names.list[[1L]])
  if (h.depth > 1L) {
    for (i in 2L:h.depth) {
      if(is.null(remaining.nodes)) break
      next_hier_info = findNextHier(Smat[remaining.nodes,])
      names.list[[i]] <- next_hier_info$nodes
      remaining.nodes <- next_hier_info$remaining.nodes
    }
    names(names.list) <- paste("Level", 0:(h.depth-1))
  }
  return(names.list)
}

findNextHier = function(Smat)
{
  #arrange with largest number of summands first (i.e. largest number of bottom TS)
  Smat[order(-rowSums(Smat)),]
  all.names = rownames(Smat)

  if(max(rowSums(Smat))>1) # check if only leave nodes left
    level.list = c(all.names[1]) else return(list(remaining.nodes = NULL, nodes = all.names))

  n_check = nrow(Smat)
  for(k in 2:n_check){
    c_match = 0
    #check whether the new row lies in the complement of the rows columns
    for(c_label in level.list)
    {
      if(any(Smat[c_label,]- Smat[k,]==-1)) c_match = c_match+1
    }
    #add the row if it lies in the complement for all discovered rows
    if(c_match == length(level.list)) level.list = c(level.list, all.names[k])
  }

  return(list(remaining.nodes = setdiff(all.names, level.list), nodes = level.list))
}

findNodes = function(Smat, h.levels)
{
  my_nodes = list()
  my_counts = list()


  for(k in 1:(length(h.levels)-1))
  {
    my_nodes[[k]] = list()
    my_counts[[k]] = list()
    av.nodes = h.levels[[k+1]]
    for(root_node in h.levels[[k]])
    {
      my_nodes[[k]][[root_node]] = c()
      #check if any next node is subsumed under the root node
      if(length(av.nodes)==0){
        my_counts[[k]][[root_node]] = 0
        next
      }
      for(c_node in av.nodes)
      {
        if(!any(Smat[root_node,]-Smat[c_node,]==-1))
        {
          my_nodes[[k]][[root_node]] = c(my_nodes[[k]][[root_node]], c_node)
        }
      }
      my_counts[[k]][[root_node]] = length(my_nodes[[k]][[root_node]])
      av.nodes = setdiff(av.nodes, my_nodes[[k]][[root_node]])

    }
    my_counts[[k]] = unlist(my_counts[[k]])
  }

  names(my_counts) = paste0("Level ", 0:(length(my_counts)-1))
  names(my_nodes) = paste0("Level ", 0:(length(my_nodes)-1))
  return(list(counts = my_counts, hierNodes = my_nodes))
}
