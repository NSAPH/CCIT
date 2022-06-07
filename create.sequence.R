create.sequence = function(exploration.dat, ulist.used, parms.used, f){ # add more args
  # This function calculates the sequence of candidate trees.
  
  # Fit a large tree using the user written splitting functions
  a <- rpart(
    #id ~ em1 + em2 + em3 + em4, 
    f,
    data = exploration.dat, 
    method   = ulist.used,
    parms    = parms.used,
    control  = rpart.control(cp = -10, minbucket = 1, maxsurrogate = 0, maxcompete = 0))
  
  if (dim(a$frame)[1] == 1){ # Deal with the root only tree
    
    tree.list = list(a)
    g.h.list = list(Inf)
    return(list(tree.list = tree.list, g.h.list = g.h.list))
    
  } else {
    
    # Finding which variables are leaf nodes
    is.leaf <- (a$frame$var == "<leaf>")
    
    # A function adapted from the partykit package that identifies the rows of the frame 
    # which correspond to child nodes of row i in frame matrix
    rpart.kids <- function(i, is.leaf) {
      if (is.leaf[i]) return(NULL)
      else return(c(i + 1L, 
                    which((cumsum(!is.leaf[-(1L:i)]) + 1L) == cumsum(is.leaf[-(1L:i)]))[1L] + 1L + i))
    }
    
    # Finding goodness of the split
    a$frame$split.stat = 0
    a$frame$split.stat[!is.leaf] = a$splits[, 3]
    
    # Calculating the g(h) parameter for each non-terminal node
    g.h = rep(0, nrow(a$frame))
    for(i in 1:nrow(a$frame)){
      if(is.leaf[i]){g.h[i] = Inf} else{
        # Find all kids of node i
        kids.i = i
        stop.loop = FALSE
        while(stop.loop == FALSE){
          kids.old = kids.i
          for(j in 1:length(kids.i)){
            kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf)))
          }
          if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
          # Calculating g.h for node i
          g.h[i] = sum(a$frame$split.stat[kids.i])/sum(a$frame$split.stat[kids.i] != 0)
        }
      }
    }
    
    # Adding g.h to frame
    a$frame$g.h = g.h
    
    # Start pruning
    # First tree is the large tree
    tree.list = list(a)
    g.h.list = list(0)
    # pruned.node = list(NULL)
    stop.prune = FALSE
    k = 1
    
    while(stop.prune == FALSE){
      tree.used = tree.list[[k]]
      # Calculating the g(h) parameter for each non-terminal node
      tree.used$frame$g.h = rep(0, nrow(tree.used$frame))
      is.leaf.prune <- (tree.used$frame$var == "<leaf>")
      # Setting splitting statistics for new terminal nodes to 0
      tree.used$frame$split.stat[is.leaf.prune] = 0
      
      # Calculating the g(h) function for each non-terminal node
      for(i in 1:nrow(tree.used$frame)){
        if(is.leaf.prune[i]){tree.used$frame$g.h[i] = Inf} else{
          # Find all kids of node i
          kids.i = i
          stop.loop = FALSE
          while(stop.loop == FALSE){
            kids.old = kids.i
            for(j in 1:length(kids.i)){
              kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
            }
            if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
            tree.used$frame$g.h[i] = sum(tree.used$frame$split.stat[kids.i])/sum(tree.used$frame$split.stat[kids.i] != 0)
          }
        }
      }
      
      # Finding the value which minimizes g(h) (among internal nodes)
      to.prune = which.min(tree.used$frame$g.h)
      #print(tree.used$frame[to.prune,])
      # Finding the minimum g.h value
      g.h.min = min(tree.used$frame$g.h)
      
      # Find all kids of node to.prune
      kids.i = to.prune
      stop.loop = FALSE
      while(stop.loop == FALSE){
        kids.old = kids.i
        for(j in 1:length(kids.i)){
          kids.i = unique(c(kids.i, rpart.kids(kids.i[j], is.leaf.prune)))
        }
        if(length(kids.i) == length(kids.old)){stop.loop = TRUE}          
      }
      
      
      # Finding number of splits to prune
      split.to.prune = length(kids.i[which(!is.leaf.prune[kids.i])])
      
      # Creating the new splits and frames for new tree
      splits.new = tree.used$splits[-c(sum(!is.leaf.prune[1:to.prune]):(sum(!is.leaf.prune[1:to.prune]) + split.to.prune - 1)), ]
      frame.new = tree.used$frame[-setdiff(kids.i, to.prune), ]
      
      # Changing all nodes that were internal nodes and are now terminal node to terminal nodes
      frame.new$var[to.prune] =  "<leaf>"
      
      tree.new = tree.used
      tree.new$frame = frame.new
      if(class(splits.new) == "matrix"){tree.new$splits = splits.new}
      if(class(splits.new) == "numeric"){
        tree.new$splits = matrix(splits.new, nrow = 1)
        colnames(tree.new$splits) = colnames(tree.used$splits)
      }
      
      # Changing the terminal node for $where in rpart object
      tree.new$where = tree.used$where
      tree.new$where[tree.new$where %in% kids.i] = to.prune
      tree.new$where[tree.new$where > max(kids.i)] = tree.new$where[tree.new$where > max(kids.i)] - length(kids.i) + 1
      tree.new$where = as.integer(tree.new$where)
      
      k = k+1
      # Add tree and lambda to the list
      tree.list[[k]] <- tree.new
      g.h.list[[k]] <- g.h.min
      # pruned.node[[k]] <- 
      
      if(sum(tree.new$frame$var == "<leaf>") == 1){stop.prune = TRUE}  
    }
    return(list(tree.list = tree.list, g.h.list = g.h.list))
  }
}
