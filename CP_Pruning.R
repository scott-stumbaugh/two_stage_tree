# ---------- PRUNING + COMPLEXITY PARAMETER FUNCTIONS ------------
# ****************************************************************

# find next immediate set of nodes according to node numbering scheme
next_nodenum.fun <- function(a){c(2*a, 2*a + 1)}


# specify node and find all subsequent branch nodes
find_child_nodes.fun <- function(b, max_node){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # returns a vector of all node numbers corresponding to child nodes of a given parent
  # node (b) and the max node number in a tree (max_node)
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # find max power for base 2 binary recursive splitting
  n_max <- ceiling((log(max_node+1)-log(b+1))/log(2))
  
  # loop through all levels of binary splitting
  nodenum.vec <- c()
  for(n in 1:n_max){
    nodenum.vec <- c(nodenum.vec, 2^(n)*b + 0:(2^(n)-1))
  }
  
  return(nodenum.vec[nodenum.vec<=max_node])
}


# find a subtree by pruning
subtree.fun <- function(tree_model, nodenum){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # returns a subtree based on pruning a larger tree at a specified node
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # find branches to remove
  max_node <- max(tree_model$node_num)
  prune_branch.indx <- find_child_nodes.fun(nodenum, max_node)
  
  # pruned tree
  subtree <- tree_model[!(tree_model$node_num %in% prune_branch.indx), ]
  
  # make node where pruning occurred a leaf
  subtree[subtree$node_num==nodenum, 'leaf_flag'] <- 1
  
  return(subtree)
}


# find size of tree defined by number of terminal leaves
subtree_size.fun <- function(tree_model, node_prune){sum(subtree.fun(tree_model, node_prune)$leaf_flag)}


# find error sum of squared for a given tree model
subtree_sse.fun <- function(tree_model, node_prune){
  sum(subtree.fun(tree_model, node_prune)$sse*subtree.fun(tree_model, node_prune)$leaf_flag)
}


# find best subtree for every size |T|
opt_subtrees.fun <- function(tree_model){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # returns optimal subtree (min SSE criteria) for every subtree size <= full tree size
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # find all non-terminal (parent) nodes
  branch_node.vec <- tree_model[tree_model$leaf_flag==0, 'node_num']
  
  # find all subtree sizes associated with with pruning at each non-terminal node
  t.vec <- sapply(branch_node.vec, subtree_size.fun, tree_model=mytree.ff)
  sse.vec <- sapply(branch_node.vec, subtree_sse.fun, tree_model=mytree.ff)
  
  # combine node trimming location, size, and sse for all subtrees
  subtrees.df <- data.frame('pruned_node'=branch_node.vec, 'subtree_size'=t.vec, 'sse'=sse.vec)
  
  # find min SSE by subtree size
  best_subtrees.df <- aggregate(.~subtree_size, subtrees.df, min)
  
  # ...and don't forget to include the full exhaustive tree
  ft_max_node <- max(mytree.ff$node_num)
  ft_t <- sum(mytree.ff$leaf_flag)
  ft_sse <- sum(mytree.ff$leaf_flag*mytree.ff$sse)
  best_subtrees.df <- rbind(best_subtrees.df, data.frame('subtree_size'=ft_t, 'pruned_node'=ft_max_node, 'sse'=ft_sse))
  
  return(best_subtrees.df)
}


cp_tbl.fun <- function(tree_model){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # returns table of pruned subtree references and corresponding complexity parameters
  # for nested sequence of subtrees minimizing the cost function
  # 
  # define cost function C(T) = |T| + alpha*SSE(T)
  # where |T| is size of tree T (num of terminal leaves), alpha is cp, and SSE(T) is error sum of square for tree model T
  # solve for minimizing cost function across all subtrees using alpha = (|T_0| - |T|)/(SSE(T) - SSE(T_0)) for subtree T_0 and tree T
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # define subset of best candidate subtrees to minimize cost over
  cp.df <- opt_subtrees.fun(tree_model)
  
  # start with largest tree denoted by t
  max_size <- max(cp.df$subtree_size)
  indx_pointer <- which(cp.df$subtree_size==max_size)
  
  # iteratively search over small trees given larger tree
  cp_tbl.df <- data.frame(); cp_opt <- NA
  repeat{
    # identify largest tree in current iteration
    t <- cbind(cp.df[indx_pointer, ], 'alpha'=cp_opt)
    t_size <- t$subtree_size
    t_sse <- t$sse
    
    # results from previous iteration added to collective results
    cp_tbl.df <- rbind(cp_tbl.df, t)
    
    # define smaller subtree candidates (denoted t0) to search for minimizing cost function
    subtrees_temp <- cp.df[cp.df$subtree_size<t_size, ]
    
    # break if no more subtree candidates to search over
    if(nrow(subtrees_temp)==0){break}
    
    # compute cp over all candidate subtrees
    cp_temp <- c()
    for(i in 1:nrow(subtrees_temp)){
      t0_size <- subtrees_temp$subtree_size[i]
      t0_sse <- subtrees_temp$sse[i]
      
      alpha0 <- (t0_size - t_size)/(t_sse - t0_sse)
      
      cp_temp[i] <- alpha0
    }
    cp_opt <- max(cp_temp)
    indx_pointer <- ifelse(cp_opt>0, which.max(cp_temp), NA)
    
    # break if no more valid cost minimizing subtree candidates
    if(is.na(indx_pointer)){break}
  }
  
  # clean up
  cp_tbl.df <- cp_tbl.df[order(cp_tbl.df$subtree_size), c('pruned_node', 'subtree_size', 'sse', 'alpha')]
  
  return(cp_tbl.df)
}