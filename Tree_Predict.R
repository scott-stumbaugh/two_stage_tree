# ---------- PREDICTION ------------
# **********************************

# single value prediction
tree_pred_val.fun <- function(tree_model, x, nodenum){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # returns a single prediction given a single observation and tree model
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # traverse tree and predict when arriving at leaf
  repeat{
    # retrieve values from specified node
    node <- tree_model[tree_model$node_num==nodenum, ]
    
    # if we have arrived at a leaf then predict
    if(node$leaf_flag==1){
      pred <- node$mean
      return(pred)
    }
    
    # retrieve child node
    leftson <- tree_model[tree_model$node_num==next_nodenum.fun(nodenum)[1], ]
    
    # use child nodes to evaluate next direction in the tree
    splitvar <- leftson$split_var
    splitval <- leftson$split_val_num
    splitfac <- leftson$split_val_f
    
    if(!is.na(splitval)){ # x is continuous or discrete numeric
      nodenum <- ifelse(x[, splitvar]<=splitval, 2*nodenum, 2*nodenum + 1)
    }else{  # x is categorical factor
      nodenum <- ifelse(x[, splitvar]==splitfac, 2*nodenum, 2*nodenum + 1)
    }
  }
}


# predict entire data set
tree_pred.fun <- function(tree_model, x.dat){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # returns prediction on entire data set given a tree model
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # format data
  x.dat <- data.frame(x.dat)
  n <- nrow(x.dat)
  
  # loop through each row of data and compute prediction
  pred.vec <- c()
  for(i in 1:n){
    x_i <- x.dat[i, ]
    pred.vec[i] <- tree_pred_val.fun(tree_model=tree_model, x=x_i, 1)
  }
  
  return(pred.vec)
}