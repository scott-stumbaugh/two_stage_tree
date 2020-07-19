# required package(s)
library(rpart)




# ---------- TREE BUILDING FUNCTIONS  ------------
# ************************************************

# ----- total sum of squares for a vector of values
sstotal.fun <- function(y){
  ybar <- mean(y)
  return(sum((y-ybar)^2))
}


# ----- error sum of squares from rpart tree model
sse.fun <- function(tree.model){
  return(sum(residuals(tree.model)^2))
}


# ----- results from all splits 
splits.fun <- function(d, nvals, sstotal, minsplit, cp){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # iterate through every covariate and every value to create all initial splits
  # resulting in a left and right data then call rpart on left and right data to
  # get best second stage tree...this is essentially the main workhorse for the
  # two stage splitting functionality
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  splits.df <- data.frame()
  for(i in 1:(length(d)-1)){  # loops through covariates of data frame (excluding response, the last column)
    
    d.sort <- d[order(d[, i]), ]    # sort data by values in ith covariate
    x.name <- names(d.sort)[i]
    x.temp <- d.sort[, i]
    y.temp <- d.sort[, length(d)]
    
    # inner loops splits the data based on jth val in ith covariate using <= and > logic in the case of numeric data type
    if(!is.factor(x.temp)){   # x is continuous or discrete numeric
      x.unique <- unique(x.temp)
      
      for(j in 1:max((length(x.unique)-1), 1)){
        x.left <- d.sort[x.temp<=x.unique[j], -length(d)]   # exclude response
        y.left <- y.temp[x.temp<=x.unique[j]]
        
        # dat.left <- d.sort[1:j, ]
        # y.left <- y.temp[1:j]
        
        x.right <- d.sort[x.temp>x.unique[j], -length(d)]   # exclude response
        y.right <- y.temp[x.temp>x.unique[j]]
        
        # dat.right <- d.sort[(j+1):nvals, ]
        # y.right <- y.temp[(j+1):nvals]
        
        # second stage tree building 
        tree.left <- rpart(y.left ~ ., data=x.left, minsplit = minsplit, maxdepth = 1, cp=cp)
        xvar.left <- ifelse(tree.left$frame[1, 'var']=='<leaf>', NA, as.character(tree.left$frame[1, 'var']))
        xval.left <- ifelse(is.null(tree.left$splits[1, 'index']), NA, tree.left$splits[1, 'index'])
        n.left <- ifelse(tree.left$frame[1, 'var']=='<leaf>', NA, tree.left$frame[1, 'n'])
        
        tree.right <- rpart(y.right ~ ., data=x.right, minsplit = minsplit, maxdepth = 1, cp=cp)
        xvar.right <- ifelse(tree.right$frame[1, 'var']=='<leaf>', NA, as.character(tree.right$frame[1, 'var']))
        xval.right <- ifelse(is.null(tree.right$splits[1, 'index']), NA, tree.right$splits[1, 'index'])
        n.right <- ifelse(tree.right$frame[1, 'var']=='<leaf>', NA, tree.right$frame[1, 'n'])
        
        sse.left <- sse.fun(tree.left)
        sse.right <- sse.fun(tree.right)
        
        # treatment sum of squares...optimal split will maximize sst
        sst <- sstotal - (sse.left + sse.right)
        
        # assign split value to place holder column for numeric x splitting vals
        split_val_num <- x.unique[j]
        split_val_f <- NA
        
        df <- data.frame(
          x.name, split_val_num, split_val_f, nvals, xvar.left, xval.left, n.left
          ,xvar.right, xval.right, n.right, sst, stringsAsFactors=F
        )
        
        splits.df <- rbind(splits.df, df)
        splits.df <- unique(splits.df)
      }
    } else{   # x is categorical factor
      x.unique <- unique(x.temp)
      
      for(j in 1:length(x.unique)){
        x.left <- d.sort[x.temp!=x.unique[j], -length(d)]
        y.left <- y.temp[x.temp!=x.unique[j]]
        
        x.right <- d.sort[x.temp==x.unique[j], -length(d)]
        y.right <- y.temp[x.temp==x.unique[j]]
        
        # second stage tree building 
        tree.left <- rpart(y.left ~ ., data=x.left, minsplit=minsplit, maxdepth=1, cp=cp)
        xvar.left <- ifelse(tree.left$frame[1, 'var']=='<leaf>', NA, as.character(tree.left$frame[1, 'var']))
        xval.left <- ifelse(is.null(tree.left$splits[1, 'index']), NA, tree.left$splits[1, 'index'])
        n.left <- ifelse(tree.left$frame[1, 'var']=='<leaf>', NA, tree.left$frame[1, 'n'])
        
        tree.right <- rpart(y.right ~ ., data=x.right, minsplit=minsplit, maxdepth=1, cp=cp)
        xvar.right <- ifelse(tree.right$frame[1, 'var']=='<leaf>', NA, as.character(tree.right$frame[1, 'var']))
        xval.right <- ifelse(is.null(tree.right$splits[1, 'index']), NA, tree.right$splits[1, 'index'])
        n.right <- ifelse(tree.right$frame[1, 'var']=='<leaf>', NA, tree.right$frame[1, 'n'])
        
        sse.left <- sse.fun(tree.left)
        sse.right <- sse.fun(tree.right)
        
        # treatment sum of squares...optimal split will maximize sst
        sst <- sstotal - (sse.left + sse.right)
        
        # assign split value to place holder column for factor x splitting vals
        split_val_num <- NA
        split_val_f <- x.unique[j]
        
        df <- data.frame(
          x.name, split_val_num, split_val_f, nvals, xvar.left, xval.left, n.left
          ,xvar.right, xval.right, n.right, sst, stringsAsFactors=F
        )
        
        splits.df <- rbind(splits.df, df)
        splits.df <- unique(splits.df)
      }
    }
  }
  
  # clean up results from loop
  names(splits.df) <- c(
    'var', 'xval_num', 'xval_f', 'n', 'left_var', 'left_xval', 'left_n'
    ,'right_var', 'right_xval', 'right_n', 'sst'
  )
  
  # drops possible splits that do not have at least minsplit observations at 2nd level splitting
  rmv <- (splits.df[, 'left_n']<minsplit | is.na(splits.df[, 'left_n']) | splits.df[, 'right_n']<minsplit | is.na(splits.df[, 'right_n']))
  splits.df <- splits.df[!rmv, ]
  
  return(splits.df)
}


# ----- main tree building work function ANOVA method
build_tree.fun <- function(x.dat, y, nodenum=NULL, splitvar, splitvalnum, splitvalfac, splittxt, depth=0, minsplit, maxdepth, cp){
  # ++++++++++++++++++++++++++++++    function notes    ++++++++++++++++++++++++++++++
  # this function performs all the book keeping in regards to the current node including
  # checks to see if current portion of the tree is a terminal leaf, also calls the splitting
  # function and manages the general recursive structure/logic
  # ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  # initialize variables
  leaf.flag <- NULL
  d <- data.frame(cbind(x.dat, y))
  
  # key values for current node
  nvals <- nrow(d)
  nodemean <- mean(y)
  sstotal <- sstotal.fun(y) # total sum of squares for a node is also the sqaured residual component of model
  
  # determine if we are at the root node or in body of tree
  nodenum <- ifelse(is.null(nodenum), 1, nodenum)
  if(nodenum==1){
    splitvar <- NA
    splitvalnum <- NA
    splitvalfac <- NA
    splittxt <- 'root'
  }
  
  # check if current node is a leaf
  leaf.flag <- ifelse(nvals < minsplit | depth >= maxdepth, 1, 0)
  
  # assemble values to add to main tree structure
  node.df <- data.frame(
    nodenum, splitvar, splitvalnum, splitvalfac, splittxt, nvals, nodemean, sstotal, depth, leaf.flag, stringsAsFactors=F
  )
  print(node.df)
  
  # breaks function and no further splitting if we are at a leaf, otherwise find next split and continue
  if(leaf.flag==1){return(node.df)}
  
  # assess all potential next splitting move given results from all two stage splits
  splits.df <- splits.fun(d=d, sstotal=sstotal, nvals=nvals, minsplit=minsplit, cp=cp)
  
  # if there were no optimal splits (i.e. empty dataframe) then we are at a leaf, otherwise continue with best splitting results
  if(is.data.frame(splits.df) && nrow(splits.df)==0){
    node.df$leaf.flag <- 1
    return(node.df)
  }
  
  # find split that maximizes sum of squares treatment
  bestrow.indx <- which.max(splits.df[, 'sst'])
  bestrow <- splits.df[bestrow.indx, ]
  fac.flag <- is.na(bestrow$xval_num)   # T/F flag for if splitting covarite is a factor
  
  # use best split to partition data
  split.indx <- ifelse(
    rep(!fac.flag, nvals)
    ,x.dat[, names(x.dat)==bestrow$var] <= bestrow$xval_num
    ,as.character(x.dat[, names(x.dat)==bestrow$var]) == as.character(bestrow$xval_f)
  )
  
  x.leftson <- x.dat[split.indx, ]
  y.leftson <- y[split.indx]
  
  x.rightson <- x.dat[!split.indx, ]
  y.rightson <- y[!split.indx]
  
  # recursive call passing values to right and left sons...row append current node, left node, and right node results
  tree.df <- rbind(
    node.df, rbind(
      build_tree.fun(
        x.dat=x.leftson
        ,y=y.leftson
        ,nodenum=nodenum*2 + 0
        ,splitvar=bestrow$var
        ,splitvalnum=bestrow$xval_num
        ,splitvalfac=bestrow$xval_f
        ,splittxt=ifelse(!fac.flag, paste0(bestrow$var, '<=', bestrow$xval_num), paste0(bestrow$var, '=', bestrow$xval_f))
        ,depth=depth + 1
        ,minsplit=minsplit
        ,maxdepth=maxdepth
        ,cp=cp
      )   # close tree build call to left son
      ,build_tree.fun(
        x.dat=x.rightson
        ,y=y.rightson
        ,nodenum=nodenum*2 + 1
        ,splitvar=bestrow$var
        ,splitvalnum=bestrow$xval_num
        ,splitvalfac=bestrow$xval_f
        ,splittxt=ifelse(!fac.flag ,paste0(bestrow$var, '>', bestrow$xval_num) ,paste0(bestrow$var, '<>', bestrow$xval_f))
        ,depth=depth + 1
        ,minsplit=minsplit
        ,maxdepth=maxdepth
        ,cp=cp
      )   # close tree build call to right son
    )   # close inner rbind - rbind(leftson, rightson)
  )   # close outer rbind - rbind(current node, inner rbind)
  
  return(tree.df)
}


# ----- wrapper function to formalize call to tree building
two_stage_rpart <- function(formula, data, minsplit=10, maxdepth=30, cp=0.005){
  mf <- model.frame(formula, data)
  tree <- build_tree.fun(x.dat=mf[, -1], y=mf[, 1], minsplit=minsplit, maxdepth=maxdepth, cp=cp)
  names(tree) <- c(
    'node_num', 'split_var', 'split_val_num', 'split_val_f', 'split_txt', 'n', 'mean', 'sse', 'depth', 'leaf_flag'
  )
  tree <- tree[order(tree$node_num), ]
  tree$split_val_f <- as.factor(tree$split_val_f)
  return(tree)
}