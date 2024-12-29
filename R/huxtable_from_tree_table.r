library(huxtable)

huxtable_from_tree_table <- function(tree, tree_table, use.names=TRUE, spacer=NULL) {
  
  # validate tree
  
  root <- validate_tree(tree)
  
  # display column is tree depth + 1
  
  tree_table$column <- dfs(tree, root, mode="in", order=FALSE, dist=TRUE)$dist[tree_table$id] + 1
  
  # create matrix for huxtable data
  
  width <- max(tree_table$column)
  m <- matrix(nrow=nrow(tree_table), ncol=width)
  for (c in 1:width) {
    rows = which(tree_table$column == c)
    
    if (!is.null(spacer))
      m[rows, 1:(c-1)] = spacer
    
    m[rows, c] <- if (use.names)
      paste(tree_table[rows, "id"], tree_table[rows, "name"], sep=" ")
    else
      tree_table[rows, "id"]
    
  }
  
  # create huxtable
  
  ht <- hux(m)
  
  # set column spans
  
  for (c in 1:(width - 1)) {
    rows <- which(tree_table$column == c)
    ht <- set_colspan(ht, rows, c, width - c + 1)
  }
  
  ht
}

