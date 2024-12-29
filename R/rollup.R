#
# rollup(tree, update, ds, validate_ds)
#
#     tree:         igraph dag
#     ds:           dataset to update
#     update:       called depth-first at each node as update(ds, parent_key, child_keys)
#     validate_ds:  called as validate_ds(tree, ds)

rollup <- function(tree, ds, update, validate_ds, validate_tree = default_validate_tree) {
  root <- validate_tree(tree)
  validate_ds(tree, ds)
  Reduce(
    f = function(s, v) update(s, names(V(tree))[v], names(neighbors(tree, v, "in"))),
    x = dfs(tree, root, mode="in", order=FALSE, order.out=TRUE)$order.out,
    init = ds
  )
}

# 
# update_rollup(tree, vertex, update, ds)
#
#     tree:         igraph dag
#     vs:           vertex to update from
#     update:       called for successive ancestors of vertex as update(ds, parent_key, child_keys)
#     ds:           dataset to update

update_rollup <- function(tree, vertex, update, ds) {
  if (degree(tree, vertex, mode="out") > 0) stop("update_rollup() on non-leaf")
  Reduce(
    f = function(s, v) update(s, names(V(tree))[v], names(neighbors(tree, v, "in"))),
    x = na.omit(as.vector(dfs(tree, vertex, mode="out", unreachable=FALSE, order=TRUE)$order))[-1],
    init = ds
  )
}

#
# default_validate_tree(tree)
#

default_validate_tree <- function(tree) {
  if (girth(tree, circle = FALSE)$girth != Inf) stop("graph is cyclic") # girth() changed in igraph 1.6.0
  if (any(which_loop(tree))) stop("graph contains loops")
  if (any(which_multiple(tree))) stop("graph contains multiple edges")
  if (!is_connected(tree)) stop("graph is disconnected")
  if (!is_directed(tree)) stop("graph is undirected")
  roots <- which(degree(tree, mode = "out") == 0)
  if (length(roots) > 1) stop("graph contains multiple roots")
  roots[1]
}

#
# validate_table(tree, table)
#

validate_table <- function(tree, table, get_keys, get_prop, op=function(x) is.numeric(x)) {
  tree_ids <- names(V(tree))
  table_ids <- get_keys(table)
  if (!setequal(tree_ids, table_ids)) stop("mismatched ids")
  leaves <- names(which(degree(tree, mode = "in") == 0))
  if (any(sapply(leaves, FUN=function(l) !op(get_prop(table, l)))))
    stop (paste("leaf with invalid value"))
}
