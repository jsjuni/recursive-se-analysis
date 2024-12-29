library(igraph)

make_initial_tree <- function(name) {
    graph <- make_empty_graph()
    add_vertices(graph, 1, name = name)
}

add_levels <- function(graph, depth, nchild, child_name) {
    for (level in 1:depth) {
        graph <- add_level(graph, nchild, child_name)
    }
    graph
}

add_level <- function(graph, nchild, child_name) {
    dv <- degree(graph, mode = "out")
    leaves <- dv[dv == 0]
    for (p_name in names(leaves)) {
        graph <- add_children(graph, p_name, nchild, child_name)
    }
    graph
}

add_children <- function(graph, p_name, nchild, child_name) {
    for (ci in 1:sample(nchild, 1)) {
        graph <- add_child(graph, p_name, ci, child_name)
    }
    graph
}

add_child <- function(graph, p_name, index, child_name) {
    c_name <- child_name(p_name, index)
    graph <- add_vertices(graph, 1, name = c_name)
    vs = V(graph)
    graph <- add_edges(graph, c(vs[name = p_name], vs[name = c_name]))
    graph
}

is_acyclic <- function(graph) {
    as.integer(girth(graph)[1]) == 0
}

is_tree <- function(graph) {
    is_directed(graph) && is_acyclic(graph) && is_connected(graph)
}

nodes <- function(k, h) if (k == 1) h else (k^(h+1) - 1) / (k - 1)

generate_test_trees <- function(klimit = 10, hlimit = 20, nlimit = 10000, callback) {
	f0 <- data.frame(expand.grid(1:klimit, 1:hlimit))
	n <- apply(f0, 1, FUN = function(x) nodes(x[1], x[2] + 1))
	f1 <- cbind(f0, n)
	colnames(f1) <- c('k', 'h', 'n')
	f <- f[f$n < nlimit,]
	for (r in 1:nrow(f)) {
		t <- add_levels(make_initial_tree('C'), f[r, 'h'], f[r, 'k'], function(n, i) sprintf("%s.%d", n, i))
		callback(t)
	}
}
