source("make-test-tree.R")
library(RUnit)

test.make.test.tree <- function() {
	generate_test_trees(callback = function(t) {
		checkTrue(is_dag(t))
		checkTrue(is_connected(t))
	})
}