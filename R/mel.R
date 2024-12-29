#
# rollup_sum(id, x)

rollup_sum <- function(id, x) sum(x)

#
# rollup_sum_with_margin(id, m, x)
#
# sums x and applies margin m
# example combining operator for update_prop
#

rollup_sum_with_margin <- function(id, m, x) (1 + m) * sum(x)

#
# rollup_rss(id, x)
#
# calculates RSS of x

rollup_rss <- function(id, x) norm(as.matrix(x), type='F')

#
# rollup_sum_with_margin_policy

rollup_sum_with_margin_policy <- function(id, p, get, ds, x) (1 + p(id, get, ds)) * sum(x)