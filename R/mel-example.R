library(igraph)
#library(RUnit)
#source("make-test-tree.R")

depth <- 7
nchildren <- 2:5
name_prefix <- "C"
mass_set <- 1:10
s_m_set <- 1:10 / 5
mar_cat_set <- c("A", "A", "A", "B", "B", "C")
use_cached_structures = TRUE
export_vertex_list = FALSE

hierarchical_number <- function(p_name, index) {
	sprintf("%s.%d", p_name, index)
}

print("create initial tree")
if (use_cached_structures) {
  pbs <- graph_from_edgelist(as.matrix(read.csv("el-ok.csv")))
} else {
  pbs0 <- make_initial_tree(hierarchical_number(name_prefix, 1))
  print("expand tree")
  pbs <- add_levels(pbs0, depth - 1, nchildren, hierarchical_number)
}
if (!is_tree(pbs)) stop("not a valid tree")
print(sprintf("%d components", length(V(pbs))))

print("create mass table")
leaves <- which(degree(pbs, mode="out") == 0)
non_leaves <- which(degree(pbs, mode="out") > 0)
if (use_cached_structures) {
  mt <- read.csv("mt-ok.csv")
} else {
  id <- names(V(pbs))
  m1 <- c(NA)
  m2 <- c(NA)
  s_m <- c(NA)
  mt <- data.frame(id, m1, m2, s_m, stringsAsFactors = FALSE)
  mt[leaves, "m1"] <- sample(mass_set, length(leaves), replace=TRUE)
  mt[leaves, "m2"] <- mt[leaves, "m1"]
  mt[leaves, "s_m"] <- sample(s_m_set, length(leaves), replace=TRUE)
  mt[leaves, "m3"] <- mt[leaves, "m1"]
  mt[non_leaves, "mar.cat"] <- sample(mar_cat_set, length(non_leaves), replace=TRUE)
}

print(sprintf("%d leaf components", length(leaves)))
total_mass <- sum(mt[leaves, "m1"])
print(sprintf("total mass %f", total_mass))

#
# update_m1(target, children, df)
#
# specialization of update_prop with "m1" and sum
#

update_m1 <- function(target, sources, df) {
	update_prop(df, df_get_by_id, df_set_by_id, target, sources, "m1",
	            function(x) rollup_sum(target, x))
}

#
# update_m2(target, children, df)
#
# specialization of update_prop with "m2" and sum_with_margin(.1, x)
#

margin <- 0.05
update_m2 <- function(target, children, df) {
	update_prop(df, df_get_by_id, df_set_by_id, target, children, "m2",
	            function(x) rollup_sum_with_margin(target, margin, x))
}

#
# update_s_m(target, children, df)
#
# specialization of update_prop with "s_m" and Frobenius norm (RSS)

update_s_m <- function(target, children, df) {
	update_prop(df, df_get_by_id, df_set_by_id, target, children, "s_m",
	            function(x) rollup_rss(target, x))
}

#
# update_m3(target, children, df)
#
# specialization of update_prop with "m3" and sum_with_margin_policy
#

policy_table <- data.frame(
  cat = c("A", "B", "C"),
  margin = c(.01, .02, .05)
)

margin_policy <- function(id, get, ds) {
  cat <- get(ds, id, "mar.cat")
  x <- policy_table[policy_table$cat == cat, "margin"]
}

update_m3 <- function(target, children, df) {
  update_prop(df, df_get_by_id, df_set_by_id, target, children, "m3",
              function(x) rollup_sum_with_margin_policy(target, margin_policy, df_get_by_id, df, x))
}

update_masses <- function(target, child, df) {
	df1 <- update_m1(target, child, df)
	df2 <- update_m2(target, child, df1)
	df3 <- update_s_m(target, child, df2)
	df4 <- update_m3(target, child, df3)
	df4
}

validate_mt <- function(tree, table) {
  validate_table(tree, table, df_get_ids, df_get_by_id, 'm1')
  validate_table(tree, table, df_get_ids, df_get_by_id, 'm2')
  validate_table(tree, table, df_get_ids, df_get_by_id, 's_m')
  validate_table(tree, table, df_get_ids, df_get_by_id, 'm3')
}

# damage pbs to check validation

pbs_disc <- pbs
pbs_disc <- delete_edges(pbs_disc, c(V(pbs_disc)[name = "C.1"], V(pbs_disc)[name = "C.1.1"])) # disconnected
pbs_cyc <- pbs
pbs_cyc <- add_edges(pbs_cyc, c(V(pbs_cyc)[name = "C.1.1"], V(pbs_cyc)[name = "C.1.2.1"])) # cyclic
pbs_mr <- pbs
pbs_mr <- add_edges(pbs_mr <- add_vertices(pbs_mr, 1, name = "C.2"),
                 c(V(pbs_mr)[name = "C.2"], V(pbs_mr)[name = "C.1.1"])) # multiple roots

# damage mt to check validation

mt_extra <- data.frame(rbind(as.matrix(mt), c("X", 0, 0, 0, 0, 0))) # extra entry in mt
mt_missing <- mt[1:(nrow(mt) - 1),] # missing entry in mt
mt_na <- mt; mt_na[nrow(mt_na),]$m1 <- NA # leaf without mass value

print("compute masses; starting timer")
s_time <- proc.time()
mt1 <- rollup(pbs, update_masses, mt, validate_mt)
print("produce mel")
mel <- mt1[order(mt1$id),]
print(mel[1:10,])
print(proc.time() - s_time)

print("update masses from leaf; starting timer")
leaf = names(leaves)[1]
mel2 <- mel
mel2 <- df_set_by_id(mel2, leaf, "m1", df_get_by_id(mel2, leaf, "m1") + 5)
mel2 <- df_set_by_id(mel2, leaf, "m2", df_get_by_id(mel2, leaf, "m2") + 5)
mel2 <- df_set_by_id(mel2, leaf, "m3", df_get_by_id(mel2, leaf, "m3") + 5)
s_time <- proc.time()
mel2 <- update_rollup(pbs, V(pbs)[leaf], update_masses, mel2)
print(mel2[1:10,])
print(proc.time() - s_time)

if (export_vertex_list) {
  children <- sapply(names(V(pbs)), FUN=function(x) paste(names(neighbors(pbs, x, "out")), collapse=","))
  vl <- data.frame(id=names(V(pbs)), level=as.vector(distances(pbs, v=1)),children)
  leaves <- which(!is.na(mt$m1))
  vl[leaves,'m1'] <- mt[leaves,'m1']
  write.csv(vl[order(vl$id),], "vl.csv")
}