library(igraph, quietly=TRUE, warn.conflict=FALSE)
library(huxtable)
library(parallel)

source("../R/rollup.R")
source("../R/update_prop.R")
source("../R/huxtable_from_tree_table.R")
source("../R/mass-props.R")

# load sbs table and edges

sbs_table <- read.csv("../mt-names.csv", stringsAsFactors = FALSE)
sbs_table <- sbs_table[order(sbs_table$key), c("id", "name")]
sbs_edges <- read.csv("../el-ok.csv", stringsAsFactors = FALSE)

# sort sbs table in id order

# construct graph from edges; add any vertices named in the table
# but not part of an edge

sbs <- graph_from_edgelist(as.matrix(sbs_edges[, c("child", "parent")]), directed = TRUE)
sbs <- sbs + vertices(setdiff(sbs_table$id, names(V(sbs))))

# create huxtable from sbs and sbs_table

ht <- huxtable_from_tree_table(sbs, sbs_table)

# add leaf masses

leaves <- names(V(sbs)[degree(sbs, mode="in") == 0])
sbs_table[is.element(sbs_table$id, leaves), "mass"] <- runif(length(leaves), min=1.0, max=3.0)

# configure helper functions for rollup

update_mass <- function(parent_id, child_ids, df) {
  update_df_prop_by_id(df, parent_id, child_ids, "mass")
}
validate_mass <- function(tree, df) {
  validate_df_with_id(tree, df, "mass")
}

# roll up masss

sbs_table_out1 <- rollup(sbs, update_mass, sbs_table, validate_mass)

# combine indented sbs and mass numbers

mht <- cbind(ht, hux(sbs_table$mass))
mht <- set_number_format(mht, everywhere, ncol(mht), "%8.2f")

# add other mass properties

saved_mass <- sbs_table[is.element(sbs_table$id, leaves), "mass"]
sbs_table <- read.csv("../mt-names.csv", stringsAsFactors = FALSE)
sbs_table <- sbs_table[order(sbs_table$key), c("id", "name")]
nr <- nrow(sbs_table)
sbs_table <- cbind(
  sbs_table,
  mass = NA,
  Cx = NA,
  Cy = NA,
  Cz = NA,
  Ixx = NA,
  Iyy = NA,
  Izz = NA,
  Ixy = NA,
  Ixz = NA,
  Iyz = NA,
  Ipoint = NA,
  Iconv = NA
)

# https://physics.stackexchange.com/questions/348944/what-is-the-problem-of-having-an-inertia-tensor-not-satisfying-the-triangle-ineq

random_mass_props <- function() {
  mp <- list()
  
  mp$mass <- runif(1, min=0.001, max=1.0)
  
  mp$centroid <- runif(3, min=-1.0, max=1.0)
  
  # principal moments of a rectangular prism of random dimensions
  
  pm <- as.vector(
    matrix(
      c(0, 1, 1,
        1, 0, 1,
        1, 1, 0), nrow=3) %*% runif(3, min=.001, max=1.0)
    )
  
  # rotate through random x, y, z angles
  
  angles <- runif(3, min=-pi, max=pi)
  R <- rotation_from_angles(angles)
  
  mp$inertia_tensor <- list(
    point= FALSE,
    matrix = R %*% diag(pm) %*% t(R),
    conv = '-'
  )
  
  # random products of inertia sign convention
  
  as_it_conv(mp, sample(c('+', '-'), 1))
}

for (id in leaves) {
  sbs_table <- df_set_mass_props(sbs_table, id, random_mass_props())
}
sbs_table[is.element(sbs_table$id, leaves), "mass"] <- saved_mass

# microbenchmark::microbenchmark(rollup(sbs, update_mass_props, sbs_table, validate_mass_props_table), times=50)
# Unit: seconds
#                                                                  expr      min       lq     mean  median       uq      max neval
#  rollup(sbs, update_mass_props, sbs_table, validate_mass_props_table) 2.668145 2.866308 3.176349 3.08931 3.433663 4.010564    50

sbs_table_out2 <- rollup(sbs, update_mass_props, sbs_table, validate_mass_props_table)

# compare with example from https://github.jpl.nasa.gov/Mass-Properties/excel_mass_properties_tools/blob/master/simplified_MP_sheet.xlsx
# uses positive integral convention
# note: test data fails inertia tensor validation so we call with a vacuous validator

test_table <- read.table("mass-props.txt", header=TRUE, colClasses=c("character", "character", rep("numeric", 10)))
test_edges <- read.table("mp_edges.txt", header=TRUE, colClasses=c("character", "character"))
test_tree <- graph_from_edgelist(as.matrix(test_edges[, c("pid", "cid")]), directed=TRUE)
test_table_out <- rollup(test_tree,
       function(parent_key, child_keys, ds) {
         update_prop(ds, parent_key, child_keys, df_set_mass_props, df_get_mass_props,
                     get_l = NULL, op = function(...) combine_mass_props(..., it_conv='+'))
         },
       test_table,
       function(tree, df) TRUE)

test_ref <- as.numeric(read.csv("simplified_MP_sheet.csv", header=FALSE, skip=9)[1, 11:20])
r <- round(test_table_out[1, 3:12], 2) - test_ref  # Div 35 spreadsheet rounds to 2 places
r_norm <- norm(as.matrix(r), "F")

# monte carlo simulation

# monte_carlo_sample_rollup
# take n samples of a top-level rollup by perturbing all leaves and performing the rollup
#
# n:        number of samples
# update:   property update method as used in rollup()
# tree:     tree as used in rollup()
# table:    property table as used in rollup()
# get:      property getter as used in rollup()
# set:      property setter as used in rollup()
# perturb   property perturbation method, called as perturb(property), returns perturbed property
# row:      row number of top-level element in table
# cores:    number of cores to use for parallel execution

monte_carlo_sample_rollup <- function(n, update, tree, table, get, set, perturb, row, cores=2) {
  
  # cache loop invariants
  
  t <- cbind(table) # copy
  leaves <- names(V(tree)[degree(tree, mode="in") == 0])
  
  # roll up n perturbations in parallel, collect specified rows
  
  do.call(rbind, mclapply(
    1:n,
    FUN = function(i) {
      for (leaf in leaves) {
        t <- set(t, leaf, perturb(get(t, leaf)))
      }
      
      rollup(tree, update, t, function(t, s) TRUE)[row, ]
    }, mc.cores=cores
  ))
  
}

s_10 <- monte_carlo_sample_rollup(
  10,
  update_mass_props,
  sbs, 
  sbs_table,
  df_get_mass_props,
  df_set_mass_props,
  function(mp) perturb_mass_props(mp,
                                  mass_sd=0.1, mass_min=0.001,
                                  Cx_sd=0.1, Cy_sd=0.1, Cz_sd=0.1,
                                  x_sd=0.1, x_min=0.001,
                                  y_sd=0.1, y_min=0.001,
                                  z_sd=0.1, z_min=0.001,
                                  psi_sd=0.1, theta_sd=0.1, phi_sd=0.1),
  1,
  detectCores())

# threats and opportunities

overrides <- read.csv("overrides.csv", strip.white=TRUE, stringsAsFactors=FALSE)

# override_mass_props
#
# ol:       overrides list with elements m and b
# ds:       dataset for input parameters (not used in this example)
# id:       element id in dataset (not used in this example)
# mp:       mass properties to be overridden

override_mass_props <- function(ol, ds, id, mp) {
  
  # get initial values
  
  mass <- mp$mass
  it <- mp$inertia_tensor
  
  # apply overrides
  
  for (row in which(ol$id == id)) {
    ov <- ol[row, ]
    new_mass <- mass * ov$m + ov$b
    if (ov$scale_it & mass > 0.0) 
      if (!it$point) 
        it$matrix <- (new_mass / mass) * it$matrix
    mass <- new_mass
  }
  
  list(mass=mass, centroid=mp$centroid, inertia_tensor=it, it_conv=mp$it_conv)
}

# rollup with one override

overrides3 <- overrides[overrides$name == "Element 1 Descope",]
sbs_table_out3 <-
  rollup(sbs,
         function(parent_key, child_keys, ds)
           update_mass_props(
             parent_key,
             child_keys,
             ds,
             override =
               function(ds, id, mp)
                 override_mass_props(overrides3, ds, id, mp)
           ),
         sbs_table,
         validate_mass_props_table)

# rollup with all assembly overrides

overrides4 <- overrides[grepl("Assembly", overrides$name), ]
sbs_table_out4 <-
  rollup(sbs,
         function(parent_key, child_keys, ds)
           update_mass_props(
             parent_key,
             child_keys,
             ds,
             override =
               function(ds, id, mp)
                 override_mass_props(overrides4, ds, id, mp)
           ),
         sbs_table,
         validate_mass_props_table)
