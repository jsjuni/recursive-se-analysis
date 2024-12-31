library(tidyverse)

source("R/mass-props.R")

el <- as.matrix(read_csv("el-ok.csv", show_col_types = FALSE))[, c("child", "parent")]
g <- graph_from_edgelist(el)
leaf_ids <- names(V(g)[degree(g, mode="in") == 0])

pt <- read.csv("mt-names.csv")
leaves <- which(is.element(pt$id, leaf_ids))
pt$POIconv <- "-"
pt[leaves, "POIconv"] <- sample(c("-", "+"), length(leaves), replace = TRUE)

Rxyz <- function(a, b, c) {
  sina <- sin(a)
  cosa <- cos(a)
  sinb <- sin(b)
  cosb <- cos(b)
  sinc <- sin(c)
  cosc <- cos(c)
  matrix(c(
    cosa * cosb, cosa * sinb * sinc - sina * cosc, cosa * sinb * cosc + sina * sinc,
    sina * cosb, sina * sinb * sinc + cosa * cosc, sina * sinb * cosc - cosa * sinc,
    -sinb, cosb * sinc, cosb * cosc
  ), nrow=3, byrow=TRUE)
}

random_mass_props <- function() {
  mp <- list()
  
  mp$mass <- runif(1, min=0.001, max=1.0)
  
  mp$center_mass <- runif(3, min=-100.0, max=100.0)
  
  # principal moments of a cuboid of specified mass and random dimensions
  
  xyz <- list("x", "y", "z")
  pm <- mp$mass / 12. * as.vector(matrix(c(0, 1, 1,
                                           1, 0, 1,
                                           1, 1, 0), nrow=3) %*% runif(3, min=5.0, max=10.0))
  
  # rotate through random x, y, z angles
  
  angle <- runif(3, min=-pi, max=pi)
  R <- Rxyz(angle[1], angle[2], angle[3])
  mp$inertia <- R %*% diag(pm) %*% t(R)
  dimnames(mp$inertia) <- list(xyz, xyz)
  
  mp$point <- FALSE

  mp
}

for (id in leaf_ids) {
  pt <- df_set_mass_props(pt, id, random_mass_props())
}

pt[leaves, "σ_mass"] <- sd(pt[leaves, "mass"])

pt[leaves, "σ_Cx"] <- sd(abs(pt[leaves, "Cx"]))
pt[leaves, "σ_Cy"] <- sd(abs(pt[leaves, "Cy"]))
pt[leaves, "σ_Cz"] <- sd(abs(pt[leaves, "Cz"]))

pt[leaves, "σ_Ixx"] <- sd(pt[leaves, "Ixx"])
pt[leaves, "σ_Iyy"] <- sd(pt[leaves, "Iyy"])
pt[leaves, "σ_Izz"] <- sd(pt[leaves, "Izz"])
pt[leaves, "σ_Ixy"] <- sd(abs(pt[leaves, "Ixy"]))
pt[leaves, "σ_Ixz"] <- sd(abs(pt[leaves, "Ixz"]))
pt[leaves, "σ_Iyz"] <- sd(abs(pt[leaves, "Iyz"]))

write_tsv(pt, file = "mp-input.tsv")

# test to make sure it's ok
rs <- rollup(g, pt, update_mass_props_and_unc, validate_mass_props_table)

