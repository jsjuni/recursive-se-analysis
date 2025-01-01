suppressPackageStartupMessages(library(magrittr))

# getter for mass props aggregate

df_get_mass_props <- function(df, id) {
  poi_conv <- df_get_by_id(df, id, "POIconv")
  list(
    mass = df_get_by_id(df, id, "mass"),
    center_mass = sapply(c(x = "Cx", y = "Cy", z = "Cz"), FUN=function(p) df_get_by_id(df, id, p)),
    inertia = {
      xyz <- list("x", "y", "z")
      it <- matrix(data = rep.int(0, 9), nrow = 3, dimnames = list(xyz, xyz))
      it["x", "x"] <- df_get_by_id(df, id, "Ixx")
      it["y", "y"] <- df_get_by_id(df, id, "Iyy")
      it["z", "z"] <- df_get_by_id(df, id, "Izz")
      poi_factor <- if (poi_conv == '-') 1 else -1
      it["x", "y"] <- it["y", "x"] <- poi_factor * df_get_by_id(df, id, "Ixy")
      it["x", "z"] <- it["z", "x"] <- poi_factor * df_get_by_id(df, id, "Ixz")
      it["y", "z"] <- it["z", "y"] <- poi_factor * df_get_by_id(df, id, "Iyz")
      it
    },
    point = df_get_by_id(df, id, "Ipoint")
  )
}

df_get_mass_props_and_unc <- function(df, id) {
  r <- df_get_mass_props(df, id)
  r$σ_mass <- df_get_by_id(df, id, "σ_mass")
  r$σ_center_mass <- sapply(c(x = "σ_Cx", y = "σ_Cy", z = "σ_Cz"), FUN=function(p) df_get_by_id(df, id, p))
  r$σ_it <- {
    xyz <- list("x", "y", "z")
    σ_it <- matrix(data = rep.int(0, 9), nrow = 3, dimnames = list(xyz, xyz))
    σ_it["x", "x"] <- df_get_by_id(df, id, "σ_Ixx")
    σ_it["y", "y"] <- df_get_by_id(df, id, "σ_Iyy")
    σ_it["z", "z"] <- df_get_by_id(df, id, "σ_Izz")
    σ_it["x", "y"] <- σ_it["y", "x"] <- df_get_by_id(df, id, "σ_Ixy")
    σ_it["x", "z"] <- σ_it["z", "x"] <- df_get_by_id(df, id, "σ_Ixz")
    σ_it["y", "z"] <- σ_it["z", "y"] <- df_get_by_id(df, id, "σ_Iyz")
    σ_it
  }
  r
}

# setter for mass props aggregate

df_set_mass_props <- function(df, id, v) {
  m <- v$inertia
  poi_factor <- if (v$poi_conv == "-") 1 else -1
  df %>% df_set_by_id(id, "mass", v$mass) %>%
  
    df_set_by_id(id, "Cx", v$center_mass[1]) %>%
    df_set_by_id(id, "Cy", v$center_mass[2]) %>%
    df_set_by_id(id, "Cz", v$center_mass[3]) %>%
    
    df_set_by_id(id, "Ixx", m["x", "x"]) %>%
    df_set_by_id(id, "Iyy", m["y", "y"]) %>%
    df_set_by_id(id, "Izz", m["z", "z"]) %>%
    df_set_by_id(id, "Ixy", poi_factor * (m["x", "y"] + m["y", "x"]) / 2.0) %>%
    df_set_by_id(id, "Ixz", poi_factor * (m["x", "z"] + m["z", "x"]) / 2.0) %>%
    df_set_by_id(id, "Iyz", poi_factor * (m["y", "z"] + m["z", "y"]) / 2.0) %>%

    df_set_by_id(id, "POIconv", v$poi_conv) %>%
    df_set_by_id(id, "Ipoint", v$point)
}

df_set_mass_props_and_unc <- function(df, id, v) {
  df %>% df_set_mass_props(id, v) %>%
    
    df_set_by_id(id, "σ_mass", v$σ_mass) %>%
    
    df_set_by_id(id, "σ_Cx", v$σ_center_mass[1]) %>%
    df_set_by_id(id, "σ_Cy", v$σ_center_mass[2]) %>%
    df_set_by_id(id, "σ_Cz", v$σ_center_mass[3]) %>%
    
    df_set_by_id(id, "σ_Ixx", v$σ_it["x", "x"]) %>%
    df_set_by_id(id, "σ_Iyy", v$σ_it["y", "y"]) %>%
    df_set_by_id(id, "σ_Izz", v$σ_it["z", "z"]) %>%
    df_set_by_id(id, "σ_Ixy", v$σ_it["x", "y"]) %>%
    df_set_by_id(id, "σ_Ixz", v$σ_it["x", "z"]) %>%
    df_set_by_id(id, "σ_Iyz", v$σ_it["y", "z"])
}

# mass properties combiner

combine_mass_props <- function(vl) {
  
  r <- list()
  
  # sum of masses
  
  r$mass <- Reduce(`+`, Map(f = function(v) v$mass, vl))
  
  # mass-weighted sum of centers of mass
  
  r$center_mass <- Reduce(`+`, Map(f = function(v) v$mass * v$center_mass, vl)) / r$mass 
  
  # parallel axis theorem
  # https://en.wikipedia.org/wiki/Parallel_axis_theorem#Moment_of_inertia_matrix
  # d_ss2 is [d]^2 computed using the identities given
  r$inertia <- Reduce(
    `+`,
    Map(
      f  = function(v) {
        d <- r$center_mass - v$center_mass
        ddt <- outer(d, d)
        d_ss2 <- ddt - sum(diag(ddt)) * diag(3)
        if (v$point) -v$mass * d_ss2 else v$inertia - v$mass * d_ss2
      },
      vl
    )
  )
  
  # aggregate is a point mass iff all parts are point masses at the same center
  
  r$point <- all(vl$point) && isTRUE(all.equal(r$center_mass, vl$center_mass))
  
  r
}

combine_mass_props_and_unc <- function(vl) {
  
  r <- combine_mass_props(vl)
  
  # mass uncertainty
  
  r$σ_mass = sqrt(Reduce(`+`, Map(f = function(v) v$σ_mass^2, vl)))
  
  # center of mass uncertainty
  
  r$σ_center_mass = sqrt(Reduce(`+`, Map(
    f = function(v) {
      (v$mass * v$σ_center_mass) ^ 2 +
        (v$σ_mass * (v$center_mass - r$center_mass)) ^ 2
    },
    vl
  ))) / r$mass
  
  # inertia tensor uncertainty
  
  r$σ_it = sqrt(Reduce(`+`, Map(
    f = function(v) {
      
      d <- r$center_mass - v$center_mass
      
      P <- outer(d, v$σ_center_mass)
      p <- diag(P)
      diag_1 <- diag(c(p['x'] + 2 * p['y'], p['y'] + 2 * p['x'], p['z'] + 2 * p['x']))
      diag_2 <- diag(c(p['x'] + 2 * p['z'], p['y'] + 2 * p['z'], p['z'] + 2 * p['y']))
      
      Q <- outer(d, d)
      diag_3 <- sum(diag(Q)) * diag(3)
      
      v$σ_it^2 + (v$mass * (P - diag_1))^2 + (v$mass * (t(P) - diag_2))^2 + (v$σ_mass * (Q - diag_3))^2
    },
    vl
  )))
  
  # result
  
  r
}

df_override_mass_props <- function(ds, id, v) {
  v$poi_conv <- df_get_by_id(ds, id, "POIconv")
  v
}

sigma_moment_sq <- function(sigma_moment, m, sigma_m, d, sigma_c_m) {
  sigma_moment^2 +
    (2 * m * d[1] * sigma_c_m[1])^2 + 
    (2 * m * d[2] * sigma_c_m[2])^2 +
    ((d[1]^2 + d[2]^2) * sigma_m)^2
}

sigma_product_sq <- function(sigma_product, m, sigma_m, d, sigma_c_m) {
  sigma_product^2 +
    (d[1] * m * sigma_c_m[2])^2 +
    (d[1] * d[2] * sigma_m)^2 +
    (d[2] * m * sigma_c_m[1])^2
  
}

combine_mass_props_and_uncX <- function(vl) {
  
  r <- combine_mass_props(vl)
  
  # mass uncertainty
  
  r$σ_mass = sqrt(Reduce(`+`, Map(f = function(v) v$σ_mass^2, vl)))
  
  # center of mass uncertainty
  
  r$σ_center_mass = sqrt(Reduce(`+`, Map(
    f = function(v) {
      (v$mass * v$σ_center_mass) ^ 2 +
        (v$σ_mass * (v$center_mass - r$center_mass)) ^ 2
    },
    vl
  ))) / r$mass
  
  # inertia tensor uncertainty
  
  xyz <- list("x", "y", "z")
  r$σ_it = sqrt(Reduce(`+`, Map(
    f = function(v) {
      it0 = matrix(nrow = 3, ncol = 3, dimnames = list(xyz, xyz))
      d <- v$center_mass - r$center_mass
      
      # moments
      
      it1 <- Reduce(
        f = function(t, ss) {
          t[ss[1], ss[1]] = sigma_moment_sq(v$σ_it[ss[1], ss[1]], v$mass, v$σ_mass, d[ss[2:3]], v$σ_center_mass[ss[2:3]])
          t
        },
        x = list(c("x", "y", "z"), c("y", "x", "z"), c("z", "x", "y")),
        init = it0
      )
      
      # products
      
      Reduce(
        f = function(t, ss) {
          t[ss[1], ss[2]] = t[ss[2], ss[1]] = sigma_product_sq(v$σ_it[ss[1], ss[2]], v$mass, v$σ_mass, d[ss[1:2]],  v$σ_center_mass[ss[1:2]])
          t
        },
        x = list(c("x", "y"), c("x", "z"), c("y", "z")),
        init = it1
      )
    },
    vl
  )))
  
  # result
  
  r
}

# mass properties updater

update_mass_props <- function(ds, parent_key, child_keys) {
  update_prop(
    ds,
    target = parent_key,
    sources = child_keys,
    set = df_set_mass_props,
    get = df_get_mass_props,
    combine = combine_mass_props,
    override = df_override_mass_props
  )
}

update_mass_props_and_unc <- function(ds, parent_key, child_keys) {
  update_prop(
    ds,
    target = parent_key,
    sources = child_keys,
    set = df_set_mass_props_and_unc,
    get = df_get_mass_props_and_unc,
    combine = combine_mass_props_and_unc,
    override = df_override_mass_props
  )
}

# validate mass property set

validate_mass_props <- function(mp) {
  
  # ensure mass is numeric and positive.

  if (is.na(mp$mass)) stop("mass missing")
  if (!is.numeric(mp$mass)) stop("mass non-numeric")
  if (mp$mass <= 0.) stop("mass non-positive")
  
  # ensure center of mass is numeric.
  
  if (any(is.na(mp$center_mass))) stop("center of mass element missing")
  if (any(!is.numeric(mp$center_mass))) stop("center of mass element non-numeric")
  
  # ensure inertia tensor point mass indicator is logical
  
  if (!is.logical(mp$point)) stop("invalid inertia tensor point mass indicator")
  
  if (!mp$point) {
    
    # ensure inertia tensor elements for non-point-masses are numeric.
    
    if (any(is.na(mp$inertia))) stop("inertia tensor element missing")
    if (any(!is.numeric(mp$inertia))) stop("inertia tensor element non-numeric")
    
    # ensure inertia tensor is positive definite.
    
    ev <- eigen(mp$inertia, symmetric=TRUE, only.values=TRUE)$values
    if (any(ev <= 0.)) stop("inertia tensor not positive definite")
    
    # ensure principal moments obey triangle inequalities
    
    if (any(c(
      ev[1] >= ev[2] + ev[3],
      ev[2] >= ev[1] + ev[3],
      ev[3] >= ev[1] + ev[2]
    ))) stop("inertia tensor violates triangle inequalities")
    
  }
  
  TRUE
}

# validate mass properties table

validate_mass_props_table <- function(tree, df) {
  validate_table(tree, df, df_get_ids, df_get_mass_props, validate_mass_props)
}

# perturb mass with a sample of mean zero and given sd, clamped below at min

perturb_mass <- function(mass, sd=0.0, min=0.0) {
  max(rnorm(1, mass, sd), min)
}

# perturb center of mass with three samples of mean zero and given sds

perturb_center_mass <- function(center_mass, x_sd=0.0, y_sd=0.0, z_sd=0.0) {
  center_mass + c(rnorm(1, 0.0, x_sd), rnorm(1, 0.0, y_sd), rnorm(1, 0.0, z_sd))
}

# equivalent unit mass cuboid dimensions and rotation

cuboid_from_inertia_tensor <- {
  minv <- matrix(c(-1, 1, 1,
                   1,-1, 1,
                   1, 1,-1), byrow = TRUE, nrow = 3) * 6
  
  function(it) {
    eg <- eigen(it, symmetric=TRUE, only.values = FALSE)
    
    list(dims = sqrt(as.vector(minv %*% eg$values)),
         rotation = eg$vector)
  }
}

# inertia tensor of unit mass cuboid with specified dimensions and rotation

inertia_tensor_from_cuboid <- {
  m <- matrix(c(0, 1, 1,
                1, 0, 1,
                1, 1, 0), byrow = TRUE, nrow = 3) / 12
  xyz <- list("x", "y", "z")
  
  function(dims, R) {
    i <- R %*% diag(as.vector(m %*% matrix(dims^2, ncol = 1))) %*% t(R)
    dimnames(i) <- list(xyz, xyz)
    i
  }
}

# After Arvo 1991
# https://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
#
#   *  R A N D _ R O T A T I O N      Author: Jim Arvo, 1991                  *
#   *                                                                         *
#   *  This routine maps three values (x[0], x[1], x[2]) in the range [0,1]   *
#   *  into a 3x3 rotation matrix, M.  Uniformly distributed random variables *
#   *  x0, x1, and x2 create uniformly distributed random rotation matrices.  *
#   *  To create small uniformly distributed "perturbations", supply          *
#   *  samples in the following ranges                                        *
#   *                                                                         *
#   *      x[0] in [ 0, d ]                                                   *
#   *      x[1] in [ 0, 1 ]                                                   *
#   *      x[2] in [ 0, d ]                                                   *
#   *                                                                         *
#   * where 0 < d < 1 controls the size of the perturbation.  Any of the      *
#   * random variables may be stratified (or "jittered") for a slightly more  *
#   * even distribution.                                                      *

random_rotation <- function(x0, x1, x2) {
  theta <- 2.0 * pi * x0
  phi   <- 2.0 * pi * x1
  z     <- 2.0 * x2
  
  r <- sqrt(z)
  Vx <- sin(phi) * r
  Vy <- cos(phi) * r
  Vz <- sqrt(2.0 - z)
  
  st <- sin(theta)
  ct <- cos(theta)
  Sx <- Vx * ct - Vy * st
  Sy <- Vx * st + Vy * ct
  
  matrix(
    data = c(
      Vx * Sx - ct, Vx * Sy - st, Vx * Vz,
      Vy * Sx + st, Vy * Sy - ct, Vy * Vz,
      Vz * Sx,      Vz * Sy,      1.0 - z
    ),
    nrow = 3,
    byrow = TRUE
  )
}

# perturb inertia tensor

perturb_inertia_tensor <- function(it,
                                   x_sd=0.0, x_min=0.0, y_sd=0.0, y_min=0.0, z_sd=0.0, z_min=0.0,
                                   r_sd = 0.0) {
  # get equivalent cuboid
  
  cuboid <- cuboid_from_inertia_tensor(it)
  dims <- cuboid$dims
  R <- cuboid$rotation
  
  # perturb cuboid dimensions
  
  dims_p <- c(max(rnorm(1, dims[1], x_sd), x_min),
              max(rnorm(1, dims[2], y_sd), y_min),
              max(rnorm(1, dims[3], z_sd), z_min))
  
  # perturb cuboid rotation
  
  r3 = runif(3) * c(r_sd, 1, r_sd)
  R_p <- R %*% random_rotation(r3[1], r3[2], r3[3])
  
  # construct tensor
  
  inertia_tensor_from_cuboid(dims_p, R_p)
}

# perturb mass properties

perturb_mass_props <- function(mp,
                               mass_sd=0.0, mass_min=0.0,
                               Cx_sd=0.0, Cy_sd=0.0, Cz_sd=0.0,
                               x_sd=0.0, x_min=0.0, y_sd=0.0, y_min=0.0, z_sd=0.0, z_min=0.0,
                               r_sd=0.0) {
  list(
    mass = perturb_mass(mp$mass, mass_sd, mass_min),
    center_mass = perturb_center_mass(mp$center_mass, Cx_sd, Cy_sd, Cz_sd),
    inertia = perturb_inertia_tensor(mp$inertia,
                                     x_sd, x_min, y_sd, y_min, z_sd, z_min,
                                     r_sd),
    point = mp$point
  )
}
