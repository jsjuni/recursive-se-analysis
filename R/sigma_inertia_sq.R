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

sigma_inertia_sq_1 <- function(sigma_inertia, m, sigma_m, d, sigma_cm) {
  it0 = matrix(nrow = 3, ncol = 3, dimnames = list(c('x', 'y', 'z'), c('x', 'y', 'z')))
  
  # moments
  
  it1 <- Reduce(
    f = function(t, ss) {
      t[ss[1], ss[1]] = sigma_moment_sq(sigma_inertia[ss[1], ss[1]], m, sigma_m, d[ss[2:3]], sigma_c_m[ss[2:3]])
      t
    },
    x = list(c('x', 'y', 'z'), c('y', 'x', 'z'), c('z', 'x', 'y')),
    init = it0
  )
  
  # products
  
  Reduce(
    f = function(t, ss) {
      t[ss[1], ss[2]] = t[ss[2], ss[1]] = sigma_product_sq(sigma_inertia[ss[1], ss[2]], m, sigma_m, d[ss[1:2]], sigma_c_m[ss[1:2]])
      t
    },
    x = list(c('x', 'y'), c('x', 'z'), c('y', 'z')),
    init = it1
  )
}

sigma_inertia_sq_2 <- function(sigma_inertia, m, sigma_m, d, sigma_c_m) {

  P <- outer(d, sigma_c_m)
  p <- diag(P)
  diag_1 <- diag(c(p['x'] + 2 * p['y'], p['y'] + 2 * p['x'], p['z'] + 2 * p['x']))
  diag_2 <- diag(c(p['x'] + 2 * p['z'], p['y'] + 2 * p['z'], p['z'] + 2 * p['y']))
  
  Q <- outer(d, d)
  diag_3 <- sum(diag(Q)) * diag(3)

  sigma_inertia^2 + (m * (P - diag_1))^2 + (m * (t(P) - diag_2))^2 + (sigma_m * (Q - diag_3))^2
}

sigma_inertia_sq_3 <- function(sigma_inertia, m, sigma_m, d, sigma_c_m) {
  
  # intermediate computations
  
  si <- c(sigma_inertia[1,1], sigma_inertia[2,2], sigma_inertia[3,3], sigma_inertia[1,2], sigma_inertia[1,3], sigma_inertia[2,3])
  P <- (m * outer(d, sigma_c_m)) ^ 2     # $(w_i (x_i - \bar{x}) \sigma_{x_i})^2$, ...
  Q <- sigma_m * outer(d, d)             # $\sigma_{w_i} (x_i - \bar{x}) (y_i - \bar{y})$, ...
  
  si ^ 2 +
    c(
      
      # moments
      
      4 * c(
        P['y', 'y'] + P['z', 'z'],        # Ixx
        P['x', 'x'] + P['z', 'z'],        # Iyy
        P['x', 'x'] + P['y', 'y']         # Izz
      ) + c(
        Q['y', 'y'] + Q['z', 'z'],        # Ixx
        Q['x', 'x'] + Q['z', 'z'],        # Iyy
        Q['x', 'x'] + Q['y', 'y']         # Izz
      ) ^ 2,
      
      # products
      
      c(
        P['x', 'y']  + P['y', 'x'],       # Ixy
        P['x', 'z']  + P['z', 'x'],       # Ixz
        P['y', 'z']  + P['z', 'y']        # Iyz
      ) + c(
        Q['x', 'y'],                      # Ixy
        Q['x', 'z'],                      # Ixz
        Q['y', 'z']                       # Iyz
      ) ^ 2
    )
}

# test data

m <- 2
sigma_m <- 3
d <- c(x = 4, y = 5, z = 6)
sigma_c_m <- c(x = 7, y = 8, z = 9)
sigma_inertia <- matrix(
  nrow = 3,
  ncol = 3,
  data = c(
     10, 11, 12,
     11, 13, 14,
     12, 14, 15
  ),
  dimnames = list(c('x', 'y', 'z'), c('x', 'y', 'z'))
)

s1 <- sigma_inertia_sq_1(sigma_inertia, m, sigma_m, d, sigma_c_m)
s2 <- sigma_inertia_sq_2(sigma_inertia, m, sigma_m, d, sigma_c_m)
s3 <- sigma_inertia_sq_3(sigma_inertia, m, sigma_m, d, sigma_c_m)
