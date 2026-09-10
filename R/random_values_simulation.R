# [export]
matrnorm <- function(n, p, seed = NULL) {
  .Deprecated(new = "Rnorm,mat(n, p)", msg = "Seed argument in package rangen can be changed via `setSeed`",package = "rangen")
  if ( !is.null(seed) ) zigg::zsetseed(seed)
  matrix( zigg::zrnorm(n * p), ncol = p )
}


# [export]
racg <- function(n, sigma, seed = NULL) {
  .Deprecated(new = "racg(n, sigma, seed)", msg = "Just moved to another package.",package = "rangen")
  ## n is the sample size,
  ## mu is the mean vector and
  ## sigma is the covariance matrix
  ## sigma does not have to be of full rank
  p <- dim(sigma)[1]
  if ( !is.null(seed) ) zigg::zsetseed(seed)
  x <- rangen::Rnorm.mat(n, p)
  x <- x %*% chol(sigma)
  x / sqrt( Rfast::rowsums(x^2) )
}


# [export]
rbing <- function(n, lam) {
  .Deprecated(new = "rbing(n, lam)", msg = "Just moved to another package.",package = "rangen")
  .Call(Rfast_rbing, n, lam)
}


######### Simulation using any symmetric A matrix
# [export]
rbingham <- function(n, A) {
  .Deprecated(new = "rbingham(n, A)", msg = "Just moved to another package.",package = "rangen")
  p <- dim(A)[2] ## dimensionality of A
  eig <- eigen(A)
  lam <- eig$values ## eigenvalues
  V <- eig$vectors ## eigenvectors
  lam <- lam - lam[p]
  lam <- lam[-p]
  x <- Rfast::rbing(n, lam)
  ## the x contains the simulated values
  tcrossprod(x, V) ## simulated data
}


# [export]
rmvlaplace <- function(n, lam, mu, G, seed = NULL) {
  .Deprecated(new = "rmvlaplace(n, lam, mu, G, seed)", msg = "Just moved to another package.",package = "rangen")
  ## n is the sample size
  ## lam is the parameter of the exponential distribution
  ## m is the mean vector
  ## G is a d x d covariance matrix with determinant 1
  if (summary(det(G))[1] == 1) {
    y <- paste("The determinant of the covariance matrix is not 1.")
  } else {
    d <- length(mu) ## dimensionality of the data
    z <- rexp(n, lam)
    if (!is.null(seed)) zigg::zsetseed(seed)
    x <- rangen::Rnorm.mat(n, d)
    y <- sqrt(z) * x %*% chol(G) + rep(mu, rep(n, d)) ## the simulated sample
  }
  y
}


# [export]
rmvnorm <- function(n, mu, sigma, seed = NULL) {
  .Deprecated(new = "rmvnorm(n, mu, sigma, seed)", msg = "Just moved to another package.",package = "rangen")
  p <- length(mu)
  if (!is.null(seed)) zigg::zsetseed(seed)
  x <- rangen::Rnorm.mat(n, p)
  x %*% chol(sigma) + rep(mu, rep(n, p))
}


# [export]
rmvt <- function(n, mu, sigma, v, seed = NULL) {
  .Deprecated(new = "rmvt(n, mu, sigma, v, seed)", msg = "Just moved to another package.",package = "rangen")
  ## n is the sample size
  ## mu is the mean vector
  ## sigma is the covariance matrix
  ## sigma does not have to be of full rank
  ## v is the degrees of freedom
  p <- length(mu)
  if (!is.null(seed)) zigg::zsetseed(seed)
  x <- rangen::Rnorm.mat(n, p)
  w <- sqrt(v / rchisq(n, v))
  w * x %*% chol(sigma) + rep(mu, rep(n, p))
}


# [export]
Rnorm <- function(n, m = 0, s = 1, seed = NULL) {
  .Deprecated(new = "Rnorm(n, m, s)", msg = "Seed argument in package rangen can be changed via `setSeed`",package = "rangen")
  if (!is.null(seed)) zigg::zsetseed(seed)
  if (m == 0 & s == 1) {
    x <- zigg::zrnorm(n)
  } else if (m == 0 & s != 1) {
    x <- zigg::zrnorm(n) * s
  } else if (m != 0 & s == 1) {
    x <- zigg::zrnorm(n) + m
  } else {
    x <- zigg::zrnorm(n) * s + m
  }
  x
}


# [export]
rvmf <- function(n, mu, k, parallel = FALSE, cores = 0, seed = rangen::new.Seed()) {
  # rotation <- function(a, b) {
    # p <- length(a)
    # ab <- sum(a * b)
    # ca <- a - b * ab
    # ca <- ca / sqrt(sum(ca^2))
    # A <- b %*% t(ca)
    # A <- A - t(A)
    # theta <- acos(ab)
    # diag(p) + sin(theta) * A + (cos(theta) - 1) * (b %*%
      # t(b) + ca %*% t(ca))
  # }
  # d <- length(mu)
  # if (k > 0) {
    # mu <- mu / sqrt(sum(mu^2))
    # ini <- c(numeric(d - 1), 1)
    # d1 <- d - 1
    # v1 <- rangen::Rnorm.mat(n, d1) ##  matrix( zigg::zrnorm(n * d1), ncol = d1 )
    # v <- v1 / sqrt(Rfast::rowsums(v1^2))
    # b <- (-2 * k + sqrt(4 * k^2 + d1^2)) / d1
    # x0 <- (1 - b) / (1 + b)
    # m <- 0.5 * d1
    # ca <- k * x0 + (d - 1) * log(1 - x0^2)
    # w <- as.vector(rvmf_h(n, ca, d1, x0, m, k, b))
    # S <- cbind(sqrt(1 - w^2) * v, w)
    # if (isTRUE(all.equal(ini, mu, check.attributes = FALSE))) {
      # x <- S
    # } else if (isTRUE(all.equal(-ini, mu, check.attributes = FALSE))) {
      # x <- -S
    # } else {
      # A <- rotation_R(ini, mu)
      # x <- tcrossprod(S, A)
    # }
  # } else {
    # x1 <- rangen::Rnorm.mat(n, d) ## matrix( zigg::zrnorm(n * d), ncol = d )
    # x <- x1 / sqrt(Rfast::rowsums(x1^2))
  # }
  .Call(Rfast_rvmf,n,mu,k,parallel, cores, seed)
}


# rvmf <- function(n, mu, k) {
#   ## n is the sample size
#   ## mu is the mean direction and
#   ## k is the concentration parameter
#   ## n is the sample size
#   d <- length(mu)  ## the dimensions
#   if (k > 0) {
#     mu <- mu / sqrt( sum(mu^2) )  ## the mean direction
#     ini <- c(numeric(d - 1), 1)  ## mean direction is set to (0, ..., 0, 1)
#     d1 <- d - 1
#     v1 <- matrix( zigg::zrnorm(n * d1), ncol = d1)
#     v <- v1 / sqrt( rowsums(v1^2) )
#     b <- ( -2 * k + sqrt(4 * k^2 + d1^2) ) / d1
#     x0 <- (1 - b)/(1 + b)
#     w <- numeric(n)
#     m <- 0.5 * d1
#     ca <- k * x0 + (d - 1) * log(1 - x0^2)
#
#     for (i in 1:n) {
#       ta <-  -1000
#       u <- 1
#       while ( ta - ca < log(u) ) {
#         z <- rbeta(1, m, m)
#         u <- runif(1)
#         w[i] <- ( 1 - (1 + b) * z ) / ( 1 - (1 - b) * z )
#         ta <- k * w[i] + d1 * log(1 - x0 * w[i])
#       }
#     }
#     S <- cbind(v, w)
#     A <- rotation(ini, mu)  ## calculate the rotation matrix
#     ## in order to rotate the initial mean direction from ini to mu
#     x <- tcrossprod(S, A)  ## the x has direction mu
#   } else {  ## uniform distribution
#     ## requires MASS if k = 0
#     x1 <- matrix( zigg::zrnorm(n * d), ncol = d )
#     x <- x1 / sqrt( Rfast::rowsums(x1^2) )
#   }
#   x
# }


# [export]
rvonmises <- function(n, m, k, rads = TRUE, parallel = FALSE, cores = 0, seed = rangen::new.Seed()) {
  #  if (!rads) m <- m / 180 * pi ## turn the degrees into radians
  #  mu <- c(cos(m), sin(m))
  #  if (k > 0) { ## draw from a von Mises distribution
  #    x <- Rfast::rvmf(n, mu, k) ## sample from the von Mises distribution
  #    u <- (atan(x[, 2] / x[, 1]) + 
  #    pi * I(x[, 1] < 0))
  #    u <- u %% (2 * pi) ## u is in radians
  #  } else {
  #    u <- runif(n, 0, 2 * pi)
  #  } ## draw from a von Mises distribution
  #  if (!rads) u <- u * pi / 180 ## should the data be in degrees?
  #  u
  .Call(Rfast_rvonmises, n,m,k,rads,parallel, cores, seed)
}
