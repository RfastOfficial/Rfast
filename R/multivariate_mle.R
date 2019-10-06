#[export]
invdir.mle <- function(x, tol = 1e-09) {

  n <- dim(x)[1]
  p <- dim(x)[2]
  zx <- t( log(x) )
  sx2 <- sum( log1p( Rfast::rowsums(x) ) )
  com <- c( Rfast::rowsums(zx) - sx2, -sx2 )
  a <- Rfast::colmeans(x)   
  b <- Rfast::colVars(x, suma = n * a)
  D <- p + 1 
  aD <- 0.5 * ( mean(a)^2 + mean(a) ) / mean(b) + 1
  a1 <- abs( c( a * (aD - 1), aD) ) / 2
  phi <- sum(a1)
  f1 <- n * digamma( phi ) - n * digamma(a1) + com
  f2 <- matrix(n * trigamma(phi), D, D)
  diag(f2) <- diag(f2) - n * trigamma(a1)
  a2 <- a1 - solve(f2, f1)
  i <- 2
  while ( sum( abs(a2 - a1) ) > tol ) {   
    i <- i + 1
    a1 <- a2
    phi <- sum(a1)
    f1 <- n * digamma( phi ) - n * digamma(a1) + com
    f2 <- matrix(n * trigamma(phi), D, D)
    diag(f2) <- diag(f2) - n * trigamma(a1)
    a2 <- a1 - solve(f2, f1)
  }
  
  phi <- sum(a2)
  lik <-  n * lgamma( phi ) - n * sum( lgamma(a2) ) +
      sum( zx * (a2[1:p] - 1) ) - phi * sx2   
  list(iters = i, loglik = lik, param = a2 )
} 


#[export]
diri.nr2 <- function(x, type = 1, tol = 1e-07) {

  dm <- dim(x)
  n <- dm[1]  ## the sample size
  p <- dm[2]  ## dimensionality

  if (type == 1) {

    m <- Rfast::colmeans(x)
    zx <- t( Log(x) )
    down <-  - sum( m * ( Rfast::rowmeans( zx ) - log(m) ) )
    sa <- 0.5 * (p - 1) / down  ## initial value for precision
    a1 <- sa * m  ## initial values
    gm <- rowsums(zx)
    z <- n * Digamma( sa )
    g <- z - n * Digamma(a1) + gm
    qk <-  - n * Trigamma(a1)
    b <- sum(g / qk) / ( 1/z - sum(1 / qk) )
    a2 <- a1 - (g - b)/qk
    while ( sum( abs( a2 - a1 ) ) > tol ) {
      a1 <- a2
      z <- n * digamma( sum(a1) )
      g <- z - n * Digamma(a1) + gm
      qk <-  - n * Trigamma(a1)
      b <- sum(g / qk) / ( 1/z - sum(1 / qk) )
      a2 <- a1 - (g - b) / qk
    }
    loglik <- n * Lgamma( sum(a2) ) - n * sum( Lgamma(a2) ) + sum( zx * (a2 - 1) )
    if ( is.null(colnames(x)) ) {
      names(a2) <- paste("X", 1:p, sep = "")
    } else  names(a2) <- colnames(x)
    res <- list(loglik = loglik, param = a2)

  } else {
    zx <- t( Log(x) )
    ma <- Rfast::rowmeans(zx)
    m <- Rfast::colmeans(x)
    down <- -sum(m * (ma - log(m)))
    sa <- 0.5 * (p - 1)/down
    a1 <- sa * m
    f <- ma - Digamma(a1) + digamma(sa)
    der <-  - Trigamma(a1) + trigamma(sa)
    a2 <- a1 - f/der
    a <- .Call(Rfast_diri_nr_type2, a1, a2, ma, p, tol)
    loglik <- n * lgamma(sum(a)) - n * sum(Lgamma(a)) + sum(zx * (a - 1))
    if (is.null(colnames(x))) {
      names(a) <- paste("X", 1:p, sep = "")
    } else names(a) <- colnames(x)
    res <- list(loglik = loglik, param = a)
  }

  res
}


#[export]
mvlnorm.mle <- function(x) {
  dm <- dim(x)
  d <- dm[2]
  n <- dm[1]
  y <- Rfast::Log(x)  ## transform the data to the whole of R^d
  m1 <- Rfast::colmeans(y)  ## mean vector of y
  sigma <- crossprod(y)/n - tcrossprod(m1)
  a <- n * d * log(2 * pi) + n * log(det(s)) + n * d - sum(y)
  
  s1 <- diag(sigma)
  m <- exp( m1 + 0.5 * s1 )  ## mean vector of x
  
  m2 <- outer(m1, m1, "+")
  s2 <- outer(s1, s1, "+")
  s <- exp( m2 + 0.5 * s2 ) * ( exp(sigma) - 1 ) 

  list(loglik = -0.5 * a, mu = m1, sigma = sigma, m = m, s = s)
}


#[export]
mvnorm.mle <- function(x) {
   m <- Rfast::colmeans(x)
   dm <- dim(x)
   n <- dim(x)[1]
   d <- dm[2]
   s <- crossprod(x)/n - tcrossprod(m)
   a <-  n * d * log(2 * pi) + n * log( det(s) ) + n * d
   list(loglik = - 0.5 * a, mu = m, sigma = s)
}


#[export]
mvt.mle <- function(x, v = 5, tol = 1e-07){
  ## x contains the data
  ## v is the degrees of freedom, set to 5 by default
  dm <- dim(x)
  p <- dm[2]   ;    n <- dm[1]  ## dimensions
  m <- Rfast::colmeans(x)  ## initial parameters
  y <- Rfast::eachrow(x, m, oper = "-")
  R <- crossprod(y)/(n - 1)
  y <- NULL
  if (v != 1 ) R <- abs( v - 1 ) / v  * R     
  con <- n * lgamma( (v + p)/2 ) - n * lgamma(v/2) - 0.5 * n * p * log(pi * v)
  ### step 1
  wi <- (v + p) / ( v + Rfast::mahala(x, m, R) )  ## weights
  y <- sqrt(wi) * ( Rfast::eachrow(x, m, oper = "-" ) )
  sumwi <- sum(wi)
  R <- crossprod(y) / sumwi   ## scatter estimate
  m <- Rfast::colsums(wi * x) / sumwi  ## location estimate
  dis <- Rfast::mahala(x, m, R)
  el1 <-  - n * log( det(R) ) - (v + p) * sum( log1p(dis/v) ) 
  ### step 2
  wi <- (v + p) / ( v + dis )  ## weights
  y <- sqrt(wi) * ( Rfast::eachrow(x, m, oper = "-" ) ) 
  sumwi <- sum(wi)
  R <- crossprod(y) / sumwi  ## scatter estimate 
  m <- Rfast::colsums(wi * x) / sumwi  ## location estimate
  dis <- Rfast::mahala(x, m, R)
  el2 <-  - n * log( det(R) ) - (v + p) * sum( log1p(dis/v) ) 
  ## Step 3 and above
  i <- 2
  while ( el2 - el1 > tol ) { ## 1e-06 is the tolerance level 
    ## between two successive values of the log-likelihood
    i <- i + 1
    el1 <- el2
    wi <- (v + p) / ( v + dis) ## updated weights
    y <- sqrt(wi) * ( Rfast::eachrow(x, m, oper = "-" ) ) 
    sumwi <- sum(wi)
    R <- crossprod(y) / sumwi  ## updated scatter estimate
    m <- Rfast::colsums(wi * x) / sumwi  ## updated location estimate
    dis <- Rfast::mahala(x, m, R)
    el2 <-  - n * log( det(R) )- (v + p) * sum( log1p(dis/v) )  
  }  ## updated log-likelihood 

  list(iters = i, loglik = 0.5 * el2 + con, location = m, scatter = R) 
}





