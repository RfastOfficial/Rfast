#[export]
beta.mle <- function(x, tol = 1e-09) {
  n <- length(x)
  sly1 <- sum( log(x) ) / n
  sly2 <- sum( log(1 - x) ) / n
  sy <- sum(x)
  sy2 <- sum(x^2)
  iniphi <- (sy - sy2)/(sy2 - sy^2/n) * (n - 1)/n
  a <- sy * iniphi/n
  b <- iniphi - a
  phi <- a + b 
  lik1 <-  - n * lbeta(a, b) + (a - 1) * n * sly1 + (b - 1) * sly2 * n
  dera <- sly1 - digamma(a) + digamma(phi)
  derb <- sly2 - digamma(b) + digamma(phi)
  derab <- trigamma(phi)
  dera2 <-  - trigamma(a) + derab
  derb2 <-  - trigamma(b) + derab 
  anew <- c(a, b) - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
  a <- anew[1]     ;   b   <- anew[2]
  phi <- a + b
  lik2 <-  - n * lbeta(a, b) + (a - 1) * n * sly1 + (b - 1) * sly2 * n
  
  i <- 2
  while (lik2 - lik1 > tol) {
    i <- i + 1
	lik1 <- lik2
    dera <- sly1 - digamma(a) + digamma(phi)
    derb <- sly2 - digamma(b) + digamma(phi)
    derab <- trigamma(phi)
    dera2 <-  - trigamma(a) + derab
    derb2 <-  - trigamma(b) + derab 
    anew <- anew - c(derb2 * dera - derab * derb, -derab * dera + dera2 * derb)/(dera2 * derb2 - derab^2)  
	a <- anew[1]    ;   b <- anew[2]
    phi <- a + b
    lik2 <-  - n * lbeta(a, b) + (a - 1) * n * sly1 + (b - 1) * sly2 * n
  }
  loglik <-  lik2
  names(anew) <- c("alpha", "beta")
  list(iters = i, loglik = loglik, param = anew)
}

  
#[export]
hsecant01.mle <- function(x, tol = 1e-09) {
  
  sy1 <- sum( log(x) )  ;   sy2 <-  sum( log( 1 - x) )
  com <-  - 0.5 * sy1 - 0.5 * sy2
  comp <- sy1 / pi - sy2 / pi
  n <- length(x)
  a <-   - 0.5 * pi   ;   b <- 0.5 * pi
  ratio <- 2 / (sqrt(5) + 1)
  x1 <- b - ratio * (b - a)
  x2 <- a + ratio * (b - a)
  f1 <-  n * log( cos(x1) ) + x1 * comp
  f2 <-  n * log( cos(x2) ) + x2 * comp

  while (abs(b - a) > tol) {
    if (f2 < f1) {
      b <- x2
      x2 <- x1
      f2 <- f1
      x1 <- b - ratio * (b - a)
      f1 <- n * log( cos(x1) ) + x1 * comp
    } else {
      a <- x1
      x1 <- x2
      f1 <- f2
      x2 <- a + ratio * (b - a)
      f2 <- n * log( cos(x2) ) + x2 * comp
    }
  }
  theta <- (x1 + x2) / 2
  loglik <- n * log( cos(theta) / pi ) + com + theta * comp
  list(loglik = loglik, theta = theta)
}

#hsecant01 <- function(y, tol = 1e-09) {
#  sy1 <- sum( log(y) )  ;   sy2 <-  sum( log( 1 - y) )
#  com <-  - 0.5 * sy1 - 0.5 * sy2
#  comp <- sy1 / pi - sy2 / pi
#  n <- length(y)
#  f <-  function(theta, n, comp)  n * log( cos(theta)/pi ) + theta * comp
#  mod <- optimise(f, c(-pi/2, pi/2), n = n, comp = comp, maximum = TRUE, tol = tol)
#  loglik <- mod$objective + com
#  list(loglik = mod$objective + com, theta = mod$maximum)
#}
  

#[export]
ibeta.mle <- function(x, tol = 1e-09) {
 
  if ( all(x > 0  &  x < 1) ) {   
    res <- Rfast::beta.mle(x)
    mes <- "Regular beta distribution was fitted."

  } else {
    n <- length(x)
    z <- x[ x > 0 & x < 1 ]
    T2 <- sum( log(z) )
    T3 <- sum( log(1 - z) ) 
    if ( min(x) == 0 ) {
      T1 <- sum( x == 0 )
      mes <- "Zero inflated beta was fitted."
    } else {
      T1 <- sum( x == 1 )
      mes <- "One inflated beta was fitted."
    }
    a <- T1 / n
    t23 <- T2 - T3
    f <- n - T1
    ini <- Rfast::beta.mle(z)
    phi1 <- sum(ini$param)
    m1 <- ini$param[1]/phi1
    m1phi <- m1 * phi1
    m2phi <- (1 - m1) * phi1
    derm <- phi1 * f * ( digamma(m2phi) - digamma( m1phi ) ) + phi1 * t23 
    derphi <- f * ( digamma(phi1) - m1 * digamma(m1phi) - (1 - m1) * digamma(m2phi) ) + m1 * t23  + T3
    derm2 <-  - f * phi1^2 * ( trigamma(m1phi) + trigamma(m2phi) )
    derphi2 <- f * ( trigamma(phi1) - trigamma(m1phi) * m1^2 - trigamma(m2phi) * (1 - m1)^2 )
    dermphi <- f * ( - trigamma(m1phi) * m1phi - digamma(m1phi) + trigamma(m2phi) * m2phi + digamma(m2phi) ) + t23
    aold <- c(m1, phi1)
    anew <- aold - c( derphi2 * derm - dermphi * derphi, - dermphi * derm + derm2 * derphi ) / ( derm2 * derphi2 - dermphi^2 )
    i <- 2
    while ( sum( abs(aold - anew) ) > tol ) {
      i <- i + 1
      m1 <- anew[1]    ;   phi1 <- anew[2]
      aold <- anew
      m1phi <- m1 * phi1
      m2phi <- (1 - m1) * phi1
      derm <- phi1 * f * ( digamma(m2phi) - digamma( m1phi ) ) + phi1 * t23 
      derphi <- f * ( digamma(phi1) - m1 * digamma(m1phi) - (1 - m1) * digamma(m2phi) ) + m1 * t23  + T3
      derm2 <-  - f * phi1^2 * ( trigamma(m1phi) + trigamma(m2phi) )
      derphi2 <- f * ( trigamma(phi1) - trigamma(m1phi) * m1^2 - trigamma(m2phi) * (1 - m1)^2 )
      dermphi <- f * ( - trigamma(m1phi) * m1phi - digamma(m1phi) + trigamma(m2phi) * m2phi + digamma(m2phi) ) + t23
      anew <- aold - c( derphi2 * derm - dermphi * derphi, - dermphi * derm + derm2 * derphi ) / ( derm2 * derphi2 - dermphi^2 )
    }
    m <- anew[1]   ;   phi <- anew[2]
    param <- c(a, m, phi1)
    names(param) <- c("Propoprtion", "mean", "precision")
    loglik <- T1 * log(a) + f * log(1 - a) + f * ( lgamma(phi) - lgamma(m * phi) - lgamma( (1 - m) * phi ) ) + T2 * (m * phi - 1) + T3 * ( (1 - m) * phi - 1 )
    res <- list(iters = i, loglik = loglik, param = param)
  } 
  
  list(mes <- mes, res <- res)
}


#[export]
logitnorm.mle <- function(x) {
  n <- length(x)
  lx1 <- log(x)
  lx2 <- log(1 - x)
  y <- lx1 - lx2 
  sy <- sum(y)
  m <- sy / n
  s <- ( sum(y^2) - n * m^2 ) / n
  loglik <- sum( dnorm(y, m, sqrt(s), log = TRUE ) ) - sy
  param <- c(m, n * s / (n - 1) )
  names(param) <- c("mean", "unbiased variance")
  list(loglik = loglik, param = param)
}




