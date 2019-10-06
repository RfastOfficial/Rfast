#[export]
cauchy.mle <- function (x, tol = 1e-09) {
    n <- length(x)
    m <- Rfast::Median(x)
    es <- 0.5 * ( Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4) )
    logs <- log(es)
    y <- x - m
    y2 <- y^2
    lik1 <- n * logs - sum(log(es^2 + y2))
    down <- 1/(es^2 + y2)
    down2 <- down^2
    derm <- 2 * sum(y * down)
    ders <- n - 2 * es^2 * sum(down)
    derm2 <- 2 * sum( (y2 - es^2) * down2 )
    ders2 <-  - 2 * es^2 * ( derm2 + 2 * es^2 * sum( down2 ) )
    derms <-  - 4 * es^2 * sum(y * down2)
    m <- m - ( ders2 * derm - derms * ders ) / (derm2 * ders2 - derms^2)
    logs <- logs - ( -derms * derm + derm2 * ders ) / (derm2 * ders2 - derms^2)
    y <- x - m
    y2 <- y^2
    es <- exp(logs)
    lik2 <- n * logs - sum( log(es^2 + y2) )
    i <- 2
    while (lik2 - lik1 > tol) {
        i <- i + 1
        lik1 <- lik2
        down <- 1/(es^2 + y2)
        down2 <- down^2
        derm <- 2 * sum(y * down)
        ders <- n - 2 * es^2 * sum(down)
        derm2 <- 2 * sum( (y2 - es^2) * down2 )
        ders2 <-  - 2 * es^2 * ( derm2 + 2 * es^2 * sum( down2 ) )
        derms <-  - 4 * es^2 * sum(y * down2)
        m <- m - ( ders2 * derm - derms * ders ) / (derm2 * ders2 - derms^2)
        logs <- logs - ( -derms * derm + derm2 * ders ) / (derm2 * ders2 - derms^2)
        y <- x - m
        y2 <- y^2
        es <- exp(logs)
        lik2 <- n * logs - sum( log(es^2 + y2) )
    }
    param <- c(m, es)
    names(param) <- c("location", "scale")
    list(iters = i, loglik = lik2 - n * log(pi), param = param)
}


#[export]
ct.mle <- function(x, tol = 1e-09) {

  n <- length(x)
  f <- 0.5 * n
  f2 <- 0.25 * n 
  x2 <- x^2
  s <- sum(x2)/n
  v <- 2 * s / abs(s - 1)
  y <- x2 / v
  lik1 <- n * lgamma(0.5 * v + 0.5) - n * 0.5 * log(v * pi) - n * lgamma(0.5 * v) - 
            0.5 * (v + 1) * sum( log1p(y) )
  ra <- 1 + y
  der <- f * digamma(0.5 * v + 0.5) - f / v - f * digamma(0.5 * v) -
       0.5 * sum( log(ra) ) + (v + 1) / 2 * sum( y / ra) / v
  der2 <- f2 * trigamma(0.5 * v + 0.5) + f / v^2 - f2 * trigamma(0.5 * v) + 
        sum( y / ra ) / v - 0.5 * (v + 1) * sum( ( 2 * y + y^2) / (v + x2 )^2 )   
  v <- v - der / der2
  y <- x2 / v
  lik2 <- n * lgamma(0.5 * v + 0.5) - n * 0.5 * log(v * pi) - n * lgamma(0.5 * v) - 
            0.5 * (v + 1) * sum( log1p(y) )
  i <- 2 

  while ( lik2 - lik1 > tol  &  v > 0) {
    i <- i + 1 
    lik1 <- lik2 
    ra <- 1 + y
    der <- f * digamma(0.5 * v + 0.5) - f / v - f * digamma(0.5 * v) -
          0.5 * sum( log(ra) ) + (v + 1) / 2 * sum( y / ra) / v
    der2 <- f2 * trigamma(0.5 * v + 0.5) + f / v^2 - f2 * trigamma(0.5 * v) + 
          sum( y / ra ) / v - 0.5 * (v + 1) * sum( ( 2 * y + y^2) / (v + x2 )^2 )   
    v <- v - der / der2
    y <- x2 / v
	lik2 <- n * lgamma(0.5 * v + 0.5) - n * 0.5 * log(v * pi) - n * lgamma(0.5 * v) - 
            0.5 * (v + 1) * sum( log1p(y) )
  }   
  list(iters = i, nu = v, loglik = lik2)
}  


#[export]
gumbel.mle <- function(x, tol = 1e-09) {

  n <- length(x) 
  m <- sum(x) / n
  x2 <- x^2
  s1 <- sqrt( 6 * sum(x2) / n - 6 * m^2 ) / pi
  y <- exp( - (x - m) / s1 )
  sy <- sum(y)
  co <- sum(x * y)
  f <- s1 - m + co / sy
  f2 <- 1 + ( sum(x2 * y) * sy - co^2 ) / s1^2 / sy^2
  s2 <- s1 - f / f2  
  i <- 2
  while ( abs(s1 - s2) > tol  )  {
    i <- i + 1
    s1 <- s2
    y <- exp( - (x - m) / s1 )
    sy <- sum(y)
    co <- sum(x * y)
    f <- s1 - m + co / sy
    f2 <- 1 + ( sum(x2 * y) * sy - co^2 ) / s1^2 / sy^2
    s2 <- s1 - f / f2
  } 
    
  y2 <- exp( - x / s2 )
  sy2 <- sum(y2)
  m <-  - s2 * log( sy2 / n )
  y <- (x - m) / s2
  lik <-  - n * log(s2) - sum(y) - sum( exp(-y) )
  param <- c(m, s2)
  names(param) <- c("location", "scale")
  list(iters = i, loglik = lik, param = param)
}


#[export]
laplace.mle <- function(x) {
  n <- length(x)
  m <- Rfast::Median(x)
  b <- sum( abs(x - m) ) / n
  param <- c(m, b)
  names(param) <- c("location", "scale")
  loglik <-  - n * log(2 * b) - n 
  list(loglik = loglik, param = param)
}


#[export]
logistic.mle <- function(x, tol = 1e-07) {
 
  n <- length(x)
  m <- sum(x) / n    
  exps <- sqrt(3) * sd(x) / pi
  y <- ( x - m ) / exps
  y2 <- 0.5 * y
  tanhy2 <- tanh(y2)
  sechy22 <- 1 / cosh(y2)^2
  derm <- sum( tanhy2 ) / exps
  derm2 <-  - 0.5 * sum( sechy22 ) / exps^2 
  ders <-  - n + sum( y * tanhy2 )
  ders2 <-  - ders - n - sum( y2^2 * sechy22 )
  derms <-  - derm - sum( sechy22 * y2) / exps
  aold <- c(m, log(exps) )
  anew <- aold - c( ders2 * derm - derms * ders, - derms * derm + derm2 * ders) / (derm2 * ders2 - derms^2)
  
  i <- 2
  while ( sum( abs(aold - anew) ) > tol ) {
    i <- i + 1
    aold <- anew  
    m <- anew[1]    
    exps <- exp( anew[2] ) 
    y <- (x - m ) / exps
    y2 <- 0.5 * y
    tanhy2 <- tanh(y2)
    sechy22 <- 1 / cosh(y2)^2
    derm <- sum( tanhy2 ) / exps
    derm2 <-  - 0.5 * sum( sechy22 ) / exps^2 
    ders <-  - n + sum( y * tanhy2 )
    ders2 <-  - ders - n - sum( y2^2 * sechy22 )
    derms <-  - derm - sum( sechy22 * y2) / exps
    anew <- aold - c( ders2 * derm - derms * ders, - derms * derm + derm2 * ders) / (derm2 * ders2 - derms^2)
  }
  
  anew[2] <- exp(anew[2])
  names(anew) <- c("location", "scale")
  lik <-  - n * log(4 * exps) + sum( log(sechy22) )
  list(iters = i, param = anew, loglik = lik)
}  


#[export]
normal.mle <- function(x) {
    n <- length(x)
    m <- sum(x)/n
    s <- (sum(x^2) - n * m^2)/(n - 1)
    loglik <-  - 0.5 * n * ( log(2 * pi) + log(s) ) - 0.5 * (n - 1)
    param <- c(m, s)
    names(param) <- c("mean", "unbiased variance")
    list(loglik = loglik, param = param)
}


#[export]
tmle <- function(x, v = 5, tol = 1e-08) {
  
  n <- length(x)   ;    m <- sum(x) / n
  s <- ( sum(x^2) - n * m ) / ( n - 1)
   ## m and s are the initial parameters
  f <- n / 2
  up <- (v + 1) / 2
  con <- n * lgamma( up ) - n * lgamma(v/2) - f * log(pi * v)
  ### step 1
  y <- (x - m)^2
  wi <- 1 / ( v + y / s )  ## weights
  sumwi <- sum(wi) 
  s <-  sum(wi * y) / sumwi  ## scatter estimate
  m1 <- sum(wi * x) / sumwi  ## location estimate
  y <- (x - m1)^2
  z <- v + y / s
  ### step 2  
  wi <- 1 / z  ## weights
  sumwi <- sum(wi) 
  s <-  sum(wi * y) / sumwi  ## scatter estimate
  m2 <- sum(wi * x) / sumwi  ## location estimate
  y <- (x - m2)^2
  z <- v + y / s
  ## Step 3 and above
  i <- 2

  while ( abs(m1 - m2) > tol ) { ## 1e-06 is the tolerance level 
    ## between two successive values of the log-likelihood
    i <- i + 1
    m1 <- m2
    wi <- 1 / z  ## weights
    sumwi <- sum(wi) 
    s <-  sum(wi * y) / sumwi  ## scatter estimate
    m2 <- sum(wi * x) / sumwi  ## location estimate
    y <- (x - m2)^2
    z <- v + y / s
  }  
  
  loglik <-  - f * log( s ) - up * sum( log( z / v ) )
  param <- c(m2, s)
  names(param) <- c("location", "scatter")
  list(iters = i, loglik = loglik + con, param = param) 
}


#[export]
wigner.mle <- function(x, tol = 1e-09) {
  n <- length(x)
  r <- max( abs(x ) )
  down <- r^2 - x^2
  list(loglik = n * log(2 / pi / r^2 ) + 0.5 * sum( log(down[down != 0]) ), R = r^2 )
}




