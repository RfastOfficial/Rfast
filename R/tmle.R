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
  list(iters = i, loglik = loglik + con,  param = param) 
}





