#[export]
betaprime.mle <- function(x, tol = 1e-09) {

  n <- length(x)
  slx <- sum( log(x) )
  slx2 <- sum( log1p(x) )
  m <- sum(x) / n    ;   s <- sum(x^2) - m^2
  b <- ( m^2 + m ) / s + 2 
  a <- abs( m * b - m )
  lik1 <- (a - 1) * slx - (a + b) * slx2 - n * lbeta(a, b)
  
  dera <- n * digamma(a + b) - n * digamma(a) + slx - slx2
  derb <- n * digamma(a + b) - n * digamma(b) - slx2 
  derab <- n * trigamma(a + b)
  dera2 <- derab - n * trigamma(a)
  derb2 <- derab - n * trigamma(b)
  anew <- c(a, b) - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  a <- anew[1]      ;      b <- anew[2]
  lik2 <- (a - 1) * slx - (a + b) * slx2 - n * lbeta(a, b)
  
  i <- 2
  while ( lik2 - lik1 > tol ) {
    i <- i + 1
	lik1 <- lik2
    a <- anew[1]      ;      b <- anew[2]
    dera <- n * digamma(a + b) - n * digamma(a) + slx - slx2
    derb <- n * digamma(a + b) - n * digamma(b) - slx2 
    derab <- n * trigamma(a + b)
    dera2 <- derab - n * trigamma(a)
    derb2 <- derab - n * trigamma(b)
    anew <- anew - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
    a <- anew[1]      ;      b <- anew[2]
    lik2 <- (a - 1) * slx - (a + b) * slx2 - n * lbeta(a, b)
  }  
    
   names(anew) <- c("alpha", "beta")
   loglik <- (a - 1) * slx - (a + b) * slx2 - n * lbeta(a, b)
   list(iters = i, loglik = loglik, param = anew)   
}


#[export]
chisq.mle <- function(x, tol = 1e-09) {

  n <- length(x)
  f1 <- 0.5 * n 
  f2 <- 0.5 * f1
  com <-  - f1 * log(2)
  slx2 <- 0.5 * sum( log(x) ) 
  sx <- sum(x)
  v <- sx / n
  lik1 <- v * com - n * lgamma(0.5 * v) + 2 * (0.5 * v - 1) * slx2 
  der <-  com - f1 * digamma(0.5 * v) + slx2
  der2 <-  - f2 * trigamma(0.5 * v)
  v <- v - der / der2
  i <- 2
  lik2 <-  v * com - n * lgamma(0.5 * v) + 2 * (0.5 * v - 1) * slx2
  while ( lik2 - lik1 > tol ) {
    i <- i + 1 
	lik1 <- lik2
    der <-  com - f1 * digamma(0.5 * v) + slx2
    der2 <-  - f2 * trigamma(0.5 * v)
    v <- v - der / der2
    lik2 <-  v * com - n * lgamma(0.5 * v) + 2 * (0.5 * v - 1) * slx2
  }
  
  list(iters = i , loglik = lik2 - 0.5 * sx, nu = v)
}


#[export]
exp2.mle <- function(x) {
  a <- min(x)
  n <- length(x)
  b <- sum(x)/n - a
  param <- c(a, b)
  names(param) <- c("a", "b")
  loglik <-  - n * log(b) - 1/n
  list(loglik = loglik, param = param)
}  


#[export]
expmle <- function (x) {
    n <- length(x)
    lambda <- sum(x)/n
    loglik <-  - n * log(lambda) - n
    list(loglik = loglik, lambda = lambda)
}


#[export]
foldnorm.mle <- function(x, tol = 1e-09) {
  n <- length(x)
  m <- sum(x) / n
  x2 <- x^2
  sx2 <- sum(x2)
  a <- sx2 / n - m^2  
  y <- m * x / a
  tanhy <- tanh(y)
  coshy2 <- 1 / cosh(y)^2
  derm <-  - n * m / a + sum(x * tanhy ) / a
  ders <-  - n / 2 / a + (sx2 + n * m^2) / 2 / a^2 - sum( y * tanhy ) / a
  derm2 <-  - n/ a + sum(x2 * coshy2 ) / a^2
  ders2 <- n / 2 / a^2 - (sx2 + n * m^2) / a^3  + 2 * sum( y * tanhy ) / a^2 + sum(y^2 * coshy2 ) / a^2
  derms <- n * m / a^2 - sum(x * tanhy ) / a^2 - sum(y * x * coshy2 ) / a^2
  aold <- c(m, a)
  anew <- aold - c( ders2 * derm - derms * ders, - derms * derm + derm2 * ders) / (derm2 * ders2 - derms^2)
  
  i <- 2
  while ( sum( abs(anew - aold) ) > tol ) {
    i <- i + 1
    m <- anew[1]    ;   a <- anew[2]
    aold <- anew
    y <- m * x / a
    tanhy <- tanh(y)
    coshy2 <- 1 / cosh(y)^2
    derm <-  - n * m / a + sum(x * tanhy ) / a
    ders <-  - n / 2 / a + (sx2 + n * m^2) / 2 / a^2 - sum( y * tanhy ) / a
    derm2 <-  - n/ a + sum(x2 * coshy2 ) / a^2
    ders2 <- n / 2 / a^2 - (sx2 + n * m^2) / a^3  + 2 * sum( y * tanhy ) / a^2 + sum(y^2 * coshy2 ) / a^2
    derms <- n * m / a^2 - sum(x * tanhy ) / a^2 - sum(y * x * coshy2 ) / a^2
    anew <- aold - c( ders2 * derm - derms * ders, - derms * derm + derm2 * ders) / (derm2 * ders2 - derms^2)
  }
  
  m <- anew[1]   ;  a <- anew[2]
  names(anew) <- c("mean", "sigma squared")
  loglik <- n/2 * log( 2 / pi / a) - 0.5 * n * m^2 / a - 0.5 * sum(x2) / a + sum( log( cosh( y ) ) )
  list(iters = i, loglik = loglik, param = anew)
}
 

#[export]
gammamle <- function(x, tol = 1e-09) {
  n <- length(x)
  m <- sum(x)/n
  slx <- sum( log(x) ) / n
  s <- log(m) - slx
  a1 <- 3 - s + sqrt( (s-3)^2 + 24 * s )
  a1 <- a1 / (12 * s)
  a2 <- a1 - ( log(a1) - digamma(a1) - s) / (1/a1 - trigamma(a1) )
  i <- 2
  while ( abs(a2 - a1) > tol) {
    i <- i + 1
    a1 <- a2
    a2 <- a1 - ( log(a1) - digamma(a1) - s) / (1/a1 - trigamma(a1) )
  }
  b <- a2 / m
  loglik <-  - b * n * m + (a2 - 1) * n * slx + n * a2 * log(b) - n * 
        lgamma(a2)
  param <- c(a2, b) 
  names(param) <- c("shape", "rate")
  list(iters = i, loglik = loglik, param = param)
}

#old_gammamle <- function(x, tol = 1e-09) {
#  n <- length(x)
#  sx <- sum(x)
#  sx2 <- sum(x^2)
#  slx <- sum( log(x) )
  # b <-  sx / sx2     ;     a <- b * sx / n         
  # dera <- slx + n * log(b) - n * digamma(a)
  # dera2 <-  - n * trigamma(a)
  # derb2 <-  - n * a / b^2
  # derb <-  - sx  - b * derb2 
  # derab <-  n / b
  # aold <- c(a, b)
  # anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  # i <- 2
  # while ( sum( abs(anew - aold) ) > tol ) {
    # i <- i + 1
    # aold <- anew
    # a <- anew[1]     ;      b <- anew[2] 
    # dera <- slx + n * log(b) - n * digamma(a)
    # dera2 <-  - n * trigamma(a)
    # derb2 <-  - n * a / b^2
    # derb <-  - sx  - b * derb2 
    # derab <-  n / b
    # anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  # }
  # a <- anew[1]    ;     b <- anew[2]
  # loglik <-  - b * sx + (a - 1) * slx + n * a * log(b) - n * lgamma(a)
  # names(anew) <- c("shape", "scale")
  # list(iters = i, loglik = loglik, param = anew) 
# }


#[export]
halfnorm.mle <- function(x) {
  n <- length(x)
  s <- sqrt( sum(x^2) / n )
  loglik <- n / 2 * log(2 / s / pi) - n / 2
  list(loglik = loglik, sigma.squared = s)
}
 
 
#[export]
invgauss.mle <- function(x) {
  n <- length(x)
  sx <- sum(x)
  sx2 <- sum(1/x)
  m <- sx / n
  lambda <- 1 / ( sx2 / n - 1 / m )
  loglik <- n * 0.5 * log(lambda / 2 / pi) - 1.5 * sum( log(x) ) - lambda / 2 / m^2 * (-sx + m^2 * sx2)
  param <- c(m, lambda)
  names(param) <- c("mu", "lambda")
  list(loglik = loglik, param = param)
}  


#[export]
lindley.mle <- function(x) {
  n <- length(x)
  sx <- sum(x)
  a <- sx/n
  b <- a - 1
  delta <- b^2 + 8 * a
  theta <-  0.5 * ( - b + sqrt(delta) ) / a
  loglik <- 2 * n * log(theta) - n * log(1 + theta) + sum( log1p(x) ) - theta * sx
  list(loglik = loglik, theta = theta)
}


#[export]
logcauchy.mle <- function(x, tol = 1e-09) {
  x <- log(x)
  a <- Rfast::cauchy.mle(x, tol = tol)$param
  m <- a[1]  ;  s <- a[2]
  loglik <- length(x) * log(s/pi) - sum(x) - sum( log( (x - m)^2 + s^2) ) 
  list(loglik = loglik, param = a)
}  


#[export]
loglogistic.mle <- function(x, tol = 1e-09) {
  y <- log(x)
  n <- length(x)
  param <- Rfast::logistic.mle(y, tol = tol)$param
  a <- exp(param[1])   ;  b <- 1/param[2]
  loglik <- n * log(b/a) + (b - 1) * sum(y) - n * (b - 1) * log(a) - 2 * sum( log1p( (x/a)^b ) )
  list(loglik = loglik, param = c(a, b) )
}


#[export]
lognorm.mle <- function (x) {
    n <- length(x)
    x <- log(x)
    sx <- sum(x)
    m <- sx/n
    s <- sum(x^2)/n - m^2
    loglik <-  - 0.5 * n * ( log(2 * pi * s) + 1 ) - sx
    param <- c(m, s)
    names(param) <- c("mean", "variance")
    list(loglik = loglik, param = param)
}


#[export]
lomax.mle <- function(x, tol = 1e-09) {
  n <- length(x)
  x2 <- x^2
  m <- sum(x) / n
  s2 <- sum(x2) / n - m^2
  expa <- abs( 2 * s2 / (s2 - m^2) )
  explam <- abs( expa - 1 ) * m
  dera2 <-  - expa * sum( log1p( x / explam) )
  dera <-  n + dera2
  com <- sum( x / (explam + x) )
  derlama <- expa * com
  derlam <-  - n + derlama + com
  derlam2 <-  - ( expa + 1)* explam * sum( x / (explam + x)^2 )
  aold <- log( c(expa, explam) )
  anew <- aold - c( derlam2 * dera - derlama * derlam, - derlama * dera + dera2 * derlam ) / ( dera2 * derlam2 - derlama^2 )
  
  i <- 2
  while ( sum( abs(anew - aold) ) > tol ) {
    i <- i + 1
    aold <- anew
    expa <- exp(anew[1])     ;      explam <- exp(anew[2]) 
    dera2 <-  - expa * sum( log1p( x / explam) )
    dera <-  n + dera2
    com <- sum( x / (explam + x) )
    derlama <- expa * com
    derlam <-  - n + derlama + com
    derlam2 <-  - ( expa + 1) * explam * sum( x / (explam + x)^2 )
    anew <- aold - c( derlam2 * dera - derlama * derlam, - derlama * dera + dera2 * derlam ) / ( dera2 * derlam2 - derlama^2 )
  } 

  a <- exp(anew[1])    ;     lam <- exp(anew[2])
  loglik <- n * log(a / lam) - (a + 1) * sum( log1p(x / lam) )  
  names(anew) <- c("shape", "scale")
  list( iters = i, loglik = loglik, param = c(a, lam) ) 
}


#[export]
maxboltz.mle <- function(x) {
  n <- length(x)
  a <- sqrt( sum(x^2) / 3 / n )
  loglik <- n/2 * log(2 / pi) +  2 * sum( log(x) ) - 1.5 * n - 3 * n * log(a)
  list(loglik = loglik, a = a)
}
  

#[export]
normlog.mle <- function(x) {
  n <- length(x)
  mx <- sum(x) /n
  m <- log( mx )
  sigma <- sum(x^2)/n - mx^2
  loglik <-  - 0.5 * n * log(2 * pi * sigma) - 0.5 * n
  param <- c(mx, m, sigma, sigma * n / (n - 1) )
  names(param) <- c("exp_mu", "mu", "biased variance", "unbiased variance")
  list(loglik = loglik, param = param)
}


#[export]
pareto.mle <- function(x) {
  n <- length(x)
  xm <- min(x) 
  com <- n * log(xm)
  slx <- sum( log(x) )
  a <- n / ( slx - com )
  param <- c(xm, a)
  names(param) <- c("scale", "shape")
  loglik <- n * log(a) + a * com - (a + 1) * slx
  list(loglik = loglik, param = param)
}


#[export]
rayleigh.mle <- function(x) {
  n <- length(x)
  sigma <- 0.5 * sum(x^2) / n
  loglik <- sum( log(x) ) - n * log(sigma) - n
  list(loglik = loglik, sigma = sigma)
}


#[export]
tobit.mle <- function(y, tol = 1e-09) {
  y1 <- y[y >0]  ;  n1 <- length(y1)
  n <- length(y)
  n0 <- n - n1
  sy12 <- sum(y1^2)
  m <- mean(y)   ;  s <- sqrt( sy12/n - m^2 )
  sy1 <- n * m
  com <- dnorm(m, 0, s) / pnorm(-m/s)
  derm <- (sy1 - n1 * m)/s^2 - n0 * com
  derm2 <-  -n1/s^2 - n0 * ( -m /s^2 * com + com^2 ) 
  ders <-  - n1 + (sy12 - 2 * m * sy1 + n1 * m^2)/s^2 + n0 * m * com
  ders2 <-  - 2 * (sy12 - 2 * m * sy1 + n1 * m^2)/s^2 + n0 * m * ( - com + m^2/s^2 * com - com^2 * m ) 
  derms <-  - 2 * (sy1 - n1 * m)/s^2 - n0 * ( - com + m^2/s^2 * com - com^2 * m ) 
  aold <- c(m, log(s))
  anew <- aold - c( ders2 * derm - derms * ders, - derms * derm + derm2 * ders ) / ( derm2 * ders2 - derms^2 )
  i <- 2
  while ( sum( abs(aold - anew) ) > tol ) {
    i <- i + 1
    aold <- anew   
    m <- anew[1]     ;    s <- exp( anew[2] )
    com <- dnorm(m, 0, s) / pnorm(-m/s, 0, 1)
    derm <- (sy1 - n1 * m)/s^2 - n0 * com
    derm2 <-  - n1/s^2 - n0 * ( -m /s^2 * com + com^2 ) 
    ders <-  - n1 + (sy12 - 2 * m * sy1 + n1 * m^2)/s^2 + n0 * m * com
    ders2 <-  - 2 * (sy12 - 2 * m * sy1 + n1 * m^2)/s^2 + n0 * m * ( - com + m^2/s^2 * com - com^2 * m ) 
    derms <-  - 2 * (sy1 - n1 * m)/s^2 + n0 * ( com - m^2/s^2 * com + com^2 * m ) 
    anew <- aold - c( ders2 * derm - derms * ders, - derms * derm + derm2 * ders ) / ( derm2 * ders2 - derms^2 )
  }
  s <- exp(anew[2])
  loglik <-  - 0.5 * n1 * log(2 * pi * s^2) - 0.5 * ( sy12 - 2 * m * sy1 + n1 * m^2 ) / s^2 + n0 * pnorm(-m/s, log.p = TRUE) 
  param <- c(anew[1], s)
  names(param) <- c("location", "scale")
  list(iters = i, loglik = loglik, param = param)
} 


#[export]
weibull.mle <- function(x, tol = 1e-09, maxiters = 100) {
  .Call(Rfast_weibull_mle, x, tol,maxiters)
}

