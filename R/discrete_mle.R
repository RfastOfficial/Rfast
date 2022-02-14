#[export]
betabinom.mle <- function(x, N, tol = 1e-07) {
  x1 <- x    ;    x2 <- N - x
  z1 <- Rfast::sort_unique(x1)
  n1 <- as.vector( table(x1) )
  z2 <- Rfast::sort_unique(x2) 
  n2 <- as.vector( table(x2) )
  n <- length(x1)
  m1 <- sum(x1) / n    ;    m2 <- sum(x1^2) / n
  down <- N * m2 / m1 - N * m1 - N +m1
  a <- (N * m1 - m2) / down
  b <- (N - m1) * (N - m2/m1) / down
  a <- abs(a)    ;    b <- abs(b)
  co <-  - n * digamma(N + a + b) + n * digamma(a + b)
  dera <- sum( digamma(z1 + a) * n1 ) + co - n * digamma(a)
  derb <- sum( digamma(z2 + b) * n2 ) + co - n * digamma(b)
  derab <-  - n * trigamma(N + a + b) + n * trigamma(a + b)
  dera2 <- sum( trigamma(z1 + a) * n1 ) + derab - n * trigamma(a)
  derb2 <- sum( trigamma(z2 + b) * n2 ) + derab - n * trigamma(b) 
  aold <- c(a, b)
  anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  
  i <- 2
  while ( sum( abs(aold - anew) ) > tol ) {
    i <- i + 1
    aold <- anew     
    a <- anew[1]     ;      b <- anew[2] 
    co <-  - n * digamma(N + a + b) + n * digamma(a + b)
    dera <- sum( digamma(z1 + a) * n1 ) + co - n * digamma(a)
    derb <- sum( digamma(z2 + b) * n2 ) + co - n * digamma(b)
    derab <-  - n * trigamma(N + a + b) + n * trigamma(a + b)
    dera2 <- sum( trigamma(z1 + a) * n1 ) + derab - n * trigamma(a)
    derb2 <- sum( trigamma(z2 + b) * n2 ) + derab - n * trigamma(b) 
    anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  }

   a <- anew[1]    ;     b <- anew[2]
   loglik <- n * Rfast::Lgamma(N + 1) - sum( Rfast::Lgamma(z1 + 1) * n1 ) - sum( Rfast::Lgamma(z2 + 1) * n2 ) +
          sum( Rfast::Lgamma(z1 + a) * n1 ) + sum( Rfast::Lgamma(z2 + b) * n2 ) - n * lgamma(N + a + b) - n * lbeta(a, b)
   pa <- c(a, b) 
   names(pa) <- c( "alpha", "beta" )
   list(iters = i, param = c(a, b), loglik = loglik)
}


#[export]
betageom.mle <- function(x, tol = 1e-07) {

  n <- length(x)
  m1 <- sum(x) / n
  m2 <- sum(x^2) / n
  b <- abs( 2 * (m2 - m1^2) / (m2 - m1 - 2 * m1^2) )
  a <- abs( m1 * (b - 1) )
  ya <- a + x
  y <- ya + b + 1

  com <- n * digamma(a + b) - sum( Digamma(y) )
  dera <- sum( Rfast::Digamma(ya) ) + com - n * digamma(a)
  derb <- n * digamma(b + 1) + com - n * digamma(b)
  derab <- n * trigamma(a + b) - sum( Rfast::Trigamma( y ) )
  dera2 <- sum( Rfast::Trigamma(ya) ) + derab - n * trigamma(a)
  derb2 <- n * trigamma(b + 1) + derab - n * trigamma(b)
  
  aold <- c(a, b)
  anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  
  i <- 2
  while ( sum( abs(anew - aold) ) > tol ) {
    i <- i + 1
    aold <- anew
    a <- aold[1]    ;    b <- aold[2]
    ya <- a + x
    y <- ya + b + 1

    com <- n * digamma(a + b) - sum( Digamma(y) )
    dera <- sum( Rfast::Digamma(ya) ) + com - n * digamma(a)
    derb <- n * digamma(b + 1) + com - n * digamma(b)
    derab <- n * trigamma(a + b) - sum( Rfast::Trigamma( y ) )
    dera2 <- sum( .Call("Rfast_Trigamma", PACKAGE = "Rfast", ya) ) + derab - n * trigamma(a)
    derb2 <- n * trigamma(b + 1) + derab - n * trigamma(b)
    anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  }
  
  list(iters = i, anew = anew)

}
  
#[export]
binom.mle <- function(x, N = NULL, tol = 1e-07) {
  
  if ( !is.null(N) ) {
    if ( length(N) == 1 ) {
	  p <- sum(x) / length(x) / N 
    } else  p <- sum(x) / sum(N)
    loglik <- sum( dbinom(x, N, p, log = TRUE) )
    res <- list(loglik = loglik, prob = p)   
  } else {
   
    n <- length(x)
    k <- max(x)
    sx <- sum(x)
    p <- sx / n / k
    x1 <- x - 1 
    y <- k - x1
    derk <- n * digamma(k + 1) - sum( digamma(y) ) + n * log(1 - p)
    derp <- sx / p + (sx - n * k) / (1 - p)
    derk2 <- n * trigamma(k + 1) - sum( trigamma(y) )
    derp2 <-  - sx / p^2 + (sx - n * k) / (1 - p)^2
    derkp <-  - n / (1 - p)
    aold <- c(k, p)
    anew <- aold - c( derp2 * derk - derkp * derp, - derkp * derk + derk2 * derp ) / ( derk2 * derp2 - derkp^2 )
     
    i <- 2
    while ( sum( abs(aold - anew) ) > tol ) {
      i <- i + 1
      aold <- anew
      k <- aold[1]  ;  p <- aold[2]
      y <- k - x1
      derk <- n * digamma(k + 1) - sum( digamma(y) ) + n * log(1 - p)
      derp <- sx / p + (sx - n * k) / (1 - p)
      derk2 <- n * trigamma(k + 1) - sum( trigamma(y) )
      derp2 <-  - sx / p^2 + (sx - n * k) / (1 - p)^2
      derkp <-  - n / (1 - p)
      anew <- aold - c( derp2 * derk - derkp * derp, - derkp * derk + derk2 * derp ) / ( derk2 * derp2 - derkp^2 )    
    }
    k <- round( anew[1] )   ;   p <- anew[2]
    loglik <- sum( dbinom(x, k, p, log = TRUE) )
    names(anew) <- c("Number of trials", "probability")
    res <- list(iters = i, loglik = loglik, param = anew)
  }   
  
  res
}

  
#[export]
borel.mle <- function(x) {
  n <- length(x)
  sx <- sum(x) 
  m <- 1 - n / sx
  loglik <-  - sx + n + sum( (x - 1) * log(m * x) ) - sum( Rfast::Lgamma(x + 1) )
  list(loglik = loglik, m = m)
}

  
#[export]
dirimultinom.mle <- function(x, tol = 1e-07) {
  dm <- dim(x)
  p <- dm[2]  ## dimensionality
  n <- dm[1]  ## sample size
  rs <- sum(x[1, ])
  a1 <- Rfast::colmeans(x) / rs

  x <- t(x)
  y <- x + a1
  sa <- sum(a1)
  lik1 <- n * lgamma( sa ) - sum( Rfast::Lgamma( rs + sa ) ) - n * sum( lgamma( a1 ) ) + sum( Rfast::Lgamma( y ) )
  f <- n * digamma(sa) - sum( Rfast::Digamma(rs + sa) ) - n * digamma(a1) + Rfast::rowsums( Rfast::Digamma(y) )
  f2 <- matrix(n * trigamma(sa) - sum( Rfast::Trigamma(rs + sa) ), p, p)
  diag(f2) <- diag(f2) - n * trigamma(a1) + Rfast::rowsums( Rfast::Trigamma(y) ) 
  a2 <- a1 - solve(f2, f)
  sa <- sum(a2)
  y <- x + a2
  lik2 <- n * lgamma( sa ) - sum( lgamma( rs + sa ) ) - n * sum( lgamma( a2 ) ) + sum( Rfast::Lgamma( y ) )

  i <- 2
  while ( lik2 - lik1 > tol ) {
    i <- i + 1
    lik1 <- lik2
    a1 <- a2
    f <- n * digamma(sa) - sum( Rfast::Digamma(rs + sa) ) - n * digamma(a1) + Rfast::rowsums( Rfast::Digamma(y) )
    f2 <- matrix(n * trigamma(sa) - sum( Rfast::Trigamma(rs + sa) ), p, p)	
    diag(f2) <-  diag(f2) - n * trigamma(a1) + Rfast::rowsums( Rfast::Trigamma(y) ) 
    a2 <- a1 - solve(f2, f)
    sa <- sum(a2)
    y <- x + a2
    lik2 <- n * lgamma( sa ) - sum( Rfast::Lgamma( rs + sa ) ) - n * sum( lgamma( a2 ) ) + sum( Rfast::Lgamma( y ) ) 
  }
  
  list(iters = i, loglik = lik2, param = a2)  
}

  
#[export]
geom.mle <- function (x, type = 1) {
    if (type == 1) {
        sx <- sum(x) 
        n <- length(x)
        prob <- 1/(1 + sx/n)
        loglik <- n * log(prob) + sx * log(1 - prob)
    }
    else {
        n <- length(x)
        prob <- n/sum(x)
        loglik <- n * log(prob) + (n/prob - n) * log(1 - prob)
    }
    list(loglik = loglik, prob = prob)
}


#[export]
logseries.mle <- function(x, tol = 1e-09) {
  
  n <- length(x)
  sx <- sum(x)
  m <- sx / n
  p1 <- 1 / m
  p <- 1 - p1
  a1 <- log( p / p1 )
  loga1 <- log(p1)
  com <- p / loga1
  der <- m * p1 + com
  der2 <-  - m * p * p1 + p * p1 / loga1 + com^2
  a2 <- a1 - der / der2
  
  i <- 2
  while ( abs(a1 - a2) > tol ) {
    i <- i + 1
    a1 <- a2
    ea <- exp(a1)
    p <- ea / (1 + ea)
    p1 <- 1 - p
    loga1 <- log(p1)
    com <- p / loga1
    der <- m * p1 + com
    der2 <-  - m * p * p1 + p * p1 / loga1 + com^2
    a2 <- a1 - der / der2
  }
   
  list(iters = i, prob = p, loglik = sx * log(p) - sum( log(x) ) - n * log( -loga1 ) )
}


#[export]
multinom.mle <- function(x) {
  N <- sum( x[1, ] )
  p <- Rfast::colmeans(x) / N
  n <- dim(x)[1]  
  loglik <- n * lgamma(N + 1) +  sum( p * log(p) )* N * n - sum( Rfast::Lgamma(x + 1) )
  list(loglik = loglik, prob = p)
}


#[export]
ordinal.mle <- function(y, link = "logit") {
  ina <- tabulate(y)
  k <- length(ina)
  ni <- cumsum(ina)/length(y)
  if ( link == "logit" ) {
    a <- log( ni/(1 - ni) )[-k]
  } else  a <- qnorm(ni)[-k]
  loglik <- sum( ina * log( c(ni[1], diff(ni) ) ) )
  list(loglik = loglik, param = a)
}

 
#[export]
negbin.mle <- function(x, type = 1, tol = 1e-09) {

  n <- length(x)
  sx <- sum(x)
  m <- sx / n
  m2 <- sum(x^2) / n
  p <- 1 - m / (m2 - m^2)
  expr1 <- m / p - m 
  mess <- NULL
  if ( expr1 < 0 )  mess <- c("Negative estimate of number of failures. A geometric or a binomial distribution is perhaps more suitable.")
  expr1 <- abs(expr1)  
  
  if ( type == 1 ) {
    r1 <- log(expr1)
    a <- x + expr1
    f <- sum( Rfast::Digamma(a) ) * expr1 - n * Rfast::Digamma(expr1) * expr1 + 
         n * expr1 * r1 - n * expr1 * log(expr1 + m)  
    f2 <- f + sum( Rfast::Trigamma(a) ) * expr1^2 - n * Rfast::Trigamma(expr1) * expr1^2 + n * expr1 -
          n * expr1^2 / (expr1 + m)
    r2 <- r1 - f / f2
    i <- 2
    while ( abs(r1 - r2) > tol  & r2 < 15) {
      i <- i + 1
      expr1 <- exp(r2)
      r1 <- r2
      a <- x + expr1
      f <- sum( Rfast::Digamma(a) ) * expr1 - n * Rfast::Digamma(expr1) * expr1 + 
           n * expr1 * r1 - n * expr1 * log(expr1 + m)  
      f2 <- f + sum( Rfast::Trigamma(a) ) * expr1^2 - n * trigamma(expr1) * expr1^2 + n * expr1 -
            n * expr1^2 / (expr1 + m)
      r2 <- r1 - f / f2
    }

  } else {
    z <- Rfast::sort_unique(x) 
    nz <- Rfast::Table(x)  
    r1 <- log(expr1)
    a <- z + expr1
    f <- sum( Rfast::Digamma(a) * nz) * expr1 - n * Rfast::Digamma(expr1) * expr1 + 
         n * expr1 * r1 - n * expr1 * log(expr1 + m)  
    f2 <- f + sum( Rfast::Trigamma(a)* nz ) * expr1^2 - n * Rfast::Trigamma(expr1) * expr1^2 + n * expr1 -
          n * expr1^2 / (expr1 + m)
    r2 <- r1 - f / f2
    i <- 2
    while ( abs(r1 - r2) > tol  &  r2 < 15) { 
      i <- i + 1
      expr1 <- exp(r2)
      r1 <- r2
      a <- z + expr1
      f <- sum( Rfast::Digamma(a) * nz ) * expr1 - n * Rfast::Digamma(expr1) * expr1 + 
           n * expr1 * r1 - n * expr1 * log(expr1 + m)  
      f2 <- f + sum( Rfast::Trigamma(a) * nz ) * expr1^2 - n * Rfast::Trigamma(expr1) * expr1^2 + n * expr1 -
            n * expr1^2 / (expr1 + m)
      r2 <- r1 - f / f2
    }
  }
  
    expr2 <- exp(r2)
    p <- expr2 / (expr2 + m)
    param <- c( p, expr2, m )
    names(param) <- c("success probability", "number of failures", "mean")
    loglik <- sum( Rfast::Lgamma(x + expr2) ) - sum( Rfast::Lgamma(x + 1) ) - n * lgamma(expr2) +
              sx * log( 1- p) + n * expr2 * log(p)   
  list(mess = mess, iters = i, loglik = loglik, param = param)
}


#[export]
poisson.mle <- function(x) {
   n <- length(x)
   sx <- sum(x)
   lambda <- sx / n
   loglik <-  - sx + sx * log(lambda) - sum( Rfast::Lgamma(x + 1) )
   list(loglik = loglik, lambda = lambda)
}


#[export]
zip.mle <- function(x, tol = 1e-09) {
  no <- sum(x == 0)
  n <- length(x)  
  prop <- no / n
  n1 <- n - no
  x1 <- x[ x > 0 ]  
  sx <- sum(x1) 
  m <- sx / n
  s <- ( sum(x1^2) - m * sx ) / (n - 1)
  l1 <- s / m + m - 1
  fx <- m - m * exp(-l1) - l1 + prop * l1
  der <- m * exp(-l1) - 1 + prop
  l2 <- l1 - fx / der
  i <- 2
  while ( abs(l2 - l1 )> tol ) {
    i <- i + 1
    l1 <- l2
    fx <- m - m * exp(-l1) - l1 + prop * l1
    der <- m * exp(-l1) - 1 + prop
    l2 <- l1 - fx / der 
  }
  
  p <- 1 - m / l2
  loglik <- no * log( p + (1 - p) * exp(-l2) ) + n1 * log(1 - p) + sum( dpois(x1, l2, log = TRUE) )
  param <- c(l2, p)
  names(param) <- c("lambda", "pi")  
  list(iters = i, loglik = loglik, param = param)
}


#[export]
ztp.mle <- function(x, tol = 1e-09) {
  sx <- sum(x)
  n <- length(x)
  lam1 <- sx / n
  explam <- exp(lam1)
  a1 <- sx / lam1
  a2 <- n * explam / (explam - 1)
  f <- a1 - a2
  f2 <-  - a1 / lam1 + a2 / (explam - 1)
  lam2 <- lam1 - f / f2 
  i <- 2
  while ( abs( lam1 - lam2 ) > tol ) {
    i <- i + 1
    lam1 <- lam2
    explam <- exp(lam1)
    a1 <- sx / lam1
    a2 <- n * explam / (explam - 1)
    f <- a1 - a2
    f2 <-  - a1 / lam1 + a2 / (explam - 1)
    lam2 <- lam1 - f / f2 
  }

  loglik <- sx * log(lam2) - n * log( exp(lam2) - 1 ) -  sum( Lgamma(x + 1) )
  list(iters = i, loglik = loglik, lambda = lam2)
}









