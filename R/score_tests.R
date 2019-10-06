#[export]
score.betaregs <- function (y, x, logged = FALSE) {
    param <- Rfast::beta.mle(y)$param
    m1 <- digamma(param[1]) - digamma(param[2])
    z <- log(y) - log(1 - y)
    u <- Rfast::eachcol.apply( x, z - m1 )
    m2 <- trigamma(param[1]) + trigamma(param[2])
    seu <- Rfast::colsums(x^2) * m2
    stat <- u^2/seu
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}


#[export]
score.expregs <- function(y, x, logged = FALSE) {
  lam <- mean(y)
  u <- Rfast::colsums(x * y) * lam - Rfast::colsums(x) 
  vu <- Rfast::colsums(x^2) * lam^4
  stat <- u^2 / vu
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}


#[export]
score.gammaregs <- function (y, x, logged = FALSE) {
    pa <- Rfast::gammamle(y)$param
    m <- pa[1]/pa[2]
    u <- Rfast::colsums(x) - Rfast::eachcol.apply(x, y)/m
    vb <- Rfast::colsums(x^2)/pa[1]
    stat <- u^2/vb
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}


#[export]
score.geomregs <- function(y, x, logged = FALSE) {
  p <- Rfast::geom.mle(y)$prob
  sx <- Rfast::colsums(x)
  u <- (1 - p) * sx - p * Rfast::colsums( y * x)
  vb <- (1 - p ) * Rfast::colsums(x^2)
  stat <- u^2/vb
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}


#[export]
score.glms <- function(y, x, oiko = NULL, logged = FALSE ) {

 if ( is.null(oiko) ) {
   y <- as.numeric(y)
   if ( Rfast::sort_unique.length(y) == 2 ) {
       oiko <- "binomial"
   }   else oiko <- "poisson"
 } 

  n <- length(y) 
  r <- as.numeric( cor(y, x) )
  if ( oiko == "binomial" ) {
    stat <- r * sqrt(n)  
  } else  stat <- ( Var(y, std = T) / sqrt( sum(y) / n ) * sqrt(n - 1) ) * r 

  if ( logged ) {
    pvalue <- log(2) + pt( abs(stat), n - 2, lower.tail = FALSE, log.p = TRUE )
  } else  pvalue <- 2 * pt( abs(stat), n - 2, lower.tail = FALSE )
        
  cbind(stat, pvalue)
}


#[export]
score.invgaussregs <- function (y, x, logged = FALSE) {
    n <- length(y)
    m <- sum(y)/n
    lambda <- 1/( sum(1/y)/n - 1/m )
    u <- Rfast::eachcol.apply(x, m - y ) * lambda
    vu <- m^3 * Rfast::colsums(x^2)
    stat <- u^2/vu
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}


#[export]
score.multinomregs <- function(y, x, logged = FALSE) {
  n <- length(y)
  p <- dim(x)[2]
  dof <- Rfast::sort_unique.length(y) - 1
  if ( dof == 1 ) {
    res <- Rfast::score.glms(y, x, oiko = "binomial", logged = logged) 
  } else {
    m0 <- numeric(dof)
    y1 <- Rfast::design_matrix(y)[,-1]
    m <- Rfast::colmeans(y1)   
    sx <- Rfast::colsums(x)
    sx2 <- Rfast::colsums(x^2)
    vp <- diag(m) - tcrossprod(m)  
    mx <- matrix( rep( m, rep(p, dof) ), ncol = dof )
    ni <- tabulate(y)
    u <- t( rowsum( x, y ) )[, -1] - sx * mx
    stat <- Rfast::mahala(u, m0, vp ) / sx2
    pvalue <- pchisq( stat, dof, lower.tail = FALSE, log.p = logged )
    res <- cbind(stat, pvalue)
  } 
  res
}  


#[export]
score.negbinregs <- function (y, x, logged = FALSE) {
    mod <- Rfast::negbin.mle(y)
    r <- mod$param[2]
    p <- mod$param[1]
    my <- mod$param[3]
    sxy <- Rfast::eachcol.apply(x, y)
    u <- p * sxy - (1 - p) * r * Rfast::colsums(x)
    vu <- Rfast::colsums(x^2) * ( p^2 * (my + my^2/r) )
    stat <- u^2/vu
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}


#[export]
score.weibregs <- function (y, x, logged = FALSE) {
    mod <- Rfast::weibull.mle(y)
    k <- mod$param[1]
    lam <- mod$param[2]
    yk <- y^k
    u <- Rfast::eachcol.apply(x, yk)/lam^k - Rfast::colsums(x)
    vu <- Rfast::colsums(x^2)
    stat <- u^2/vu
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}


#[export]
score.ztpregs <- function (y, x, logged = FALSE) {
    a <- Rfast::ztp.mle(y)
    lam <- a$lambda
    elam <- exp(lam)
    u <- Rfast::eachcol.apply(x, y) - lam * elam/(elam - 1) * Rfast::colsums(x)
    ey <- lam * elam/(elam - 1)
    vu <- Rfast::colsums(x^2) * ( ey * (1 + lam - ey) )
    stat <- u^2/vu
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}

