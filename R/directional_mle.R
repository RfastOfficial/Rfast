#[export]
acg.mle <- function (x, tol = 1e-07) {
    p <- dim(x)[2]
    n <- dim(x)[1]
    mu <- numeric(p)
    lam1 <- Rfast::cova(x)
    maha <- 1/Rfast::mahala(x, mu, lam1)
    down <- sum(maha)
    up <- crossprod(x * maha, x)
    lam2 <- up/down
    i <- 2
    while (sum(abs(lam2 - lam1)) > tol) {
        i <- i + 1
        lam1 <- lam2
        maha <- 1/Rfast::mahala(x, mu, lam1)
        down <- sum(maha)
        up <- crossprod(x * maha, x)
        lam2 <- up/down
    }
    A <- p * lam2
    if ( is.null(colnames(x)) )  colnames(A) <- rownames(A) <- paste("X", 1:p, sep = "")
    else colnames(A) <- rownames(A) <- colnames(x)
    list(iter = i, cova = A)
}


#[export]
iag.mle <- function (x, tol = 1e-07) {
    n <- dim(x)[1]
    mod <- Rfast::vmf.mle(x)
    ka <- mod$kappa
    x23 <- crossprod(x) - n * diag(3)
    m1 <- mod$mu * sqrt(ka)
    a <- as.vector(x %*% m1)
    a2 <- a^2
    pa <- pnorm(a)
    da <- dnorm(a)
    gm <- pa + a2 * pa + a * da
    der <- 2 * (a * pa + da) * x
    sqa <- sqrt(pa) * x/sqrt(gm)
    sqa2 <- der/gm
    fm1 <- Rfast::colsums(a * x) - n * m1 + Rfast::colsums(sqa2)
    fm2 <- x23 + 2 * crossprod(sqa) - crossprod(sqa2)
    m2 <- m1 - solve(fm2, fm1)
    i <- 2
    while (sum(abs(m2 - m1)) > tol) {
        m1 <- m2
        a <- as.vector(x %*% m1)
        a2 <- a^2
        pa <- pnorm(a)
        da <- dnorm(a)
        gm <- pa + a2 * pa + a * da
        der <- 2 * (a * pa + da) * x
        sqa2 <- der/gm
        fm1 <- Rfast::colsums(a * x) - n * m1 + Rfast::colsums(sqa2)
        sqa <- sqrt(pa) * x/sqrt(gm)
        fm2 <- x23 + 2 * crossprod(sqa) - crossprod(sqa2)
        m2 <- m1 - solve(fm2, fm1)
        i <- i + 1
    }
    rl <- sum(m2^2)
    mesi <- rbind(m2, m2/sqrt(rl))
    rownames(mesi) <- c("Mean vector", "Mean direction")
    if ( is.null(colnames(x)) ) {
        colnames(mesi) <- c("X", "Y", "Z")
    } else  colnames(mesi) <- colnames(x)
    loglik <- -n * log(2 * pi) + 0.5 * sum(a2) - n/2 * rl + sum(log(gm))
    l0 <- -n * log(pi/0.25)
    param <- c(rl, loglik, l0)
    names(param) <- c("Norm of mean", "Log likelihood", "Uniform log-likelihood")
    list(iters = i, mesi = mesi, param = param)
}


#[export]
multivmf.mle <- function(x, ina, tol = 1e-07, ell = FALSE) {
    ni <- tabulate(ina)
    dm <- dim(x)
    p <- dm[2]
    n <- dm[1]
    Apk <- function(p, k) besselI(k, p/2, expon.scaled = TRUE)/besselI(k, 
        p/2 - 1, expon.scaled = TRUE)
    m1 <- rowsum(x, ina)
    Ri <- sqrt( Rfast::rowsums(m1^2) ) / ni
    m <- m1/ni/Ri
    ki <- Ri * (p - Ri^2)/(1 - Ri^2)
    g <- max(ina) 
    for (i in 1:g) {
      k1 <- ki[i]  ;  R <- Ri[i]
      if (k1 < 1e+05) {
        apk <- Apk(p, k1)
        k2 <- k1 - (apk - R)/(1 - apk^2 - (p - 1)/k1 * apk)
        while (abs(k2 - k1) > tol) {
            k1 <- k2
            apk <- Apk(p, k1)
            k2 <- k1 - (apk - R)/(1 - apk^2 - (p - 1)/k1 * apk)
        }
        ki[i] <- k2
      }  else ki[i] <- k1
    }
    loglik <- NULL
    if (ell) {
      loglik <- ni * (p/2 - 1) * log(ki) - 0.5 * ni * p * log(2 * 
        pi) - ni * (log(besselI(ki, p/2 - 1, expon.scaled = TRUE)) + 
        ki) + ki * ni * Ri
    }
    list(loglik = loglik, mi = m, ki = ki) 
}


#[export]
spml.mle <- function(x, tol = 1e-09, maxiters = 100) {
  .Call(Rfast_spml_mle, as.matrix(x), tol, maxiters)
}


#[export]
vm.mle <- function(x, tol = 1e-09) {
 
  n <- length(x)  ## sample size
  C <- sum( cos(x) ) / n 
  S <- sum( sin(x) ) / n
  mu <- ( atan(S/C) + pi * I(C < 0) ) %% (2 * pi)
  con <- sum( cos(x - mu) ) 
  R <- sqrt( C^2 + S^2 ) 
  k1 <- (1.28 - 0.53 * R^2) * tan(0.5 * pi * R)
  if ( k1 < 710 ) { 
    der <- con - n *  besselI(k1, 1, expon.scaled = TRUE) / besselI(k1, 0, expon.scaled = TRUE)
    a <- besselI(k1, 0)^2 / 2  + besselI(k1, 2) * besselI(k1, 0) / 2 - besselI(k1, 1)^2
    der2 <-  n * a / besselI(k1, 0)^2
    k2 <- k1 + der / der2
    while ( abs(k1 - k2) > tol) {
      k1 <- k2 
      der <- con - n *  besselI(k1, 1, expon.scaled = TRUE) / besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2 / 2  + besselI(k1, 2) * besselI(k1, 0) / 2 - besselI(k1, 1)^2
      der2 <-  n * a / besselI(k1, 0)^2
      k2 <- k1 + der/ der2
    }
      
  } else k2 <- k1
  
  param <- c(mu, k2) 
  names(param) <- c("mean", "concentration")
  loglik <- k2 * con - n * log(2 * pi) - n * ( log( besselI(k2, 0, expon.scaled = TRUE) ) + k2 )
  list(loglik = loglik, param = param)
}
  
  
#[export]
vmf.mle <- function (x, tol = 1e-07) {
    dm <- dim(x)
    p <- dm[2]
    n <- dm[1]
    Apk <- function(p, k) besselI(k, p/2, expon.scaled = TRUE)/besselI(k, 
        p/2 - 1, expon.scaled = TRUE)
    m1 <- Rfast::colsums(x)
    R <- sqrt(sum(m1^2))/n
    m <- m1/n/R
    k <- R * (p - R^2)/(1 - R^2)
    if (k < 1e+05) {
      lik1 <- (p/2 - 1) * log(k) - log(besselI(k, p/2 - 1, expon.scaled = TRUE)) - k + k * R
      apk <- Apk(p, k)
      k <- k - (apk - R)/(1 - apk^2 - (p - 1)/k * apk)
      lik2 <- (p/2 - 1) * log(k) - log(besselI(k, p/2 - 1, expon.scaled = TRUE)) - k + k * R
      while ( lik2 - lik1 > tol ) {
        lik1 <- lik2
        apk <- Apk(p, k)
        k <- k - (apk - R)/(1 - apk^2 - (p - 1)/k * apk)
        lik2 <- (p/2 - 1) * log(k) - log(besselI(k, p/2 - 1, expon.scaled = TRUE)) - k + k * R
      }
    }
    else k <- k
    loglik <- n* lik2 - 0.5 * n * p * log(2 * pi) 
    list(loglik = loglik, mu = m, kappa = k)
}


#[export]
wrapcauchy.mle <- function(x, tol = 1e-09) {
  n <- length(x)
  cx <- cos(x)   
  sx <- sin(x)
  cs <- cbind(cx, sx)
  sa <- Rfast::colMedians(cs)
  C <- sa[1]
  S <- sa[2]
  rho <- sqrt(C^2 + S^2)
  a <- ( atan(S/C) + pi * I(C < 0) ) %% (2 * pi) 
  mc <- 2 * rho * cos(a) / (1 + rho^2) 
  ms <- 2 * rho * sin(a) / (1 + rho^2) 
  m1 <- c(mc, ms)
  wi <- 1 / (1 - mc * cx - ms * sx)
  m2 <- Rfast::colsums(wi * cs) / sum(wi)
  i <- 2
  while ( sum( abs(m1 - m2) ) > tol ) {
    i <- i + 1
    m1 <- m2
    wi <- 1 / (1 - m1[1] * cx - m1[2] * sx)
    m2 <- Rfast::colsums(wi * cs) / sum(wi)
  }

  a <- ( atan(m2[2] / m2[1]) + pi * I(m2[1] < 0) ) %% (2 * pi) 
  k <- m2[1] / cos(a)
  rho <- ( 1 - sqrt(1 - k^2) )/abs(k)
  loglik <-  - n * log(2 * pi) + n * log(1 - rho^2) - sum( log1p(rho^2 - 2 * rho * cos(x - a)) )
  param <- c(a, rho)
  names(param) <- c("mean direction", "rho" )
  list(iters = i, loglik = loglik, param = param)
}
