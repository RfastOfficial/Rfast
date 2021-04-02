#[export]
colmaxboltz.mle <- function(x) {
    n <- dim(x)[1]
    a <- sqrt( Rfast::colsums(x^2) / (3 * n) )
    loglik <- n/2 * log(2/pi) + 2 * Rfast::colsums( Rfast::Log(x) ) - 1.5 * n - 3 * n * log(a)
    res <- cbind(a, loglik)
    colnames(res) <- c("alpha", "loglikelihood")
    res
}


#[export]
colpoisson.mle <- function(x) {
    n <- dim(x)[1]
    sx <- Rfast::colsums(x)
    loglik <-  -sx + sx * log(sx/n) - Rfast::colsums( Rfast::Lgamma(x + 1) )
    res <- cbind(sx/n, loglik)
    colnames(res) <- c("lambda", "log-likelihood")
    res
}


#[export]
colgammamle <- function (x, tol = 1e-07) {
    n <- dim(x)[1]
    m <- Rfast::colmeans(x)
    slx <- Rfast::colmeans(Rfast::Log(x))
    s <- log(m) - slx
    a1 <- 3 - s + sqrt((s - 3)^2 + 24 * s)
    a1 <- a1/(12 * s)
    a2 <- a1 - (log(a1) - Rfast::Digamma(a1) - s)/(1/a1 - Rfast::Trigamma(a1))
    while (max(abs(a2 - a1)) > tol) {
        a1 <- a2
        a2 <- a1 - (log(a1) - Rfast::Digamma(a1) - s)/(1/a1 - Rfast::Trigamma(a1))
    }
    b <- a2/m
    loglik <- -b * n * m + (a2 - 1) * n * slx + n * a2 * log(b) - 
        n * Lgamma(a2)
    res <- cbind(a2, b, loglik)
    colnames(res) <- c("shape", "scale", "log-likelihood")
    res
}


#[export]
colgeom.mle <- function (x, type = 1) {
    if (type == 1) {
        sx <- Rfast::colsums(x)
        n <- dim(x)   
        prob <- 1/(1 + sx/n )
        loglik <- n * log(prob) + sx * log(1 - prob)
    }
    else {
        n <- dim(x)[1]
        prob <- n/Rfast::colsums(x)
        loglik <- n * log(prob) + (n/prob - n) * log(1 - prob)
    }
    res <- cbind(prob, loglik)
	colnames(res) <- c("prob of success", "log-likelihood")
	res
}


#[export]
colinvgauss.mle <- function(x) {
    n <- dim(x)[1]
    sx <- Rfast::colsums(x)
    sx2 <- Rfast::colsums(1/x)
    m <- sx/n
    lambda <- 1/(sx2/n - 1/m)
    loglik <- n * 0.5 * log( 0.5 * lambda/pi) - 1.5 * colsums( Log(x) ) - 0.5 * lambda/m^2 * (-sx + m^2 * sx2)
    res <- cbind(m, lambda, loglik)
    colnames(res) <- c("mu", "lambda", "log-likelihood")
    res
}


#[export]
collaplace.mle <- function(x) {
    n <- dim(x)[1]
    m <- Rfast::colMedians(x)
    b <- Rfast::colmeans( abs( Rfast::eachrow(x, m, oper = "-") ) )
    loglik <- -n * log(2 * b) - n
    res <- cbind(m, b, loglik)
    colnames(res) <- c("location", "scale", "log-likelihood")
    res
}


#[export]
collindley.mle <- function(x) {
    n <- dim(x)[1]
    sx <- Rfast::colsums(x)
    a <- sx/n
    b <- a - 1
    delta <- b^2 + 8 * a
    theta <- 0.5 * (-b + sqrt(delta))/a
    loglik <- 2 * n * log(theta) - n * log1p(theta) + colsums (log1p(x) ) - theta * sx
    res <- cbind(theta, loglik)
    colnames(res) <- c("theta", "log-likelihood")
    res
}


#[export]
colnormal.mle <- function(x) {
  n <- dim(x)[1]
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x)
  ss <- s * (n - 1) / n
  loglik <-  - 0.5 * n * ( log(2 * pi) + log(ss) ) - 0.5 * n
  res <- cbind(m, ss, s, loglik)
  colnames(res) <- c("mean", "biased variance", "unbiased variance", "log-likelihood")
  res
}


#[export]
colnormlog.mle <- function(x) {
    dm <- dim(x)
    n <- dm[1]
    mx <- Rfast::colmeans(x)
    m <- log(mx)
    sigma <- Rfast::colmeans(x^2) - mx^2
    loglik <-  -0.5 * n * log(2 * pi * sigma) - 0.5 * n
    res <- cbind(mx, m, sigma, sigma * n/(n - 1), loglik)
    colnames(res) <- c("exp_mu", "mu", "biased variance", "unbiased variance", "log-lik")
    res
}


#[export]
colpareto.mle <- function(x) {
    n <- dim(x)[1]
    xm <- Rfast::colMins(x, value = TRUE)
    com <- n * log(xm)
    slx <- Rfast::colsums( Rfast::Log(x) )
    a <- n/(slx - com)
    loglik <- n * log(a) + a * com - (a + 1) * slx
    res <- cbind(xm, a, loglik)
    colnames(res) <- c("scale", "shape", "log-likelihood")
    res
}


#[export]
colrayleigh.mle <- function(x) {
    n <- dim(x)[1]
    sigma <- 0.5 * Rfast::colmeans(x^2)
    loglik <- Rfast::colsums( Rfast::Log(x) ) - n * log(sigma) - n
    res <- cbind(sigma, loglik)
    colnames(res) <- c("sigma", "log-likelihood")
    res
}


#[export]
colvm.mle <- function(x, tol = 1e-07) {
    n <- dim(x)[1]
    C <- Rfast::colmeans( cos(x) )
    S <- Rfast::colmeans( sin(x) )
    ep <- (C < 0)
    mu <- atan(S/C)
    mu[ep] <- mu[ep] + pi
    con <- Rfast::colsums( cos( Rfast::eachrow(x, mu, oper = "-" ) ) ) 
    R <- C^2 + S^2
    k1 <- (1.28 - 0.53 * R) * tan(0.5 * pi * sqrt(R))
    der <- con - n * besselI(k1, 1)/besselI(k1, 0)
    a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - 
        besselI(k1, 1)^2
    der2 <- n * a/besselI(k1, 0)^2
    k2 <- k1 + der/der2
    while (max(abs(k1 - k2)) > tol) {
        k1 <- k2
        der <- con - n * besselI(k1, 1)/besselI(k1, 0)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 
            0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
    }
    k2[which(is.na(k2))] <- 709
    loglik <- k2 * con - n * log(2 * pi) - n * (log(besselI(k2, 
        0, expon.scaled = TRUE)) + k2)
    res <- cbind(mu, k2, loglik)
    colnames(res) <- c("mean", "concentration", "log-likelihood")
    res
}


#[export]
colweibull.mle <- function (x, tol = 1e-09, maxiters = 100, parallel = FALSE) {
    res <- .Call("Rfast_colweibull_mle", PACKAGE = "Rfast", x, 
        tol, maxiters, parallel)
    colnames(res) <- c("shape", "scale", "log-lik")
    res
}


#[export]
colexp2.mle <- function(x) {
  a <- Rfast::colMins(x)
  n <- dim(x)[1]
  b <- Rfast::colmeans(x)/n - a
  loglik <-  - n * log(b) - 1/n
  res <- cbind(a, b, loglik)
  colnames(res) <- c("a", "b", "log-likelihood")
  res
}
  

#[export]
colexpmle <- function (x) {
    n <- dim(x)[1]
    lambda <- Rfast::colmeans(x)
    loglik <-  - n * log(lambda) - n
    res <- cbind(lambda, loglik)
    colnames(res) <- c("lambda", "log-likelihood")
    res
}
