#[export]
fish.kent <- function(x, logged = FALSE) {
  n <- dim(x)[1]
  estim <- Rfast::vmf.mle(x)
  k <- estim$kappa
  mu <-  - estim$mu
  mu[1] <- 1 + mu[1]
  i3 <- diag(3)
  P <- i3 - tcrossprod(mu)/mu[1]
  y <- tcrossprod(x, P)[, 2:3]
  lam <- eigen( crossprod(y) )$values/n
  rat <- besselI(k, 0.5, expon.scaled = TRUE)/besselI(k, 2.5, 
      expon.scaled = TRUE)
  Ta <- n * (k/2)^2 * rat * (lam[1] - lam[2])^2
  pvalue <- pchisq(Ta, 2, lower.tail = FALSE, log.p = logged)
  res <- c(Ta, pvalue)
  names(res) <- c("test", "p-value")
  res
}