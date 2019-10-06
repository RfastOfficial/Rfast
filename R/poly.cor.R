#[export]
poly.cor <- function(x, y) {
  ###########
  funa <- function(rho, a1, a2, z) {
    t2 <- ( a2 - rho * z) / sqrt(1 - rho^2)
    t1 <- ( a1 - rho * z) / sqrt(1 - rho^2)
    - sum( log( ( pnorm(t2) - pnorm(t1) ) ) )
  }
  #######
  z <- ( x - mean(x) ) / Rfast::Var(x, std = TRUE) 
  tab <- tabulate(y)
  n <- sum(tab)
  s <- length(tab)
  cuts <- qnorm( cumsum(tab)/n )[-s]
  cts <-  c(-Inf, cuts, Inf)
  a1 <- cts[y]
  a2 <- cts[y + 1]
  rho <- sum(y * z) / Var(y, std = TRUE) / n
  oop <- options(warn = -1)
  on.exit( options(oop) )
  mod <- optim(rho, funa, a1 = a1, a2 = a2, z = z, hessian = TRUE)
  if ( abs(mod$par) > 0.9999 )  mod$par <- sign(mod$par) * 0.9999
  est <- c(mod$par, 1/mod$hessian)
  z <- 0.5 * log( (1 + est[1])/(1 - est[1]) )
  stat <- z^2 * (1 - est[1]^2)^2/est[2]
  pval <- pchisq(stat, 1, lower.tail = FALSE )
  names(est) <- c("correlation", "variance")
  test <- c(stat, pval)
  names(test) <- c("statistic", "p-value")
  list(est = est, test = test)
}
  
