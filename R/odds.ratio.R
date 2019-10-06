#[export]
odds.ratio <- function(x, a = 0.05, logged = FALSE) {
  or <- x[1] * x[4] / (x[2] * x[3])
  z <- log(or) 
  s <- sqrt( sum(1/x) )
  ci <- c(z - qnorm(1 - a/2) * s, z + qnorm(1 - a/2) * s )
  stat <- abs(z)/s
  if (logged) {
    pvalue <- log(2) + pnorm(stat, lower.tail = FALSE, log.p = TRUE)
  } else  pvalue <- 2 * pnorm(stat, lower.tail = FALSE)
  res <- c(or, pvalue)
  names(res) <- c("odds ratio", "p-value")
  list(res = res, ci = exp(ci))
}