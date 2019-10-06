#[export]
poisson.anova <- function(y, ina, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  n <- sum(ni)
  si <- Rfast::group(y, ina)
  mi <- si/ni
  d1 <- sum( si * log(mi) )
  d0 <- sum(si) * log(sum( si)/n )  
  stat <- 2 * d1 - 2 * d0
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}
