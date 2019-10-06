#[export]
poisson.anovas <- function(y, ina, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  n <- sum(ni)
  si <- rowsum(y, ina)
  mi <- si/ni
  d1 <- Rfast::colsums( si * log(mi) )
  d0 <- Rfast::colsums(si) * log( Rfast::colsums(si)/n )  
  stat <- 2 * d1 - 2 * d0
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  res <- cbind(stat, pvalue)
  colnames(res) <- c("stat", "p-value")
  res
}
