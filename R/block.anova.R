#[export]
block.anova <- function(x, treat, block, logged = FALSE) {
  a <- Rfast::sort_unique.length(treat) 
  b <- Rfast::sort_unique.length(block)
  N <- length(x)
  com <- sum(x)^2/N
  sst <- sum(x^2) - com
  ssa <- sum( Rfast::group(x, treat)^2 ) / b - com
  ssb <- sum( Rfast::group(x, block)^2 ) / a - com
  dof <- (a - 1) * (b - 1) 
  mse <- (sst - ssa - ssb) / dof
  ftreat <- ssa / (a - 1)/mse
  pval <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- c(ftreat, pval)
  names(res) <- c("F-stat", "p-value")
  res
}  
