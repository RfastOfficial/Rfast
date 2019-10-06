#[export]
block.anovas <- function(x, treat, block, logged = FALSE) {
  a <- Rfast::sort_unique.length(treat)
  b <- Rfast::sort_unique.length(block)
  N <- dim(x)[1]
  com <- Rfast::colsums(x)^2/N
  sst <- Rfast::colsums(x^2) - com
  ssa <- Rfast::colsums( rowsum(x, treat)^2 ) / b - com
  ssb <- Rfast::colsums( rowsum(x, block)^2 ) / a - com
  dof <- (a - 1) * (b - 1) 
  mse <- (sst - ssa - ssb) / dof
  ftreat <- ssa / (a - 1)/mse
  pval <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- cbind(ftreat, pval)
  colnames(res) <- c("F-stat", "p-value")
  res
}  
