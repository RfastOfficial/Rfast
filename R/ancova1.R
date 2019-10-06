#[export]
ancova1 <- function(y, ina, x, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  a <- length(ni) 
  N <- length(ina)
  sy <- sum(y)
  sx <- sum(x)
  com <- sy^2/N
  com2 <- sx * sy / N
  syy <- sum(y^2) - com
  sxx <- sum(x^2) - sx^2/N
  sxy <- sum(x * y) - com2
  tyy <- sum( Rfast::group(y, ina)^2/ni ) - com
  txx <- sum( Rfast::group(x, ina)^2/ni ) - sx^2/N
  txy <- sum( Rfast::group(x, ina)/ni * Rfast::group(y, ina) ) - com2
  eyy <- syy - tyy
  exx <- sxx - txx
  exy <- sxy - txy
  b <- exy / exx
  sse <- eyy - exy^2/exx
  sse2 <- syy - sxy^2/sxx
  dof <- N - a - 1
  mse <- sse / dof
  ftreat <- (sse2 - sse)/(a - 1) / mse
  fb <- sxy^2/sxx/mse
  pvaltreat <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  pvalb <- pf(fb, 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- c(ftreat, fb, pvaltreat, pvalb)
  names(res) <- c("Ftreat", "Fbeta", "pvalue-treat", "pvalue-beta")
  res
} 
  
   