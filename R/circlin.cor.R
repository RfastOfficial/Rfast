#[export]
circlin.cor <- function(theta, x) {
  n <- length(theta)  ## sample size
  costheta <- cos(theta)
  sintheta <- sin(theta)
  rxc <- as.vector( cor(costheta, x) )  ## and cos(theta) and x correlation
  rxs <- as.vector( cor(sintheta, x) ) ## sin(theta) and x correlation
  rcs <- cor(costheta, sintheta)  ## cos(theta) and sin(theta) correlation
  R2xt <- (rxc^2 + rxs^2 - 2 * rxc * rxs * rcs)/(1 - rcs^2) 
  ## linear-circular correlation
  Ft <- (n - 3) * R2xt / (1 - R2xt) ## F-test statistic value
  pvalue <- pf(Ft, 2, n - 3, lower.tail = FALSE)
  res <- cbind(R2xt, pvalue)
  colnames(res) <- c('R-squared', 'p-value')
  res
}