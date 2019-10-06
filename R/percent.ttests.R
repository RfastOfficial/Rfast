#[export]
percent.ttests <- function(x, y, logged = FALSE) {
  n1 <- dim(x)[1]
  n2 <- dim(y)[1]
  p1 <- Rfast::colmeans(x)
  p2 <- Rfast::colmeans(y)
  j22 <- p2 * (1 - p2) * n2
  j11a <- p1 * (1 - p1) 
  j11 <- j22 + j11a * n1
  vb <- j11 / (j22 * j11 - j22^2)
  dof <- n1 + n2 - 2
  phi <- Rfast::rowsums( (t(x) - p1)^2 / j11a ) + Rfast::rowsums( (t(y) - p2)^2 /j22 ) * n2
  phi <- phi / dof
  b <- log( p2/(1 - p2) ) - log( p1/(1 - p1) )
  stat <- b^2 / (vb * phi)
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  res <- cbind(phi, stat, pval)
  colnames(res) <- c("phi", "stat", "p-value")
  res
}

