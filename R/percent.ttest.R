#[export]
percent.ttest <- function(x, y, logged = FALSE) {
  n1 <- length(x)
  n2 <- length(y)
  p1 <- sum(x)/n1
  p2 <- sum(y)/n2
  j22 <- p2 * (1 - p2) * n2
  j11 <- j22 + p1 * (1 - p1) * n1
  vb <- j11 / (j22 * j11 - j22^2)
  dof <- n1 + n2 - 2
  phi <- sum( (x - p1)^2/p1/(1 - p1) ) + sum( (y - p2)^2/p2/(1 - p2) )
  phi <- phi / dof
  b <- log( p2/(1 - p2) ) - log( p1/(1 - p1) )
  stat <- b^2 / (vb * phi)
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  res <- c(phi, stat, pval)
  names(res) <- c("phi", "stat", "p-value")
  res
}




