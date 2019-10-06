#[export]
cox.poisrat <- function(x, y, alpha = 0.05) {
  s1 <- sum(x)
  s2 <- sum(y)
  n2 <- length(x)
  n1 <- length(y)
  rat <- s1 / s2 * n2 / n1
  com <- n2 * (2 * s1 + 1) / (n1 *(2 * s2 + 1) )
  low <- com * qf(alpha/2, 2 * s1 + 1, 2 * s2 + 1)
  up <- com * qf(1 - alpha/2, 2 * s1 + 1, 2 * s2 + 1)
  res <- c(rat, low, up)
  names(res) <- c("ratio", paste(alpha/2, "%"), paste(1 - alpha/2, "%") )
  res
}
