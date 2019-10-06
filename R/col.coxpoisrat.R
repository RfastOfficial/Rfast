#[export]
col.coxpoisrat <- function(x, y, alpha = 0.05) {
  s1 <- Rfast::colsums(x)
  s2 <- Rfast::colsums(y)
  n2 <- dim(x)[1]
  n1 <- dim(y)[1]
  rat <- s1 / s2 * n2 / n1
  com <- n2 * (2 * s1 + 1) / (n1 *(2 * s2 + 1) )
  low <- com * qf(alpha/2, 2 * s1 + 1, 2 * s2 + 1)
  up <- com * qf(1 - alpha/2, 2 * s1 + 1, 2 * s2 + 1)
  res <- cbind(rat, low, up)
  colnames(res) <- c("ratio", paste(alpha/2, "%"), paste(1 - alpha/2, "%") )
  res
}
