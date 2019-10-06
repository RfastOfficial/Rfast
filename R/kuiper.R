#[export]
kuiper <- function(u) {
  u <- Rfast::Sort(u) / (2 * pi)
  n <- length(u)
  i <- (1:n)/n
  f <- sqrt(n)
  Vn <- f * ( max(u - i + 1/n ) + max(i - u) )
  m2 <- (1:50)^2
  a1 <- 4 * m2 * Vn^2
  a2 <- exp( -2 * m2 * Vn^2 )
  b1 <- 2 * ( a1 - 1 ) * a2
  b2 <- 8 * Vn / ( 3 * f ) * m2 * (a1 - 3) * a2
  pvalue <- sum(b1 - b2)
  res <- c(Vn, pvalue)
  names(res) <- c("Test", "p-value")
  res
}





