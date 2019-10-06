#[export]
watson <- function(u) {
  u <- Rfast::Sort(u) / (2 * pi)
  n <- length(u)
  i <- (1:n)/n
  Wn <- sum( ( ( u - i + 0.5/n ) - ( sum(u) / n - 0.5 ) )^2 ) + 1 / ( 12 * n )
  m <- 1:20
  pvalue <- 2 * sum( ( - 1 )^( m - 1 ) * exp(-2 * m^2 * pi^2 * Wn) )
  res <- c(Wn, pvalue)
  names(res) <- c("Test", "p-value")
  res
}
