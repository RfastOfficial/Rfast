#[export]
colwatsons <- function(u) {
  u <- Rfast::colSort(u) / (2 * pi)
  n <- dim(u)[1]
  i <- (1:n)/n
  Wn <- Rfast::rowsums( ( t(u - i + 0.5/n) - Rfast::colsums(u) / n + 0.5 )^2 ) + 1 / ( 12 * n )
  m <- 1:20
  pvalue <- 2 * Rfast::colsums( ( - 1 )^( m - 1 ) * exp( outer(-2 * m^2 * pi^2, Wn, "*") ) )
  res <- cbind(Wn, pvalue)
  colnames(res) <- c("Test", "p-value")
  res
}

