#[export]
sftests <- function(x, logged = FALSE) {
  x <- Rfast::colSort(x)
  n <- dim(x)[1]
  y <- qnorm( ( 1:n - 0.375 ) / (n + 0.25) ) 
  w <- as.vector( cor(y, x)^2 )
  ln <- log(n)
  lln <- log(ln)
  m <-  - 1.2725 + 1.0521 * ( lln - ln )
  s <-  - 0.26758 * ( lln + 2/ln ) + 1.038
  stat <- ( log(1 - w) - m ) / s
  pval <- pnorm(stat, lower.tail = FALSE, log.p = logged)
  res <- cbind(w, stat, pval)
  colnames(res) <- c("squared correlation", "statistic", "p-value")
  res

}   