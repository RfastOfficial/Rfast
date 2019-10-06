#[export]
rowcvs <- function(x, ln = FALSE, unbiased = FALSE) {
  if (ln) {
    s <- Rfast::rowVars( Rfast::Log(x) )
	cv <- sqrt( exp(s) - 1 )
  } else {
    m <- Rfast::rowsums(x)
    n <- dim(x)[2]
    x2 <- Rfast::rowsums(x^2)
    s <- (x2 - m^2/n)/(n - 1)
    cv <- n * sqrt(s) / m
  }	
  if (unbiased) cv <- (1 + 0.25 / n) * cv
  cv
}
