#[export]
colcvs <- function(x, ln = FALSE, unbiased = FALSE) {
  if (ln) {
    s <- Rfast::colVars( Log(x) )
	cv <- sqrt( expm1(s) )
  } else {
    m <- Rfast::colsums(x)
    n <- dim(x)[1]
    x2 <- Rfast::colsums(x^2)
    s <- (x2 - m^2/n)/(n - 1)
    cv <- n * sqrt(s) / m
  }	
  if (unbiased) cv <- (1 + 0.25 / n) * cv
  cv
}
