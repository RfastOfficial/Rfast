#[export]
mvkurtosis <- function(x) {
  n <- dim(x)[1]
  m <- Rfast::colmeans(x)
  s <- (crossprod(x) - n * tcrossprod(m))/(n - 1)
  sum( Rfast::mahala(x, m, s)^2 ) / n 
}