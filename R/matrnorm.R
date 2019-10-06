#[export]
matrnorm <- function(n, p) {
  matrix(  RcppZiggurat::zrnorm(n * p), ncol = p)
}