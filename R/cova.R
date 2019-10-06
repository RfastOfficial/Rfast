#[export]
cova <- function(x, center = FALSE) {
  n <- dim(x)[1]
  if ( !center ) {
    m <- sqrt(n) * Rfast::colmeans(x)
    s <- (crossprod(x) - tcrossprod(m))/(n - 1)
  } else { 
    m <- Rfast::colmeans(x)
    x <- eachrow(x, m, oper = "-")
    s <- crossprod(x)/(n - 1)
  }
  s
}


