#[export]
cora <- function(x, large = FALSE) {
    mat <- t(x) - Rfast::colmeans(x)
    mat <- mat / sqrt( Rfast::rowsums(mat^2) )
    if ( large ) {
	  Rfast::Tcrossprod(mat, mat)
	} else  tcrossprod(mat)
}


#[export]
cova <- function(x, center = FALSE, large = FALSE) {
  n <- dim(x)[1]
  if ( !center ) {
    m <- sqrt(n) * Rfast::colmeans(x) 
	if (large) {
    s <- (Rfast::Crossprod(x, x) - tcrossprod(m))/(n - 1)
	} else  s <- (crossprod(x) - tcrossprod(m))/(n - 1)
  } else { 
    m <- Rfast::colmeans(x)
    x <- eachrow(x, m, oper = "-")
	if (large) {
	  s <- Rfast::Crossprod(x, x)/(n - 1)
	}  else  s <- crossprod(x)/(n - 1)
  }
  s
}
