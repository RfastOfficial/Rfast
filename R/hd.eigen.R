#[export]
hd.eigen <- function(x, center = TRUE, scale = FALSE, k = NULL, vectors = FALSE, large = FALSE) {
    n <- dim(x)[1]
    if ( center & scale ) {
      y <- t(x) - Rfast::colmeans(x)
      y <- y / sqrt( Rfast::rowsums(y^2) ) * sqrt(n - 1)
    } else if ( center & !scale ) {
      m <- Rfast::colmeans(x) 
      y <- t(x) - m
    } else if ( !center & scale ) {
      s <- Rfast::colVars(x, std = TRUE)
      y <- t(x) / s 
    } else {
      y <- t(x)
    }
    if ( large )  {
      xx <- Rfast::Crossprod(y, y)
	} else  xx <- crossprod(y, y)	
    a <- eigen(xx )
    if ( is.null(k) )   k <- n
    L <- a$values[1:k]
	
    if (vectors) {
	FF <- Rfast::submatrix(a$vectors, 1, n, 1, k)
      vectors <- tcrossprod(y, t(FF) * L^(-0.5) )
    } else  vectors <- NULL
    list(values = L/(n - 1), vectors = vectors)
}

