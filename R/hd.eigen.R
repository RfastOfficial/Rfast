#[export]
hd.eigen <- function(x, center = TRUE, scale = FALSE, k = NULL, vectors = FALSE, large = TRUE) {
    n <- dim(x)[1]
    if (center & scale) {
        y <- t(x) - Rfast::colmeans(x)
        y <- y/sqrt(Rfast::rowsums(y^2)) * sqrt(n - 1)
        xx <- Rfast::Crossprod(y, y)
    }   else if (center & !scale) {
        m <- Rfast::colmeans(x) 
        y <- t(x) - m
        xx <- Rfast::Crossprod(y, y)
    }   else if (!center & scale) {
        s <- Rfast::colVars(x, std = TRUE)
        y <- t(x) / s 
        xx <- Rfast::Crossprod(y, y)
    }   else {
        y <- t(x)
        xx <- Rfast::Crossprod(y, y)
    }	
    a <- eigen(xx )
    if ( is.null(k) )   k <- n
    L <- a$values[1:k]
	
    if (vectors) {
	FF <- Rfast::submatrix(a$vectors, 1, n, 1, k)
      vectors <- tcrossprod(y, t(FF) * L^(-0.5) )
    } else  vectors <- NULL
    list(values = L/(n - 1), vectors = vectors)
}

