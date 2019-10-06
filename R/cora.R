#[export]
cora <- function(x) {
    mat <- t(x) - Rfast::colmeans(x)
    mat <- mat / sqrt( Rfast::rowsums(mat^2) )
    tcrossprod(mat)
}