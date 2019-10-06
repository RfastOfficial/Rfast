#[export]
spdinv <- function(A) {
  chol2inv( Rfast::cholesky(A) )
}