#[export]
poissonnb.pred <- function(xnew, m) {
  score <- tcrossprod( xnew, log(m) ) - Rfast::rowsums(m)
  Rfast::rowMaxs(score)
}