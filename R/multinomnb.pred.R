#[export]
multinomnb.pred <- function(xnew, m) {
  score <- tcrossprod( xnew, log(m) )
  Rfast::rowMaxs(score)
}