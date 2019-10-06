#[export]
Lbeta <- function(x,y) {
  .Call(Rfast_Lbeta,x,y)
}