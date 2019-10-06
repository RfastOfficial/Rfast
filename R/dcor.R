
#[export]
dcor <- function(x,y) {
  .Call(Rfast_dcor,t(x),t(y))
}