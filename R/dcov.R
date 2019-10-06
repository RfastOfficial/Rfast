
#[export]
dcov <- function(x,y) {
  .Call(Rfast_dcov,t(x),t(y))
}