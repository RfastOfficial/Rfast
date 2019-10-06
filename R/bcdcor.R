#[export]
bcdcor <- function(x,y) {
  .Call(Rfast_bcdcor,t(x),t(y))
}