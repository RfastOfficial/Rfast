#[export]
Round <- function(x, digit=0, na.rm = FALSE) {
  .Call(Rfast_Round,x,digit,na.rm)
}