#[export]
which.is <- function(x,method = "factor") {
  .Call(Rfast_which_is,x,method)
}