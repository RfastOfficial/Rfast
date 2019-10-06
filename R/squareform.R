#[export]	
squareform <- function(x) {
  .Call(Rfast_squareform_c,x)
}