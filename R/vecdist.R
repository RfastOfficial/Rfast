#[export]
vecdist <- function(x) {
  	.Call(Rfast_vecdist,x)
}