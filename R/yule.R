#[export]
yule <- function(x) {
  or <- x[1] * x[4]/(x[2] * x[3])
  ( sqrt(or) - 1 ) / ( sqrt(or) + 1 )
}
  
