#[export]
mahala <- function(x,mu,sigma,ischol = FALSE) {
  if (!is.matrix(x)) 
    x <- matrix(x, 1, length(x))
  if (!is.matrix(sigma)) 
    sigma <- as.matrix(sigma)
  .Call(Rfast_mahaCpp,x,mu,sigma,ischol)
}