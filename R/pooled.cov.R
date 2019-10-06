
#[export]
pooled.cov <- function(x, ina) {
  s <- crossprod(x)
  ni <- sqrt( tabulate(ina) )
  mi <- rowsum(x, ina)/ni
  k <- length( ni )
  denom <- dim(x)[1] - k 
  for (i in 1:k)  s <- s - tcrossprod(mi[i, ])
  s/denom
}

 


