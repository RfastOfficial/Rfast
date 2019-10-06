#[export]
dmvt <- function(x, mu, sigma, nu, logged = FALSE) {
  p <- length(mu)
  den <- lgamma( (nu + p)/2 ) - lgamma(nu/2) - 0.5 * p * log(pi * nu) - 
  0.5 * log( det(sigma) ) - 0.5 * (nu + p) * log1p( Rfast::mahala(x, mu, sigma)/nu ) 
  if ( logged ) {
    den <- den
  } else  den <- exp(den)
  den
}