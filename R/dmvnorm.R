#[export]
dmvnorm <- function(x, mu, sigma, logged = FALSE) {
  quat <-  - 0.5 * Rfast::mahala(x, mu, sigma)
  pow <- length(mu)/2
  con <- (2 * pi)^pow * sqrt( det(sigma) )
  if ( logged ) {
    den <- quat - log( con )
  } else  den <- exp(quat) / con
  den
}