################################
#### Median direction
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
#### References: Fisher, N. I. (1985). Spherical medians.
#### Journal of the Royal Statistical Society. Series B, 47(2): 342-348.
#### Fisher, N. I., Lewis, T., & Embleton, B. J. (1987).
#### Statistical analysis of spherical data. Cambridge university press.
#### Cabrera, J., & Watson, G. S. (1990). On a spherical median related distribution.
#### Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
################################
#[export]
mediandir <- function(x) {
  ## x is the directional data
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  u1 <- Rfast::colMedians(x)
  u1 <- u1 / sqrt( sum(u1^2) )
  ww <- 1 /as.vector( sqrt( abs( 1 - ( x %*% u1 )^2 ) ) )
  if ( is.infinite( max(ww) ) ) {
    u2 <- u1
  }  else  u2 <- Rfast::colsums(x * ww )
  u2 <- u2 / sqrt( sum( u2^2 ) )
  while ( sum( abs (u2 - u1 ) ) > 1e-8 ) {
    u1 <- u2
    ww <- 1/as.vector( sqrt( abs( 1 - ( x %*% u1 )^2 ) ) )
    if ( is.infinite( max(ww) ) ) {
       u2 <- u1
    }  else  u2 <- Rfast::colsums(x * ww )
    u2 <- u2 / sqrt( sum( u2^2 ) )
  }

  u2
}
