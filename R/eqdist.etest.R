#[export]
eqdist.etest <- function(y, x, R = 999) {

  #dii <- total.vecdist(x)
  ni <- length(x)
  i <- 1:ni
  x <- Rfast::Sort(x)
  sx <- sum(x)
  dii <- 2 * sum(i * x ) - (ni + 1) * sx
 
  #djj <- total.vecdist(y)
  nj <- length(y)
  j <- 1:nj
  y <- Rfast::Sort(y)
  sy <- sum(y)
  djj <- 2 * sum(j * y ) - (nj + 1) * sy
  
  s <- sx + sy  ## total sum
  n <- ni + nj  ## total sample size
  dis <- ( nj * dii / ni + ni * djj / nj ) / n
  
  z <- c(x, y)
  pdis <- numeric(R)
  for ( k in 1:R ) {
    id <- Rfast2::Sample.int(n, n)
    z <- z[id]
    # pdii <- total.vecdist( z[1:ni] )
    # pdjj <- total.vecdist( z[-(1:ni)] )  
    xp <- Rfast::Sort(z[i])
    sxp <- sum(xp)
    pdii <- 2 * sum(i * xp ) - (ni + 1) * sxp
    yp <- Rfast::Sort(z[-i])
    syp <- s - sxp
    pdjj <- 2 * sum(j * yp ) - (nj + 1) * syp
    pdis[k] <- ( nj * pdii / ni + ni * pdjj / nj ) / n
  }

  ( sum(pdis >= dis) + 1 ) / (R + 1)
}
  

#total.vecdist <- function(x) {
#  n <- length(x)
#  i <- 1:n
#  x <- Rfast::Sort(x)
#  2 * sum(i * x ) - (n + 1) * sum(x)
#}

