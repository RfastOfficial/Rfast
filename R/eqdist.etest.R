#[export]
eqdist.etest <- function(y, x, R = 999) {
  if ( !is.matrix(x) ) {
    ## univariate
    #dii is the sum of all pairwise distances of x
    ni <- length(x)
    i <- 1:ni
    x <- Rfast::Sort(x)
    sx <- sum(x)
    dii <- 2 * sum(i * x ) - (ni + 1) * sx
    #djj is the sum of all pairwise distances of y
    nj <- length(y)
    j <- 1:nj
    y <- Rfast::Sort(y)
    sy <- sum(y)
    s <- sx + sy
    djj <- 2 * sum(j * y ) - (nj + 1) * sy
    n <- ni + nj  ## total sample size
    dij <- Rfast::dista(x, y, result = "sum")
    dtot <- dii + djj + dij
    stat <- dij - nj * dii / ni - ni * djj / nj

    z <- c(x, y)
    pstat <- numeric(R)
    for ( k in 1:R ) {
      id <- Rfast2::Sample.int(n, ni)
      xp <- Rfast::Sort(z[id])
      sxp <- sum(xp)
      pdii <- 2 * sum(i * xp ) - (ni + 1) * sxp
      yp <- Rfast::Sort(z[-id])
      syp <- s - sxp
      pdjj <- 2 * sum(j * yp ) - (nj + 1) * syp
      pdij <- dtot - pdii - pdjj  
      pstat[k] <- pdij - nj * pdii / ni - ni * pdjj / nj
    }
    ## multivariate
  } else {
    nx <- dim(x)[1]  ;  ny <- dim(y)[1]
    n <- nx + ny
    stat <- Rfast::edist(x, y)
    z <- rbind(x, y) 
    for ( k in 1:R ) {
      id <- Rfast2::Sample.int(n, nx)
      pstat[k] <- Rfast::edist(z[id, ], z[-id, ])
    }
  }
  ( sum( pstat >= stat ) + 1 ) / (R + 1)
} 
  

#total.vecdist <- function(x) {
#  n <- length(x)
#  i <- 1:n
#  x <- Rfast::Sort(x)
#  2 * sum(i * x ) - (n + 1) * sum(x)
#}

