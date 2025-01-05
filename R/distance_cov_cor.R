#[export]
dcor <- function(x, y, bc = FALSE) {
  .Call( Rfast_dcor, t(x), t(y), bc )
}


#[export]
dcov <- function(x, y, bc = FALSE) {
  .Call( Rfast_dcov, t(x), t(y), bc )
}

#[export]
bcdcor <- function(x, y) {
  .Deprecated("Please use instead the dcor(..., bc = TRUE)")
}


#[export]
dcor.ttest <- function(x, y, logged = FALSE) {
  n <- dim(x)[1]
  bcr <- dcor(x ,y, bc = TRUE)$dcor
  M <- 0.5 * n * (n - 3)
  dof <- M - 1
  stat <- sqrt(M - 1) * bcr / sqrt(1 - bcr^2)
  pvalue <- pt(stat, dof, lower.tail = FALSE, log.p = logged)
  res <- c(bcr, dof, stat, pvalue)
  names(res) <- c("BC dcor", "df", "statistic", "p-value")
  res 
}


#[export]
dvar <- function (x, bc = FALSE) {
  if ( is.matrix(x) ) {
    a <- .Call(Rfast_dvar, t(x), bc)
  } else {
    n <- length(x)
    i <- 1:n
    x <- sort(x)
    sxi <- cumsum(x)
    sxn <- sxi[n]
    ai <- (2 * i - n) * x + sxn - 2 * sxi
    #D <- Rfast::Dist(x, square = TRUE, result = "sum")
    D <- n * sum(x^2) - sxn^2
    if ( bc ) {
      a <- 2 * D / ( n * (n - 3) ) - 2 / ( n * (n - 2) * (n - 3) ) * sum(ai^2) + 
           sum(ai)^2 / (n * (n - 1) * (n - 2) * (n - 3) )
    } else  a <- 2 * D/n^2 - 2/n^3 * sum(ai^2) + sum(ai)^2/n^4
      a <- sqrt(a)
  }
  a
}



#[export]
pdcor <- function(x, y, z) {
  a1 <- Rfast::dcor(x, y, bc = TRUE)$dcor
  a2 <- Rfast::dcor(x, z, bc = TRUE)$dcor 
  a3 <- Rfast::dcor(y, z, bc = TRUE)$dcor
  up <- a1 - a2 * a3
  down <- sqrt( 1 - a2^2 ) * sqrt( 1 - a3^2 )
  up / down
}
