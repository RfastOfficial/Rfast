#[export]
dcor <- function(x, y) {
  .Call( Rfast_dcor, t(x), t(y) )
}


#[export]
dcov <- function(x, y) {
  .Call( Rfast_dcov, t(x), t(y) )
}


#[export]
dcor.ttest <- function(x, y, logged = FALSE) {
  n <- dim(x)[1]
  bcr <- bcdcor(x ,y)
  M <- 0.5 * n * (n - 3)
  dof <- M - 1
  stat <- sqrt(M - 1) * bcr / sqrt(1 - bcr^2)
  pvalue <- pt(stat, dof, lower.tail = FALSE, log.p = logged)
  res <- c(bcr, dof, stat, pvalue)
  names(res) <- c("BC dcor", "df", "statistic", "p-value")
  res 
}


#[export]
dvar <- function(x) {
  if ( is.matrix(x) ) {
    a <- .Call( Rfast_dvar, t(x) )
  } else {
    n <- length(x)
    i <- 1:n
    x <- sort(x)
    sxi <- cumsum(x)
    sxn <- sxi[n]
    ai <- (2 * i - n) * x + sxn - 2 * sxi
    D <- Rfast::Dist(x, square = TRUE, result = "sum")
    a <- 2 * D/n^2 - 2/n^3 * sum(ai^2) + sum(ai)^2/n^4 
    a <- sqrt(a)
  }
  a  
}


#[export]
bcdcor <- function(x,y) {
  .Call( Rfast_bcdcor, t(x), t(y) )
}
