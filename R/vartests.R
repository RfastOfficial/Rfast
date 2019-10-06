#[export]
vartests <- function(x, ina, type = "levene", logged = FALSE) {
  ## type can be either "levene" or "bf" (Brown-Forsythe)
  ina <- as.numeric(ina)  
  k <- max(ina)
  ni <- tabulate(ina)
  n <- dim(x)[1]
  z <- t(x)

  if ( type == "levene" )  {
    for (i in 1:k) {
      xina <- x[ina == i, ]
      z[, ina == i ] <- t( xina - Rfast::colmeans(xina) )
    }

  } else {
    for (i in 1:k) {  
      xina <- x[ina == i, ]
      z[, ina == i ] <- t( xina - Rfast::colMedians( xina ) )
    }
  }

  sz2 <- Rfast::rowsums(z^2) 
  m <- matrix(nrow = k, ncol = dim(z)[1]) 
  m <- rowsum(t(z), ina)
  a <- Rfast::colsums( m^2 / ni )
  b <- Rfast::colsums( m )^2 / n 
  mst <- (a - b) / ( k - 1)
  mse <- (sz2 - a) / ( n - k)
  fa <- mst / mse
  pvalue <- pf(fa, k - 1, n - k, lower.tail = TRUE, log.p = logged)
  cbind(fa, pvalue)
}