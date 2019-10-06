#[export]
mcnemars <- function(x, y = NULL, ina, logged = FALSE) { 
  if ( is.null(y) ) {
    difa <- x[ina == 1, ] - x[ina == 2, ]	
  } else  difa <- x - y
  b1 <- Rfast::colsums( difa == 1 )
  b2 <- Rfast::colsums( difa == -1 )
  stat <- ( b1 - b2 )^2 / ( b1 + b2) 
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}