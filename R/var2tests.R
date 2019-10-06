#[export]
var2tests <- function(x, y = NULL, ina, alternative = "unequal", logged = FALSE) {
  
  if ( is.null(y) )  {
    s1 <- Rfast::colVars( x[ ina == 1, ] )
    s2 <- Rfast::colVars( x[ ina == 2, ] )
    n1 <- sum( ina == 1 )
    n2 <- length(ina) - n1
	
  } else {
    s1 <- Rfast::colVars( x )
    s2 <- Rfast::colVars( y )
    n1 <- dim(x)[1]
    n2 <- dim(y)[1]    
  }
  
  stat <- s1 / s2
    
  if ( alternative == "unequal" ) {
    if ( logged ) {
      a <- pf( stat, n1 - 1, n2 - 1, log.p = TRUE )
      pvalue <- log(2) + Rfast::rowMins( cbind(a, 1 - a), value = TRUE )  
    } else {
	  a <- pf( stat, n1 - 1, n2 - 1 )
      pvalue <- 2 * Rfast::rowMins( cbind(a, 1 - a), value = TRUE )  
	}
	  
  } else if ( alternative == "greater" ) {
    pvalue <- pf(stat, n1 - 1, n2 - 1, lower.tail = FALSE, log.p = logged)	  
  } else if ( alternative == "less" ) {
    pvalue <- pf(stat, n1 - 1, n2 - 1, log.p = logged)
  }

  cbind(stat, pvalue)
}
