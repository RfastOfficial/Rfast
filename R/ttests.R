#[export]
ttests <- function(x, y = NULL, ina, paired = FALSE, logged = FALSE, parallel = FALSE) {

  if ( !paired ) {

    if ( is.null(y) ) {
      x1 <- x[ ina == 1, ]
      x2 <- x[ ina == 2, ]
      n1 <- sum( ina == 1 )
      n2 <- length(ina) - n1
    } else {
      x1 <- x     ;    n1 <- dim(x1)[1]
	  x2 <- y     ;    n2 <- dim(x2)[1]
    }

    m1 <- Rfast::colmeans(x1, parallel = parallel)
    m2 <- Rfast::colmeans(x2, parallel = parallel)
    f1 <- Rfast::colVars(x1, suma = n1 * m1, parallel = parallel) / n1
    f2 <- Rfast::colVars(x2, suma = n2 * m2, parallel = parallel) / n2
    fac <- f1 + f2
    dof <- fac^2 / ( f1^2 / (n1 - 1) + f2^2 / (n2 - 1) )
    stat <- ( m1 - m2 ) / sqrt(fac)
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE )
    } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )  
    result <- cbind(stat, pvalue, dof)

  } else {
    n <- dim(x)[1]
    if ( is.null(y) ) {
      z <- x[ ina == 1, ] - x[ ina == 2, ]
    } else  z <- x - y    
    m <- Rfast::colmeans(z, parallel = parallel)
    s <- Rfast::colVars(z, suma = n * m, std = TRUE, parallel = parallel)
    stat <- sqrt(n) * m / s
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE )  
    } else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE )    	
    result <- cbind(stat, pvalue)
  }

  result
}
