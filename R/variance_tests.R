#[export]
var2test <- function(x, y, alternative = "unequal", logged = FALSE) {
  s1 <- Rfast::Var(x)
  s2 <- Rfast::Var(y)
  n1 <- length(x)
  n2 <- length(y)
  stat <- s1 / s2
    
  if ( alternative == "unequal" ) {
    if ( logged ) {
      a <- pf( stat, n1 - 1, n2 - 1, log.p = TRUE )
      pvalue <- log(2) + min(a, 1 - a)  
    } else {
	  a <- pf( stat, n1 - 1, n2 - 1 )
      pvalue <- 2 * min(a, 1 - a)  
	}
	  
  } else if ( alternative == "greater" ) {
    pvalue <- pf(stat, n1 - 1, n2 - 1, lower.tail = FALSE, log.p = logged)	  
  } else if ( alternative == "less" ) {
    pvalue <- pf(stat, n1 - 1, n2 - 1, log.p = logged)
  }
  
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res  
}



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



#[export]
vartest <- function(x, sigma, alternative = "unequal", logged = FALSE, conf = NULL) {

  n <- dim(x)[1]
  m <- Rfast::colsums(x)
  x2 <- Rfast::colsums(x^2) 
  s <- x2 - m^2 / n
  stat <- (n - 1) * s / sigma
  
  if ( alternative == "unequal" ) {
  
    if ( logged ) {
      a <- pchisq( stat, n - 1, log.p = TRUE )
      pvalue <- log(2) + Rfast::rowMins( cbind(a, 1 - a), value = TRUE ) 
    } else {
      a <- pchisq( stat, n - 1 )
      pvalue <- 2 * Rfast::rowMins( cbind(a, 1 - a), value = TRUE ) 	
    }
	
  } else if ( alternative == "greater" ) {
    pvalue <- pchisq(stat, n - 1, lower.tail = FALSE, log.p = logged)
	
  } else if ( alternative == "less" ) {
    pvalue <- pchisq(stat, n - 1, log.p = logged) 
  }

  res <- cbind(stat, pvalue)

  if ( !is.null(conf) ) {
    a <- 1 - conf
    a1 <- qchisq(1 - a/2, n - 1)
    a2 <- qchisq(a/2, n - 1)
    mat <- cbind(s / a1, s /a2)
    colnames(mat) <- c( paste(a/2, "%", sep = ""), paste(1-a/2, "%", sep = "") )
    res <- cbind(res, mat)
  }
  
  res
}



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