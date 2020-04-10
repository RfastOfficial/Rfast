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
