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
