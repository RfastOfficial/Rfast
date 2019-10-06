#[export]
ttest1 <- function(x, m, alternative = "unequal", logged = FALSE, conf = NULL) {

  n <- length(x) 
  xbar <- sum(x)/n
  s <- Rfast::Var(x, std = TRUE)
  stat <- sqrt(n) * ( xbar - m ) / s
  
  if ( alternative == "unequal" ) {
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE ) 
	} else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE) 
  } else if ( alternative == "greater" ) {
    pvalue <- pt( stat, lower.tail = FALSE, n - 1, log.p = logged )	
  } else if ( alternative == "less" ) {
    pvalue <- pt( stat, n - 1, log.p = logged )
  }
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")

  if ( !is.null(conf) ) {  
    a <- 1 - conf
    fac <- qt(1 - a/2, n - 1) * s/sqrt(n)
    mat <- c(xbar - fac, xbar + fac)
    names(mat) <- c( paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = "") )
    res <- list(res = res, ci = mat) 
  }

  res
}
