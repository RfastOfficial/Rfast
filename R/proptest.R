#[export]
proptest <- function(x, n, p, alternative = "unequal", logged = FALSE) {

  s <- sqrt( n * p * (1 - p) )
  stat <- (x - n * p) / s 
  if ( alternative == "unequal" ) {
    if ( !logged ) {
      pvalue <- 2 * pnorm(abs(stat), lower.tail = FALSE)
    } else pvalue <- log(2) + pnorm(abs(stat), lower.tail = FALSE, log.p = TRUE) 
  } else if ( alternative == "greater" ) {
    pvalue <- pnorm(stat, lower.tail = FALSE, log.p = logged) 
  } else if ( alternative == "less" ) {
    pvalue <- pvalue <- pnorm(stat, log.p = logged) 
  }
  
  cbind(stat, pvalue)
}

