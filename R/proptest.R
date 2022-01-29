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



#[export]
proptests <- function(x1, x2, n1, n2) {
   p1 <- x1 / n1
   p2 <- x2 / n2
   p <- (x1+ x2) / (n1 + n2)
   stat <- (p1 - p2) / sqrt( p * (1 - p) * (1/n1 + 1/n2) )
   pvalue <- 2 * pnorm( abs(stat), lower.tail = FALSE )
   cbind(p1, p2, stat, pvalue)
}   