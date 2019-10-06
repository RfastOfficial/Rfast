#[export]
odds <- function(x, y = NULL, ina, logged = FALSE) {
  
  if ( is.null(y) ) {
    difa <- 3 * x[ina == 1, ] - x[ina == 2, ]
  } else {
    difa <- 3*x - y
  }
  f <- .Call(Rfast_odds_helper,difa)
  f10 <- f[4,]
  f01 <- f[2,]
  f11 <- f[3,]
  f00 <- f[1,]
  
  ro <- f11 * f00 / (f10 * f01)
  
  stat <- log(ro)^2 / (1/f11 + 1/f10 + 1/f01 + 1/f00 )   
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
  
}