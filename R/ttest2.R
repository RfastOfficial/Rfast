
#[export]
ttest2 <- function(x, y, paired = FALSE, logged = FALSE) {

  if ( !paired ) {
    n1 <- length(x)
    n2 <- length(y)
    m1 <- sum(x)/n1
    m2 <- sum(y)/n2
    f1 <- Rfast::Var(x) / n1
    f2 <- Rfast::Var(y) / n2
    fac <- f1 + f2
    dof <- fac^2 / ( f1^2 / (n1 - 1) + f2^2 / (n2 - 1) )
    stat <- ( m1 - m2 ) / sqrt(fac)
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE )
    } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )  
    result <- c(stat, pvalue, dof)
    names(result) <- c("stat", "p-value", "dof")

  } else {
    n <- length(x)
    z <- x - y    
    m <- sum(z)/n
    s <- Rfast::Var(z, std = TRUE)
    stat <- sqrt(n) * m / s
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE )  
    } else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE )    	
    result <- c(stat, pvalue)
    names(result) <- c("stat", "p-value")
  }
  
  result
}
