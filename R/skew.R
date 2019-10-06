#[export]
skew <- function(x, pvalue = FALSE){
  n <- length(x)
  y <- x - sum(x) / n
  b1 <- n * sum( y^3 )
  nm1 <- n - 1
  b11 <- ( sum( y^2) / nm1 ) ^1.5
  skewness <- b1 / ( nm1 * (n - 2) * b11 )
  if (pvalue) {  
    vars <- 6 * n * nm1 / ( (n - 2) * (n + 1) * (n + 3) )
    stat <- skewness^2/vars
	pval <- pchisq(stat, 1, lower.tail = FALSE)
	skewness <- c(skewness, pval)
	names(skewness) <- c("skewness", "p-value")
  }
  skewness
}
