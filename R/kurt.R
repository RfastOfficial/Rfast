#[export]
kurt <- function(x, pvalue = FALSE){
  n <- length(x)
  y <- ( x - sum(x)/n )^2
  b1 <- sum(y^2)
  b11 <- sum(y)^2 
  kurtosis <- n * b1/b11
  if (pvalue) {
    vars <- 24 * n * (n -1 )^2 / ( (n - 3) * (n - 2) * (n + 3) * (n + 5) )
    stat <- kurtosis^2 / vars
    pval <- pchisq(stat, 1, lower.tail = FALSE)
	kurtosis <- c(kurtosis, pval)
	names(kurtosis) <- c("kurtosis", "p-value")
  }
  kurtosis  
}
