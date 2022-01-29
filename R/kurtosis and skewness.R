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


#[export]
kurt.test2 <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  vars1 <- 24 * n1 * (n1 -1 )^2 / ( (n1 - 3) * (n1 - 2) * (n1 + 3) * (n1 + 5) )
  vars2 <- 24 * n2 * (n2 - 1)^2 / ( (n2 - 3) * (n2 - 2) * (n2 + 3) * (n2 + 5) )
  stat <- ( Rfast::kurt(x) - Rfast::kurt(y) ) / sqrt( vars1 + vars2  )
  pval <- 2 * pnorm(abs(stat), lower.tail = FALSE) 
  res <- c(stat, pval)
  names(res) <- c("stat", "p-value")  
  res
}


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



#[export]
skew.test2 <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  vars1 <- 6 * n1 * (n1 - 1) / ( (n1 - 2) * (n1 + 1) * (n1 + 3) )
  vars2 <- 6 * n2 * (n2 - 1) / ( (n2 - 2) * (n2 + 1) * (n2 + 3) )
  stat <- ( Rfast::skew(x) - Rfast::skew(y) ) / sqrt( vars1 + vars2  )
  pval <- 2 * pnorm(abs(stat), lower.tail = FALSE) 
  res <- c(stat, pval)
  names(res) <- c("stat", "p-value")  
  res
}