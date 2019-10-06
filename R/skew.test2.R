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