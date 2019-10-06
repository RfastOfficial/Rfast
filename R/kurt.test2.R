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
