#[export]
mcnemar <- function(x, y, logged = FALSE) { 
  b1 <- sum(x == 1 & y == 0)
  b2 <- sum(x == 0 & y == 1)
  stat <- ( b1 - b2 )^2 / ( b1 + b2) 
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}