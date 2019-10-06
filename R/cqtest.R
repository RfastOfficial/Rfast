#[export]
cqtest <- function(x, treat, block, logged = FALSE) {
  k <- Rfast::sort_unique.length(treat)
  cj <- Rfast::group(x, treat)
  ri <- Rfast::group(x, block)  
  N <- sum(cj)
  up <- k * (k - 1) * sum(cj^2) - (k - 1) * N^2
  stat <- up / (k * N - sum(ri^2) )
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}