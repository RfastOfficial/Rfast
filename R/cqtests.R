#[export]
cqtests <- function(x, treat, block, logged = FALSE) {
  k <- length( Rfast::sort_unique.length(treat) )
  cj <- rowsum(x, treat)
  ri <- rowsum(x, block)  
  N <- Rfast::colsums(cj)
  up <- k * (k - 1) * Rfast::colsums(cj^2) - (k - 1) * N^2
  stat <- up / (k * N - colsums(ri^2) )
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}