#[export]
anova1 <- function (x, ina, logged = FALSE) {
  ina <- as.numeric(ina)
  k <- max(ina)
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  n <- length(x)
  sx2 <- sum(x^2)
  m <- Rfast::group(x, ina)
  a <- sum(m^2/ni)
  b <- sum(m)^2/n
  mst <- (a - b) / (k - 1)
  mse <- (sx2 - a) / (n - k)
  fa <- mst / mse
  pvalue <- pf(fa, k - 1, n - k, lower.tail = FALSE, log.p = logged)
  tab <- c(fa, pvalue)
  names(tab) <- c("F stat", "p-value")
  tab
}
