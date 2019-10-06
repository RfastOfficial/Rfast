#[export]
ftest <- function (x, ina, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  m <- Rfast::group(x, ina)/ni
  s <- Rfast::group(x^2, ina)
  s <- (s - m^2 * ni)/(ni - 1)
  w <- ni/s
  W <- sum(w)
  mesi <- sum(w * m)/W
  hi <- (1 - w/W)^2/(ni - 1)
  H <- sum(hi)
  f <- (k^2 - 1)/3/H
  stat <- sum(w * (m - mesi)^2)/(k - 1)/(1 + 2 * (k - 2)/(k^2 - 1) * H)
  pval <- pf(stat, k - 1, f, lower.tail = FALSE, log.p = logged)
  res <- c(stat, pval)
  names(res) <- c("F stat", "p-value")
  res
}
