#[export]
colpoisson.anovas <- function (y, x, logged = FALSE) {
  dm <- dim(x)
  n <- dm[1]
  sy <- sum(y)
  d0 <- 2 * sy * log(sy/n)
  stat <- numeric(dm[2])
  pvalue <- numeric(dm[2])
  phi <- numeric(dm[2])
  ######  C++
  for ( i in 1:dm[2] ) {
    ina <- x[, i]
    ni <- tabulate(ina)
    k <- length(ni)
    si <- rowsum(y, ina)
    mi <- si/ni
    d1 <- sum( si * log(mi) )
    stat[i] <- 2 * d1 - d0
    pvalue[i] <- pchisq(stat[i], k - 1, lower.tail = FALSE, log.p = logged)
  }
  ######
  res <- cbind(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}
