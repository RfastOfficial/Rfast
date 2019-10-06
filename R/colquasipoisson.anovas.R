
#[export]
colquasipoisson.anovas <- function (y, x, logged = FALSE) {
  dm <- dim(x)
  n <- dm[1]
  sy <- sum(y)
  d0 <- 2 * sy * log(sy/n)
  d1 <- numeric(dm[2])
  yi2 <- y^2
  stat <- numeric(dm[2])
  pvalue <- numeric(dm[2])
  phi <- numeric(dm[2])
  ######  C++
  for ( i in 1:dm[2] ) {
    ina <- x[, i]
    ni <- tabulate(ina)
    k <- length(ni)
    si <- Rfast::group(y, ina)
    mi <- si/ni
    d1 <- sum(si * log(mi))
    up <- 2 * d1 - d0
    yi2 <- Rfast::group(yi2, ina)/mi
    phi[i] <- (sum(yi2) - sy) / (n - k)    
    stat[i] <- up / phi[i]
    pvalue[i] <- pf(stat[i], k - 1, n - k, lower.tail = FALSE, log.p = logged)
  }
  #######
  res <- cbind(stat, pvalue, phi)
  names(res) <- c("stat", "p-value", "phi")
  res
}
