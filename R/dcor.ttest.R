#[export]
dcor.ttest <- function(x, y, logged = FALSE) {
  n <- dim(x)[1]
  bcr <- bcdcor(x ,y)
  M <- 0.5 * n * (n - 3)
  dof <- M - 1
  stat <- sqrt(M - 1) * bcr / sqrt(1 - bcr^2)
  pvalue <- pt(stat, dof, lower.tail = FALSE, log.p = logged)
  res <- c(bcr, dof, stat, pvalue)
  names(res) <- c("BC dcor", "df", "statistic", "p-value")
  res 
}
