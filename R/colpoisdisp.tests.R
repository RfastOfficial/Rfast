#[export]
colpoisdisp.tests <- function(y, alternative = "either", logged = FALSE) {
  n <- dim(y)[1]
  m <- Rfast::colmeans(y)
  up <- Rfast::colsums(y^2) - n * m^2 - n * m
  stat <- up / ( sqrt(2 * n) * m )
  if (alternative == "either") {
    if (logged) {
      pval <- log(2) + pnorm(abs(stat), lower.tail = FALSE, log.p = logged)
	} else   pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  } else if (alternative == "over") {
    pval <- pnorm(stat, lower.tail = FALSE, log.p = logged)
  } else pval <- pnorm(stat)
  res <- cbind(stat, pval)
  colnames(res) <- c("stat", "p-value")
  res
}
   
