#[export]
poisdisp.test <- function(y, alternative = "either", logged = FALSE) {
  n <- length(y)
  m <- sum(y)/n
  up <- (n - 1) * Rfast::Var(y) - n  * m
  stat <- up / sqrt(2 * n * m^2)
  if (alternative == "either") {
    if ( logged ) {   
	  pval <- log(2) + pnorm(abs(stat), lower.tail = FALSE, log.p = logged)
    } else   pval <- 2 * pnorm(abs(stat), lower.tail = FALSE)
  } else if (alternative == "over") {
    pval <- pnorm(stat, lower.tail = FALSE, log.p = logged)
  } else pval <- pnorm(stat)
  res <- c(stat, pval)
  names(res) <- c("stat", "p-value")
  res
}