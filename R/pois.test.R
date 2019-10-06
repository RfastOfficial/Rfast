#[export]
pois.test <- function(y, logged = FALSE) {
  stat <- ( Rfast::Var(y)/mean(y) - 1) * sqrt( 0.5 * (length(y) - 1) )
  pval <- 2 * pnorm(abs(stat), lower.tail = FALSE, log.p = logged)
  c(stat, pval)
}