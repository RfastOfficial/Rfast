#[export]
geom.anova <- function (y, ina, type = 1, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  if (type == 1) {   
     si <- Rfast::group(y, ina)
     pi <-  1/(1 + si/ni )
     ell1 <- sum( ni * log(pi) + si * log(1 - pi) )
     prob <- 1/(1 + sum(si)/sum(ni))
     ell0 <- 2 * sum( dgeom(y, prob, log = TRUE) )
     stat <- 2 * ell1 - ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }  else {
     si <- Rfast::group(y, ina)
     pi <- ni/si
     ell1 <- sum( ni * log(pi) + (ni/pi - ni) * log(1 - pi) )
     n <- sum(ni)
     prob <- n/sum(si)
     ell0 <- n * log(prob) + (n/prob - n) * log(1 - prob)
     stat <- 2 * ell1 - 2 * ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}

