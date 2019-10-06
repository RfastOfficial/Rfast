#[export]
geom.anovas <- function (y, ina, type = 1, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  n <- sum(ni)
  if (type == 1) {   
     ell1 <- 0
     si <- rowsum(y, ina)
     pi <-  1/(1 + si/ni )
     ell1 <- Rfast::colsums( ni * log(pi) + si * log(1 - pi) )
     prob <- 1/(1 + Rfast::colsums(si)/n )
     ell0 <- n * log(prob) + Rfast::colsums(si) * log(1 - prob)
     stat <- 2 * ell1 - 2 * ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }  else {
     si <- rowsum(y, ina)
     pi <- ni/si
     ell1 <- Rfast::colsums( ni * log(pi) + (ni/pi - ni) * log(1 - pi) )
     prob <- n/sum(si)
     ell0 <- n * log(prob) + (n/prob - n) * log(1 - prob)
     stat <- 2 * ell1 - 2 * ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }
  res <- cbind(stat, pvalue)
  colnames(res) <- c("stat", "p-value")
  res
}

    
