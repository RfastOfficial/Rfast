#[export]
colskewness <- function(x, pvalue = FALSE) {
  m <- Rfast::colmeans(x)
  y <- Rfast::eachrow(x, m, oper = "-" )
  n <- dim(x)[1]
  nm1 <- n - 1
  up <- n * Rfast::colsums(y^3)
  down <- ( Rfast::colsums(y^2) / nm1 )^1.5
  skewness <- up / ( nm1 * (n - 2) * down )
  if (pvalue) {  
    vars <- 6 * n * nm1 / ( (n - 2) * (n + 1) * (n + 3) )
    stat <- skewness^2/vars
	pval <- pchisq(stat, 1, lower.tail = FALSE)
	skewness <- cbind(skewness, pval)
	colnames(skewness) <- c("skewness", "p-value")
  }
  skewness
}

