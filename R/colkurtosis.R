#[export]
colkurtosis <- function(x, pvalue = FALSE) {
  m <- Rfast::colmeans(x)
  y <- Rfast::eachrow(x, m, oper = "-" )^2
  up <- Rfast::colmeans(y^2)
  down <- Rfast::colmeans(y)^2
  kurtosis <- up / down
  if (pvalue) {
    n <- dim(x)[1]
    vars <- 24 * n * (n -1 )^2 / ( (n - 3) * (n - 2) * (n + 3) * (n + 5) )
    stat <- kurtosis^2 / vars
    pval <- pchisq(stat, 1, lower.tail = FALSE)
	kurtosis <- cbind(kurtosis, pval)
	colnames(kurtosis) <- c("kurtosis", "p-value")
  }  
  kurtosis
}