#[export]
ginis <- function(x) {
  n <- dim(x)[1]
  x <- Rfast::colSort(x)
  g <- Rfast::eachcol.apply(x, 1:n)
  2 * g / n / Rfast::colsums(x) - (n + 1)/n
}


#[export]
gini <-function(x) {
  n <- length(x)
  x <- Rfast::Sort(x)
  g <- sum( x * (1:n) )
  2 * g/n / sum(x) - (n + 1)/n
}
