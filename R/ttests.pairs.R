#[export]
ttests.pairs <- function(x, logged = FALSE) {
  
  n <- dim(x)[1]
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, suma = n * m) / n
  fac <- outer(s, s, "+")
  down <- outer( s^2/(n - 1), s^2/(n - 1), "+" )
  dof <- fac^2 / down 
  stat <- outer(m, m, "-") / sqrt(fac)
  if ( logged ) {
    pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE )
  } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )

  if ( is.null( colnames(x) ) ) {
    rownames(stat) <- rownames(pvalue) <- rownames(dof) <- 
    colnames(stat) <- colnames(pvalue) <- colnames(dof) <- paste("Var", 1:dim(x)[2])
  } else {
    rownames(stat) <- rownames(pvalue) <- rownames(dof) <- 
    colnames(stat) <- colnames(pvalue) <- colnames(dof) <- colnames(x)
  }

  list(stat = stat, pvalue = pvalue, dof = dof)  
}