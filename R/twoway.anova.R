#[export]
twoway.anova <- function(y, x1, x2, interact = FALSE, logged = FALSE) {
  a <- Rfast::sort_unique.length(x1)
  b <- Rfast::sort_unique.length(x2)
  N <- length(x1)
  n <- N/a/b
  com <- sum(y)^2/N
  sst <- sum(y^2) - com
  ssa <- sum( Rfast::group(y, x1)^2 ) / (b * n) - com
  ssb <- sum( Rfast::group(y, x2)^2 ) / (a * n) - com
  
  if ( !interact ) {
    dof2 <- N - a - b + 1 
    mse <- (sst - ssa - ssb) / dof2
    fa <- ssa / (a - 1)/mse
    fb <- ssb / (b - 1)/mse
    pvala <- pf(fa, a - 1, dof2, lower.tail = FALSE, log.p = logged) 
    pvalb <- pf(fb, b - 1, dof2, lower.tail = FALSE, log.p = logged) 
    res <- c(fa, fb, pvala, pvalb)
    names(res) <- c("Fstat A", "Fstat B", "p-value A", "p-value B" ) 
  } else {  
    sssub <- sum( aggregate(y ~ x1 + x2, FUN = 'sum')[, -c(1:2)]^2 ) / n - com
    ssab <- sssub - ssa - ssb
    dof2 <- a * b * (n - 1)
    mse <- (sst - sssub) / dof2
    msa <- ssa / (a - 1)
    msb <- ssb / (b - 1)
    msab <- ssab / (a - 1)/(b - 1)
    fa <- msa / mse
    fb <- msb / mse
    fab <- msab / mse
    pvala <- pf(fa, a - 1, dof2, lower.tail = FALSE, log.p = logged) 
    pvalb <- pf(fb, b - 1, dof2, lower.tail = FALSE, log.p = logged) 
    pvalab <- pf(fab, (a - 1) * (b - 1), dof2, lower.tail = FALSE, log.p = logged) 
    res <- c(fa, fb, fab, pvala, pvalb, pvalab)
    names(res) <- c("Fstat A", "Fstat B", "Fstat A:B", "p-value A", "p-value B", "p-value A:B" )
  }

  res
}  

    
    

    
  