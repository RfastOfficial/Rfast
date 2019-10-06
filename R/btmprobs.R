#[export]
btmprobs <- function(x, tol = 1e-09) {
  wi <- Rfast::rowsums(x)
  p1 <- p2 <- wi / sum(x)
  D <- length(wi)
  for (i in 1:D)  p2[i] <- wi[i] / ( sum( x[i, ]/ (p1[i] + p1) ) + sum( x[, i]/ (p1[i] + p1) ) )
  p2 <- p2 / sum(p2)    
  j <- 2  
  while ( sum( abs(p2 - p1) ) > tol ) {
    j <- j + 1
    p1 <- p2
    for (i in 1:D)  p2[i] <- wi[i] / ( sum( x[i, ]/ (p1[i] + p1) ) + sum( x[, i]/ (p1[i] + p1) ) )
    p2 <- p2 / sum(p2)    
  }
  list(iters = j, probs = p2)
}


