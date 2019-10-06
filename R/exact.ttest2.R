#[export]
exact.ttest2 <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  z <- c(x, y)
  m1 <- sum(x)/n1
  m2 <- sum(y)/n2
  fac <- Rfast::Var(x)/n1 + Rfast::Var(y)/n2
  tobs <- abs(m1 - m2)/sqrt(fac)
  st <- n1 * m1 + n2 * m2
  st2 <- sum(z^2)
  perms <- Rfast::comb_n(z, n1)
  sx <- Rfast::colsums(perms)
  x2 <- Rfast::colsums(perms^2)
  v1 <- (x2 - sx^2/n1)/(n1 - 1)
  pm1 <- sx / n1

  pm2 <- (st - sx) 
  y2 <- st2 - x2
  v2 <- (y2 - pm2^2/n2)/(n2 - 1) 
  pm2 <- pm2/n2
  tp <-  abs(pm1 - pm2)/sqrt(v1/n1 + v2/n2)
  res <- c( dim(perms)[2], tobs, sum(tp > tobs) / dim(perms)[2] )
  names(res) <- c("permutations", "stat", "p-value" )
  res
}



  
    

  
  