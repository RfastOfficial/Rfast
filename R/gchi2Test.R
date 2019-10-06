#[export]
gchi2Test <- function(x, y, logged = FALSE) {
  a <- Rfast::Table(x, y, names = FALSE)
  dof <- prod( dim(a) - 1 )
  rs <- Rfast::rowsums(a) 
  cs <- Rfast::colsums(a)
  est <- outer(rs, cs, "*")/sum(a)
  chi <- sum( (a - est)^2 / est ) 
  g <- 2 * sum( a * log(a/est), na.rm = TRUE ) 
  pchi <- pchisq(chi, dof, lower.tail = FALSE, log.p = logged)
  pg <- pchisq(g, dof, lower.tail = FALSE, log.p = logged)
  res <- rbind( c(chi, pchi, dof), c(g, pg, dof) )
  colnames(res) <- c("stat", "p-value", "dof")
  rownames(res) <- c("X2", "G2")
  res
}