#[export]
el.test2 <- function(x, y, tol = 1e-07, logged = FALSE) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2

  elpa <- function(mu, x, y, n) {
    g1 <- Rfast::el.test1(x, mu, tol = tol)$info[2]
    g2 <- Rfast::el.test1(y, mu, tol = tol)$info[2]
    g1 + g2  
  }

  m1 <- sum(x)/n1
  m2 <- sum(y)/n2    
  apot <- optimise(elpa, c( min(m1, m2) - 1, max(m1, m2) + 1), x = x, y = y, n = n )
  stat <- apot$objective
  mu <- apot$minimum
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)  
  p1 <- Rfast::el.test1(x, mu, tol = tol)$p
  p2 <- Rfast::el.test1(y, mu, tol = tol)$p
  info <- c(mu, stat, pvalue)
  names(info) <- c("mean", "statistic", "p-value")
  list(info = info, p1 = p1, p2 = p2)
} 

