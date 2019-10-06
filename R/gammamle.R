#[export]
gammamle <- function(x, tol = 1e-09) {
  n <- length(x)
  m <- sum(x)/n
  slx <- sum( log(x) ) / n
  s <- log(m) - slx
  a1 <- 3 - s + sqrt( (s-3)^2 + 24 * s )
  a1 <- a1 / (12 * s)
  a2 <- a1 - ( log(a1) - digamma(a1) - s) / (1/a1 - trigamma(a1) )
  i <- 2
  while ( abs(a2 - a1) > tol) {
    i <- i + 1
    a1 <- a2
    a2 <- a1 - ( log(a1) - digamma(a1) - s) / (1/a1 - trigamma(a1) )
  }
  b <- a2 / m
  loglik <-  - b * n * m + (a2 - 1) * n * slx + n * a2 * log(b) - n * 
        lgamma(a2)
  param <- c(a2, b) 
  names(param) <- c("shape", "scale")
  list(iters = i, loglik = loglik, param = param)
}






#old_gammamle <- function(x, tol = 1e-09) {
#  n <- length(x)
#  sx <- sum(x)
#  sx2 <- sum(x^2)
#  slx <- sum( log(x) )
  # b <-  sx / sx2     ;     a <- b * sx / n         
  # dera <- slx + n * log(b) - n * digamma(a)
  # dera2 <-  - n * trigamma(a)
  # derb2 <-  - n * a / b^2
  # derb <-  - sx  - b * derb2 
  # derab <-  n / b
  # aold <- c(a, b)
  # anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  # i <- 2
  # while ( sum( abs(anew - aold) ) > tol ) {
    # i <- i + 1
    # aold <- anew
    # a <- anew[1]     ;      b <- anew[2] 
    # dera <- slx + n * log(b) - n * digamma(a)
    # dera2 <-  - n * trigamma(a)
    # derb2 <-  - n * a / b^2
    # derb <-  - sx  - b * derb2 
    # derab <-  n / b
    # anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  # }
  # a <- anew[1]    ;     b <- anew[2]
  # loglik <-  - b * sx + (a - 1) * slx + n * a * log(b) - n * lgamma(a)
  # names(anew) <- c("shape", "scale")
  # list(iters = i, loglik = loglik, param = anew) 
# }





