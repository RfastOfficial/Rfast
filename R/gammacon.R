#[export]
gammacon <- function (y, tol = 1e-08, maxiters = 50) {
  n <- length(y)
  sx <- n
  ly <- log(y)
  sy <- sum(y)
  be <- sum(ly)/n
  m <- exp(-be)
  der2 <- sy * m
  d1 <- sum(ly) + n * log(m) - der2
  der <-  - der2 + sx
  be <- be - der/der2
  m <- exp(-be)
  d2 <- sum(ly) + n * log(m) - der2
  i <- 1
  while ( abs(d2 - d1) > tol & i < maxiters) {
    i <- i + 1
    d1 <- d2
    der2 <- sy * m
    der <-  - der2 + sx
    be <- be - der/der2
    m <- exp(-be)
    d2 <- sum(ly) + n * log(m) - der2
  }
  com <- y * m
  phi <- sum( (com - 1)^2 )/(n - 1)
  devi <-  - 2 * d2 - 2 * n
  list(be = be, deviance = devi, phi = phi)
}
