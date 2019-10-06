#[export]
col.yule <- function(x, y = NULL, ina) {
  if (is.null(y)) {
    difa <- 3 * x[ina == 1, ] - x[ina == 2, ]
  } else  difa <- 3 * x - y
  f <- .Call(Rfast_odds_helper, difa)
  f10 <- f[4, ]
  f01 <- f[2, ]
  f11 <- f[3, ]
  f00 <- f[1, ]
  ro <- f11 * f00/(f10 * f01)
  ( sqrt(ro) - 1 ) / ( sqrt(ro) + 1 )   
}
