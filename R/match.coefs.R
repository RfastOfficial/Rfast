#[export]
match.coefs <- function (x, y = NULL, ina, type = "jacc") {
  if (is.null(y)) {
    difa <- 3 * x[ina == 1, ] - x[ina == 2, ]
  } else   difa <- 3 * x - y
  f <- .Call(Rfast_odds_helper, difa)
  f10 <- f[4, ]
  f01 <- f[2, ]
  f11 <- f[3, ]
  if (type == "jacc") {
    mt <- f11 / (f11 + f10 + f01)
	mt <- matrix(mt)
	colnames(mt) <- "jacc"
  } else if (type == "smc") {
    f00 <- f[1, ]
    mt <- (f11 + f00) / colsums(f)
	mt <- matrix(mt)
	colnames(mt) <- "smc"
  } else {
    f00 <- f[1, ]
    jmc <- f11 / (f11 + f10 + f01)
    smc <- (f11 + f00) / colsums(f)
	mt <- cbind(jmc, smc)
	colnames(mt) <- c("jacc", "smc")
  }
  mt  
}
