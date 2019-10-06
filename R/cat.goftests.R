
#[export]
cat.goftests <- function(x, props, type = "gsquare", logged = FALSE) {
  obs <- Rfast::colTabulate(x)
  est <- props * dim(x)[1]
  if (type == "chisquare" ) {
    sta <- 2 * obs * log( obs / est)
  } else  sta <- (obs - est)^2 / est
  stat <- Rfast::colsums(sta) 
  pvalue <- pchisq(stat, length(props) - 1, lower.tail = FALSE, log.p = logged )
  cbind(stat, pvalue)
}