#[export]
list.ftests <- function(x, logged = FALSE) {
    k <- length(x)
    p <- dim(x[[ 1 ]])[2]
    ni <- numeric(k)
    m <- matrix(0, k, p)
    s <- matrix(0, k, p)
    for ( i in 1:k) {
      y <- x[[ i ]]
      ni[i] <- dim(y)[1]
      m[i, ] <- Rfast::colmeans(y)
      s[i, ] <- Rfast::colsums(y^2)
    }
    s <- (s - m^2 * ni)/(ni - 1)
    w <- ni/s
    W <- Rfast::colsums(w)
    mesi <- Rfast::colsums(w * m)/W
    hi <- (1 - w/W)^2/(ni - 1)
    H <- Rfast::colsums(hi)
    f <- (k^2 - 1)/3/H
    stat <- Rfast::rowsums((t(w) * (t(m) - mesi)^2)/(k - 1)/(1 + 2 *
        (k - 2)/(k^2 - 1) * H))
    pvalue <- pf(stat, k - 1, f, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}
