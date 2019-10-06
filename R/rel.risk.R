#[export]
rel.risk <- function (x, a = 0.05, logged = FALSE) {
    d1 <- 1/(x[1] + x[3] )
    d2 <- 1/(x[2] + x[4] )
    rr <- x[1] * d1 / (x[2] * d2)  
    z <- log(rr)
    s <- sqrt(1/x[1] + 1/x[2] - d1 - d2)
    ci <- c(z - qnorm(1 - a/2) * s, z + qnorm(1 - a/2) * s)
    stat <- abs(z)/s
    if (logged) {
        pvalue <- log(2) + pnorm(stat, lower.tail = FALSE, log.p = TRUE)
    }
    else pvalue <- 2 * pnorm(stat, lower.tail = FALSE)
    res <- c(rr, pvalue)
    names(res) <- c("relative risk", "p-value")
    list( res = res, ci = exp(ci) )
}
