#[export]
quasipoisson.anova <- function (y, ina, logged = FALSE) {
    ni <- tabulate(ina)
	ni <- ni[ni > 0]
    k <- length(ni)
    n <- sum(ni)
    si <- Rfast::group(y, ina)
    sy <- sum(si)
    mi <- si/ni
    d1 <- sum(si * log(mi))
    d0 <- sy * log(sy/n)
    up <- ( 2 * d1 - 2 * d0 ) / (k - 1)
    yi2 <- Rfast::group(y^2, ina)/mi
    phi <- ( sum(yi2) - sy ) / (n - k)   
    stat <- up /  phi
    pvalue <- pf(stat, k - 1, n - k, lower.tail = FALSE, log.p = logged)
    res <- c(stat, pvalue, phi)
    names(res) <- c("stat", "p-value", "phi")
    res
}

