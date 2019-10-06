#[export]
kruskaltest <- function (x, ina, logged = FALSE) {
    n <- length(x)
    R <- Rfast::Rank(x)
    ni <- tabulate(ina)
	ni <- ni[ni > 0]
    up <- sum( Rfast::group(R, ina)^2/ni ) - n * (n + 1)^2/4
    down <- sum(R^2) - n * (n + 1)^2/4
    stat <- (n - 1) * up/down
    pvalue <- pchisq(stat, length(ni) - 1, lower.tail = FALSE, 
        log.p = logged)
    res <- c(stat, pvalue)
    names(res) <- c("stat", "p-value")
    res
}
