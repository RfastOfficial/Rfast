#[export]
sftest <- function (x, logged = FALSE) {
    x <- Rfast::Sort(x)
    n <- length(x)
    y <- qnorm((1:n - 0.375)/(n + 0.25))
    w <- cor(y, x)^2
    ln <- log(n)
    lln <- log(ln)
    m <- -1.2725 + 1.0521 * (lln - ln)
    s <- -0.26758 * (lln + 2/ln) + 1.038
    stat <- (log(1 - w) - m)/s
    pval <- pnorm(stat, lower.tail = FALSE, log.p = logged)
    res <- c(w, stat, pval)
    names(res) <- c("squared correlation", "statistic", "p-value")
    res
}
