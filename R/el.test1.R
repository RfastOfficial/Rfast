#[export]
el.test1 <- function (x, mu, tol = 1e-07, logged = FALSE) {
    y <- x - mu
    g <- function(lambda, y) sum(log1p(lambda * y))
    low <-  -1/max(y)
    up <-  -1/min(y)
    if (low < up) {
        mod <- optimise(g, c(low, up), y = y, tol = tol, maximum = TRUE)
        lambda <- mod$maximum
        p <- 1/(1 + lambda * y)
        p <- p/sum(p)
        stat <- 2 * mod$objective
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    }
    else {
        p <- NULL
        lambda <- 10000
        stat <- 1000
        pval <- 0
    }
    info <- c(lambda, stat, pval)
    names(info) <- c("lambda", "stat", "p-value")
    list(info = info, p = p)
}
