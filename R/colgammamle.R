#[export]
colgammamle <- function (x, tol = 1e-07) {
    n <- dim(x)[1]
    m <- Rfast::colmeans(x)
    slx <- Rfast::colmeans(Rfast::Log(x))
    s <- log(m) - slx
    a1 <- 3 - s + sqrt((s - 3)^2 + 24 * s)
    a1 <- a1/(12 * s)
    a2 <- a1 - (log(a1) - Rfast::Digamma(a1) - s)/(1/a1 - Rfast::Trigamma(a1))
    while (max(abs(a2 - a1)) > tol) {
        a1 <- a2
        a2 <- a1 - (log(a1) - Rfast::Digamma(a1) - s)/(1/a1 - Rfast::Trigamma(a1))
    }
    b <- a2/m
    loglik <- -b * n * m + (a2 - 1) * n * slx + n * a2 * log(b) - 
        n * Lgamma(a2)
    res <- cbind(a2, b, loglik)
    colnames(res) <- c("shape", "scale", "log-likelihood")
    res
}
