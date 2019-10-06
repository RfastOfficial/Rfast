#[export]
colexpmle <- function (x) {
    n <- dim(x)[1]
    lambda <- Rfast::colmeans(x)
    loglik <-  - n * log(lambda) - n
    res <- cbind(lambda, loglik)
    colnames(res) <- c("lambda", "log-likelihood")
    res
}
