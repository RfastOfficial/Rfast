#[export]
expmle <- function (x) {
    n <- length(x)
    lambda <- sum(x)/n
    loglik <-  - n * log(lambda) - n
    list(loglik = loglik, lambda = lambda)
}

 
