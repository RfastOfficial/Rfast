#[export]
poisson.nb <- function (xnew, x, ina) {
    nu <- tabulate(ina)
    m <- rowsum(x, ina)/nu
    score <- tcrossprod(log(m), xnew) - Rfast::rowsums(m)
    colMaxs(score)
}
