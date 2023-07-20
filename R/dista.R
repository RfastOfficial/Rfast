# [export]
dista <- function(xnew, x, type = "euclidean", k = 0, index = FALSE, trans = TRUE, square = FALSE, p = 0, parallel = FALSE) {
    if (type == "canberra1" || type == "canberra2") {
        .Deprecated("The replacement type is \"canberra\"", "Rfast")
    }
    if (type == "hellinger") {
        xnew <- sqrt(xnew)
        x <- sqrt(x)
    }
    x <- .Call(Rfast_dista, t(xnew), t(x), type, square, p, k, index, parallel)
    if (trans) x <- t(x)
    x
}
