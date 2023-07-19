# [export]
dista <- function(xnew, x, type = "euclidean", k = 0, index = FALSE, trans = TRUE, square = FALSE, p = 0, parallel = FALSE) {
    if (type == "canberra1" || type == "canberra2") {
        .Deprecated("The replacement type is \"canberra\"", "Rfast")
    }
    x <- .Call(Rfast_dista_index, t(xnew), t(x), type, square, p, k, index, parallel)
    if (trans) x <- t(x)
    x
}
