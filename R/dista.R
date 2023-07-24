# [export]
dista <- function(xnew, x, type = "euclidean", k = 0, index = FALSE, trans = TRUE, square = FALSE, p = 0, parallel = FALSE) {
    if (type == "canberra1" || type == "canberra2") {
        .Deprecated("The replacement type is \"canberra\"", "Rfast")
    }
    if (type == "hellinger") {
        xnew <- sqrt(xnew)
        x <- sqrt(x)
    }
    if(type != "haversine"){
        xnew <- t(xnew)
        x <- t(x)
    }
    x <- .Call(Rfast_dista, xnew, x, type, square, p, k, index, parallel)
    if (trans) x <- t(x)
    x
}

#[export]
total.dista <- function(xnew, x, type = "euclidean", k = 0, square = FALSE, p = 0, parallel = FALSE) {
    if (type == "hellinger") {
        xnew <- sqrt(xnew)
        x <- sqrt(x)
    }
    if(type != "haversine"){
        xnew <- t(xnew)
        x <- t(x)
    }
    .Call(Rfast_total_dista, xnew, x, type, square, p, k, parallel)
}