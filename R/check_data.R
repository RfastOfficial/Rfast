#[export]
check_data <- function (x, ina = NULL) {
    if (!is.matrix(x)) 
        x <- Rfast::data.frame.to_matrix(x)
    if (is.null(ina)) {
        a <- Rfast::colrange(x)
        b <- which(a == 0)
    } else {
        ina <- as.numeric( as.factor(ina) )
        k <- max(ina)
        a <- matrix(nrow = k, ncol = dim(x)[2])
        for (i in 1:k) a[i, ] <- colrange(x[ina == i, ])
        b <- which(Rfast::colsums(a == 0) > 0)
    }
    b
}
