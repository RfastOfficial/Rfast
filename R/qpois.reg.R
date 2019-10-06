#[export]
qpois.reg <- function (x, y, full = FALSE, tol = 1e-09, maxiters = 100) {
    x <- model.matrix(y ~ ., data.frame(x))
    mod <- .Call(Rfast_qpois_reg, x, y, sum(y * log(y), na.rm = TRUE), 
        tol,maxiters)
    res <- list(be = mod$be, devi = mod$deviance, varb = mod$phi * 
        spdinv(mod$L2), phi = mod$phi)
    if (full) {
        be <- mod$be
        varb <- mod$phi * spdinv(mod$L2)
        info <- cbind(be, sqrt(diag(varb)), be^2/diag(varb))
        info <- cbind(info, pchisq(info[, 3], 1, lower.tail = FALSE))
        rownames(info) <- colnames(x)
        colnames(info) <- c("Estimate", "Std. error", "Wald", 
            "p-value")
        res <- list(info = info, devi = mod$deviance, varb = varb, 
            phi = mod$phi)
    }
    res
}
