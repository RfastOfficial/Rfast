#[export]
allttests <- function (x, y = NULL, ina, logged = FALSE) {
    if (is.null(y)) {
        x1 <- x[ina == 1, ]
        x2 <- x[ina == 2, ]
        n1 <- sum(ina == 1)
        n2 <- length(ina) - n1
    }
    else {
        x1 <- x
        n1 <- dim(x1)[1]
        x2 <- y
        n2 <- dim(x2)[1]
    }
    m1 <- Rfast::colmeans(x1)
    m2 <- Rfast::colmeans(x2)
    f1 <- Rfast::colVars(x1, suma = n1 * m1)/n1
    f2 <- Rfast::colVars(x2, suma = n2 * m2)/n2
    fac <- outer(f1, f2, "+")
    down <- outer(f1^2/(n1 - 1), f2^2/(n2 - 1), "+")
    dof <- fac^2/down
    difa <- outer(m1, m2, "-")
    stat <- difa/sqrt(fac)
    if (logged) {
        pvalue <- log(2) + pt(abs(stat), dof, lower.tail = FALSE, 
            log.p = TRUE)
    }
    else pvalue <- 2 * pt(abs(stat), dof, lower.tail = FALSE)
    if (is.null(colnames(x))) {
        rownames(stat) <- rownames(pvalue) <- rownames(dof) <- colnames(stat) <- colnames(pvalue) <- colnames(dof) <- paste("Var", 
            1:dim(x)[2])
    }
    else {
        rownames(stat) <- rownames(pvalue) <- rownames(dof) <- colnames(stat) <- colnames(pvalue) <- colnames(dof) <- colnames(x)
    }
    list(stat = stat, pvalue = pvalue, dof = dof)
}
