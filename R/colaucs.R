
#[export]
colaucs <- function (group, preds) {
    ri <- Rfast::colRanks(preds) 
    n <- length(group)
    n1 <- sum(group)
    n0 <- n - n1
    s1 <- rowsum(ri, group)[2, ]
    (s1 - 0.5 * n1 * (n1 + 1))/n0/n1
}
