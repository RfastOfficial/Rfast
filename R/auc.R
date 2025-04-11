#[export]
auc <- function(group, preds) {
    ri <- Rfast::Rank(preds)
    n <- length(group)
    n1 <- sum(group)
    n0 <- n - n1
    s1 <- Rfast::group(ri, group)[2]
    (s1 - 0.5 * n1 * (n1 + 1))/n0/n1
}
