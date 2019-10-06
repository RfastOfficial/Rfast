#[export]
proptests <- function(x1, x2, n1, n2) {
   p1 <- x1 / n1
   p2 <- x2 / n2
   p <- (x1+ x2) / (n1 + n2)
   stat <- (p1 - p2) / sqrt( p * (1 - p) * (1/n1 + 1/n2) )
   pvalue <- 2 * pnorm( abs(stat), lower.tail = FALSE )
   cbind(p1, p2, stat, pvalue)
}   