#[export]
boot.ttest2 <- function(x, y, B = 999) {
   n1 <- length(x)
   n2 <- length(y)
   m1 <- sum(x)/n1
   m2 <- sum(y)/n2
   f1 <- Rfast::Var(x)/n1
   f2 <- Rfast::Var(y)/n2
   tobs <- abs(m1 - m2)/sqrt(f1 + f2)
   mc <- (m1/f1 +  m2/f2) / (1/f1 + 1/f2)
   z1 <- x - m1 + mc
   z2 <- y - m2 + mc
   R <- round( sqrt(B) )

   z1 <- matrix(sample(z1, R * n1, replace = TRUE), ncol = R)
   z2 <- matrix(sample(z2, R * n2, replace = TRUE), ncol = R)
   bm1 <- Rfast::colmeans(z1)
   bm2 <- Rfast::colmeans(z2)
   zx2 <- Rfast::colsums(z1^2)
   zy2 <- Rfast::colsums(z2^2)

   bf1 <- (zx2 - bm1^2 * n1) / ( n1 * (n1 - 1) )
   bf2 <- (zy2 - bm2^2 * n2) / ( n2 * (n2 - 1) )
   fac <- outer(bf1, bf2, "+")
   difa <- outer(bm1, bm2, "-")
   tb <- abs( difa ) / sqrt(fac)
   res <- c( tobs, ( sum(tb > tobs) + 1)/(R^2 + 1)  ) 
   names(res) <- c("stat", "bootstrap p-value")
   res
}
   




