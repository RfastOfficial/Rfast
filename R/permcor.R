#[export]
permcor <- function (x, y, R = 999) {
  a <- as.vector( .Call(Rfast_perm_cor, x, y, R) )
  names(a) <- c("cor", "p-value")  
  a
}