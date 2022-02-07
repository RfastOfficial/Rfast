#[export]
gammanb <- function(xnew = NULL, x, ina, tol = 1e-07) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  a <- matrix(0, k, d)
  b <- matrix(0, k, d)
  for (i in 1:k)  {
    res <- Rfast::colgammamle( x[ina==i, ] )[, 1:2]
    a[i, ] <- res[, 1]
    b[i, ] <- res[, 2]
  }
  rownames(a) <- rownames(b) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <-  - tcrossprod(b, xnew) + tcrossprod(a - 1, Rfast::Log(xnew)) - 
    Rfast::rowsums( Rfast::Lgamma(a) - a * Rfast::Log(b) )
	est <- Rfast::colMaxs(score)
  }
  list( a = a, b = b, est = est )
}  


#[export]
gammanb.pred <- function(xnew, a, b) {
  score <-  - tcrossprod(b, xnew) + tcrossprod(a - 1, Rfast::Log(xnew)) - 
            Rfast::rowsums( Rfast::Lgamma(a) - a * Rfast::Log(b) )
  Rfast::colMaxs(score)
}


#[export]
gaussian.nb <- function(xnew = NULL, x, ina,parallel = FALSE) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0] 
  k <- length(ni)
  con <- 2 * log( ni )
  m <- rowsum(x, ina) / ni 
  s <- ( rowsum(x^2, ina) - m^2 * ni ) / (ni - 1)
  dets <- Rfast::rowsums( log(s), parallel = parallel )
  if ( !is.null(xnew) ) {
    mat <- mat<-gaussian_nb(xnew,m,s,dets,con,k,parallel)
    est <- Rfast::colMaxs(mat, parallel = parallel)
  }
  rownames(m) <- rownames(s) <- paste("Group", 1:k)
  list(mu = m, sigma = s, ni = ni, est = est )
}


#[export]
gaussiannb.pred <- function(xnew, m, s, ni) {
  con <- 2 * log( ni )
  dets <- Rfast::rowsums( log(s) )
  xnew <- t(xnew)
  k <- dim(m)[1]
  mat <- matrix(nrow = dim(xnew)[2], ncol = k)
  for (j in 1:k)  mat[, j] <-  - Rfast::colsums( (xnew - m[j, ])^2 / s[j, ] ) - dets[j] + con[j]
  Rfast::rowMaxs(mat)
}


#[export]
geom.nb <- function (xnew, x, ina, type = 1) {
    ni <- tabulate(ina)
	ni <- ni[ni > 0]
    if (type == 1) {
        si <- rowsum(x, ina)
        prob <- 1/(1 + si/ni)
        score <- Rfast::rowsums( log(prob) ) + tcrossprod(log(1 - prob), xnew)
    }   else {
        prob <- ni/rowsum(x, ina)
        score <- Rfast::rowsums( log(prob) ) + tcrossprod(log(1 - prob), xnew)
    }
    colMaxs(score)
}


#[export]
geomnb.pred <- function(xnew, prob) {
  score <- Rfast::rowsums(log(prob)) + tcrossprod(log(1 - prob), xnew)
  Rfast::rowMaxs(score)
}  


#[export]
multinom.nb <- function(xnew, x, ina) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  x <- x / Rfast::rowsums(x)  ## normalizes the data, so that each observation sums to 1
  m <- rowsum(x, ina) / ni
  score <- tcrossprod( xnew, log(m) )
  Rfast::rowMaxs(score)
}


#[export]
multinomnb.pred <- function(xnew, m) {
  score <- tcrossprod( xnew, log(m) )
  Rfast::rowMaxs(score)
} 
  

#[export]
poisson.nb <- function (xnew, x, ina) {
    nu <- tabulate(ina)
    m <- rowsum(x, ina)/nu
    score <- tcrossprod(log(m), xnew) - Rfast::rowsums(m)
    colMaxs(score)
}


#[export]
poissonnb.pred <- function(xnew, m) {
  score <- tcrossprod( xnew, log(m) ) - Rfast::rowsums(m)
  Rfast::rowMaxs(score)
}