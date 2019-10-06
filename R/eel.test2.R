#[export]
eel.test2 <- function(x, y, tol = 1e-09, logged = FALSE) {

   funa2 <- function(x, y, n1, n2) { 
     lam1 <- 0 
     ## first sample
     elx <- 1  ## exp(lambda * x)
     xelx <- x    ## x * exp(lambda * x)
     selx <- n1  ## sum( exp(lambda * x) )
     sxelx <- sum(xelx)  ## sum( x * exp(lambda * x) )
     ### second sample
     ely <- 1  ## exp(-lambda * y)
     yely <- y   ## y * exp(lambda * y)
     sely <- n2  ## sum( exp(lambda * y) )
     syely <- sum(yely)  ## sum( y * exp(lambda * y) )
     f <- sxelx / selx - syely / sely
     der <- ( sum(x * xelx) * selx - sxelx^2 ) / selx^2 + ( sum(y * yely) * sely - syely^2 ) / sely^2
     lam2 <- lam1 - f / der
     i <- 2
     while ( abs(lam1 - lam2) > tol ) {
       i <- i + 1
       lam1 <- lam2
	   ## first sample
       elx <- exp(lam1 * x)  ## exp(lambda * x)
       xelx <- x * elx   ## x * exp(lambda * x)
       selx <- sum(elx)  ## sum( exp(lambda * x) )
       sxelx <- sum(xelx)  ## sum( x * exp(lambda * x) )
       ## second sample
       ely <- exp(-lam1 * y)  ## exp(-lambda * y)
       yely <- y * ely   ## y * exp(lambda * y)
       sely <- sum(ely)  ## sum( exp(lambda * y) )
       syely <- sum(yely)  ## sum( y * exp(lambda * y) )
       f <- sxelx / selx - syely / sely
       der <- ( sum(x * xelx) * selx - sxelx^2 ) / selx^2 + ( sum(y * yely) * sely - syely^2 ) / sely^2
       lam2 <- lam1 - f / der
     }
     list(iters = i, lam2 = lam2, p1 = elx / selx, p2 = ely / sely )
   }

  n1 <- length(x)   ;   n2 <- length(y)
  res <- try( funa2(x, y, n1, n2), silent = TRUE)

  if ( class(res) == "try-error" )  {
    p1 <- p2 <- iters <- NULL    ;   info <- c(0, 10^5, 0)
    
  } else {
    p1 <- res$p1    ;   p2 <- res$p2
    stat <-  - 2 * sum( log(n1 * p1) ) - 2 * sum( log(n2 * p2) )
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    info <- c(res$lam2, stat, pvalue)
    iters <- res$iters
  }

  names(info) <- c("lambda", "statistic", "p-value")
  list(iters = iters, info = info, p1 = p1, p2 = p2)
}
