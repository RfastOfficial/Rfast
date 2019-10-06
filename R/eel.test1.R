#[export]
eel.test1 <- function(x, mu, tol = 1e-09, logged = FALSE) {

   funa <- function(x, mu, n) { 
     lam1 <- 0 
     elx <- 1  ## exp(lambda * x)
     xelx <- x   ## x * exp(lambda * x)
     selx <- n  ## sum( exp(lambda * x) )
     sxelx <- sum(xelx)  ## sum( x * exp(lambda * x) )
     f <- sxelx / selx - mu
     der <- ( sum(x * xelx) * selx - sxelx^2 ) / selx^2
     lam2 <- lam1 - f / der
     i <- 2
     while ( abs(lam1 - lam2) > tol ) {
       i <- i + 1
       lam1 <- lam2
       elx <- exp(lam1 * x)  ## exp(lambda * x)
       xelx <- x * elx   ## x * exp(lambda * x)
       selx <- sum(elx)  ## sum( exp(lambda * x) )
       sxelx <- sum(xelx)  ## sum( x * exp(lambda * x) ) 
       f <- sxelx / selx - mu
       der <- ( sum(x * xelx) * selx - sxelx^2 ) / selx^2
       lam2 <- lam1 - f / der 
     }
     list(iters = i, lam2 = lam2, p = elx / selx )
   }
   
  n <- length(x)   
  res <- try( funa(x, mu, n), silent = TRUE)

  if ( class(res) == "try-error" )  {
    p <- iters <- NULL    ;   info <- c(0, 10^5, 0)
  } else {
    p <- res$p    
    stat <-  - 2 * sum( log(n * p) )
    pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    info <- c(res$lam2, stat, pvalue)
    iters <- res$iters
  }

  names(info) <- c("lambda", "statistic", "p-value")
  list(iters = iters, info = info, p = p)
}
