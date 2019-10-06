#[export]
mv.eeltest1 <- function(x, mu, tol = 1e-06) {
  ## x is the multivariate data
  ## xa can also be univariate data
  ## mu is the hypothesized mean
  dm <- dim(x)
  d <- dm[2]
  n <- dm[1]

   eel <- function(x, mu, n, d) {
    lam_old <- numeric(d)
    f1 <- numeric(n)
    f2 <- n
    f2a <- x
    f3 <- Rfast::colsums( f2a )
    f4 <- f3
    f <- f4 - mu
    der <-  - tcrossprod( f4 ) + crossprod(f2a, x) / f2
    lam_new <- lam_old - solve(der, f)
    i <- 2
    ## step 3 and beyond
    while ( sum( abs(lam_new - lam_old ) ) > tol )  {
      i <- i + 1
      lam_old <- lam_new
      f1 <- exp( as.vector( x %*% lam_old ) )
      f2 <- sum(f1)
      f2a <- x * f1
      f3 <- Rfast::colsums( f2a )
      f4 <- f3 / f2
      f <- f4 - mu
      der <-  - tcrossprod( f4 ) + crossprod(f2a, x) / f2
      lam_new <- lam_old - solve(der, f)
    }
    p <- f1 / f2
    stat <-  - 2 * sum( log( n * p) )
    pvalue <- pchisq(stat, d, lower.tail = FALSE)
    info <- c(stat, pvalue)
    names(info) <- c("statistic", "p-value")
    list(p = p, lambda = lam_new, iters = i, info = info)
  }
  
  runtime <- proc.time()
  res <- try( eel(x, mu, n, d), silent = TRUE )
  runtime <- proc.time() - runtime
  res$runtime <- runtime

  if ( class(res) == "try-error" ) {
    res$info[1] <- 1e10
    res$info[2] <- 0
    res$p <- NA
  }

  res
}




