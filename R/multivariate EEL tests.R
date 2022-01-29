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



#[export]
mv.eeltest2 <- function(y1, y2, tol = 1e-07, R = 0) {
  ## y1 and y2 are the two matrices containing the multivariate (or univariate data)
  ## R is the type of calibration
  ## If R = 0, the chi-square distribution is used
  ## If R = 1, the James corrected chi-square distribution is used
  ## If R = 2, the MNV F distribution is used
  dm <- dim(y1)
  n1 <- dm[1]    ;    d <- dm[2]
  n2 <- dim(y2)[2]

  eel2 <- function(x, y, n1, n2, d) {
    lam1 <- numeric(d)
    fx1 <- numeric(n1)
    fx2 <- n1
    fx2a <- x
    fx3 <- colsums( fx2a )
    fx4 <- fx3 / fx2
    fy1 <- numeric(n2)
    fy2 <- n2
    fy2a <- y
    fy3 <- colsums( fy2a )
    fy4 <- fy3 / fy2
    f <- fx4 - fy4
    der <-  - tcrossprod( fx4 ) + crossprod(fx2a, x) / fx2 - tcrossprod( fy4 ) + crossprod(fy2a, y) / fy2
    lam2 <- lam1 - solve(der, f)
    i <- 2
    while ( sum( abs(lam2 - lam1 ) ) > tol )  {
      i <- i + 1
      lam1 <- lam2
      fx1 <- exp( as.vector( x %*% lam1 ) )
      fx2 <- sum(fx1)
      fx2a <- x * fx1
      fx3 <- colsums( fx2a )
      fx4 <- fx3 / fx2
      fy1 <- exp( as.vector( - y %*% lam1 ) )
      fy2 <- sum(fy1)
      fy2a <- y * fy1
      fy3 <- colsums( fy2a )
      fy4 <- fy3 / fy2
      f <- fx4 - fy4
      der <-  - tcrossprod( fx4 ) + crossprod(fx2a, x) / fx2 - tcrossprod( fy4 ) + crossprod(fy2a, y) / fy2
      lam2 <- lam1 - solve(der, f)
    }

    p1 <- fx1 / fx2
    p2 <- fy1 / fy2
    stat <-  - 2 * sum( log( n1 * p1) ) - 2 * sum( log(n2 * p2) )
    pvalue <- pchisq(stat, d, lower.tail = FALSE)
    info <- c(stat, pvalue, d)
    names(info) <- c("statistic", "p-value", "degrees of freedom")
    list(p1 = p1, p2 = p2, lambda = lam2, iters = i, info = info)
  }

  runtime <- proc.time()
  res <- try( eel2(y1, y2, n1, n2, d), silent = TRUE )
  runtime <- proc.time() - runtime
  res$runtime <- runtime
  res$note <- paste("Chi-square approximation")

  if ( class(res) == "try-error" ) {
    res$info[1] <- 1e10
    res$info[2] <- 0
    res$p1 <- NA
    res$p2 <- NA
  }

  if ( R == 0 || res$info[1] == 1e+10 ) {
    res <- res

  } else if ( R == 1 ) {
    test <- as.numeric( res$info[1] )
    d <- dim(y1)[2]
    delta <- Rfast::james(y1, y2, R = 1)$info[3]
    stat <- as.numeric( test / delta )
    pvalue <- as.numeric( pchisq(stat, d, lower.tail = FALSE) )
    res$info[1] <- stat
    res$info[2] <- pvalue
    res$note <- paste("James corrected chi-square approximation")

  } else if ( R == 2 ) {
    test <- as.numeric( res$info[1] )
    d <- dim(y1)[2]
    dof <- Rfast::james(y1, y2, R = 2)$info[5]
    v <- dof + d - 1
    stat <- dof / (v * d) * test
    pvalue <- pf(stat, d, dof, lower.tail = FALSE)
    dof <- c(d, v - d + 1)
    res$info <- c(stat, pvalue, dof)
    names(res$info) <- c("statistic", "p-value", "numer df", "denom df")
    res$note <- paste("F approximation")
  } 
  
  res
}






