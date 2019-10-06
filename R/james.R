#[export]
james <- function(y1, y2, a = 0.05, R = 1) {
  p <- dim(y1)[2]  ## dimensionality of the data
  n1 <- dim(y1)[1]   ;   n2 <- dim(y2)[1]  ## sample sizes
  n <- n1 + n2  ## the total sample size
  ybar1 <- Rfast::colmeans(y1)  ## sample mean vector of the first sample
  ybar2 <- Rfast::colmeans(y2)  ## sample mean vector of the second sample
  dbar <- ybar2 - ybar1  ## difference of the two mean vectors
  mesoi <- rbind(ybar1, ybar2)
  rownames(mesoi) <- c("Sample 1", "Sample 2")

  if ( is.null(colnames(y1)) ) {
    colnames(mesoi) <- paste("X", 1:p, sep = "")
  } else  colnames(mesoi) <- colnames(y1)
  A1 <- Rfast::cova(y1)/n1
  A2 <- Rfast::cova(y2)/n2
  V <- A1 + A2  ## covariance matrix of the difference
  Vinv <- Rfast::spdinv(V)   ## same as solve(V), but faster
  test <- sum( dbar %*% Vinv * dbar )
  b1 <- Vinv %*% A1
  b2 <- Vinv %*% A2
  trb1 <- sum( diag(b1) )
  trb2 <- sum( diag(b2) )

  if (R <= 1) {
    ## James test
    A <- 1 + ( trb1^2/(n1 - 1) + trb2^2/(n2 - 1) ) / (2 * p)
    B <- ( sum(b1 * b1) / (n1 - 1) + sum(b2 * b2)/(n2 - 1) + 0.5 * trb1 ^ 2/ (n1 - 1) + 0.5 * trb2^2/(n2 - 1) ) / (p * (p + 2))
    x2 <- qchisq(1 - a, p)
    delta <- (A + B * x2)
    twoha <- x2 * delta  ## corrected critical value of the chi-square
    pvalue <- pchisq(test/delta, p, lower.tail = FALSE)  ## p-value of the test statistic
    info <- c(test, pvalue, delta, twoha)
    names(info) <- c("test", "p-value", "correction", "corrected.critical")
    note <- paste("James test")
    result <- list(note = note, mesoi = mesoi, info = info)

  } else if (R == 2) {
    ## MNV test
    low <- ( sum( b1^2 ) + trb1^2 ) / n1 + ( sum( b2^2 ) + trb2^2 ) / n2
    v <- (p + p^2) / low
    test <- as.numeric( (v - p + 1) / (v * p) * test )  ## test statistic
    crit <- qf(1 - a, p, v - p + 1)  ## critical value of the F distribution
    pvalue <- pf(test, p, v - p + 1, lower.tail = FALSE)  ## p-value of the test statistic
    info <- c(test, pvalue, crit, p, v - p + 1)
    names(info) <- c("test", "p-value", "critical", "numer df", "denom df")
    note <- paste("MNV variant of James test")
    result <- list(note = note, mesoi = mesoi, info = info)
  } 
  result
}
