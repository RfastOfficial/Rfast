#[export]
ttest <- function(x, m, alternative = "unequal", logged = FALSE, conf = NULL) {

  n <- dim(x)[1] 
  xbar <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  stat <- sqrt(n) * ( xbar - m ) / s
  
  if ( alternative == "unequal" ) {
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE ) 
	} else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE) 
  } else if ( alternative == "greater" ) {
    pvalue <- pt( stat, lower.tail = FALSE, n - 1, log.p = logged )	
  } else if ( alternative == "less" ) {
    pvalue <- pt( stat, n - 1, log.p = logged )
  }
  res <- cbind(stat, pvalue)

  if ( !is.null(conf) ) {  
    a <- 1 - conf
    fac <- qt(1 - a/2, n - 1) * s/sqrt(n)
    mat <- cbind(xbar - fac, xbar + fac)
    colnames(mat) <- c( paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = "") )
    res <- cbind(res, mat)  
  }

  res
}


#[export]
ttest1 <- function(x, m, alternative = "unequal", logged = FALSE, conf = NULL) {

  n <- length(x) 
  xbar <- sum(x)/n
  s <- Rfast::Var(x, std = TRUE)
  stat <- sqrt(n) * ( xbar - m ) / s
  
  if ( alternative == "unequal" ) {
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE ) 
	} else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE) 
  } else if ( alternative == "greater" ) {
    pvalue <- pt( stat, lower.tail = FALSE, n - 1, log.p = logged )	
  } else if ( alternative == "less" ) {
    pvalue <- pt( stat, n - 1, log.p = logged )
  }
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")

  if ( !is.null(conf) ) {  
    a <- 1 - conf
    fac <- qt(1 - a/2, n - 1) * s/sqrt(n)
    mat <- c(xbar - fac, xbar + fac)
    names(mat) <- c( paste(a/2, "%", sep = ""), paste(1 - a/2, "%", sep = "") )
    res <- list(res = res, ci = mat) 
  }

  res
}


#[export]
ttest2 <- function(x, y, alternative = "unequal", paired = FALSE, logged = FALSE) {

  if ( !paired ) {
    n1 <- length(x)
    n2 <- length(y)
    m1 <- sum(x)/n1
    m2 <- sum(y)/n2
    f1 <- Rfast::Var(x) / n1
    f2 <- Rfast::Var(y) / n2
    fac <- f1 + f2
    dof <- fac^2 / ( f1^2 / (n1 - 1) + f2^2 / (n2 - 1) )
    stat <- ( m1 - m2 ) / sqrt(fac)
	
	if ( alternative == "unequal" ) {
      if ( logged ) {
        pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE ) 
	  } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE) 
    } else if ( alternative == "greater" ) {
      pvalue <- pt( stat, lower.tail = FALSE, dof, log.p = logged )	
    } else if ( alternative == "less" ) {
      pvalue <- pt( stat, dof, log.p = logged )
    }
    result <- c(stat, pvalue, dof)
    names(result) <- c("stat", "p-value", "dof")

  } else {
    n <- length(x)
    z <- x - y    
    m <- sum(z)/n
    s <- Rfast::Var(z, std = TRUE)
    stat <- sqrt(n) * m / s
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE )  
    } else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE )    	
    result <- c(stat, pvalue)
    names(result) <- c("stat", "p-value")
  }
  
  result
}


#[export]
ttests <- function(x, y = NULL, ina, alternative = "unequal", paired = FALSE, logged = FALSE, parallel = FALSE) {

  if ( !paired ) {

    if ( is.null(y) ) {
      x1 <- x[ ina == 1, ]
      x2 <- x[ ina == 2, ]
      n1 <- sum( ina == 1 )
      n2 <- length(ina) - n1
    } else {
      x1 <- x     ;    n1 <- dim(x1)[1]
	  x2 <- y     ;    n2 <- dim(x2)[1]
    }

    m1 <- Rfast::colmeans(x1, parallel = parallel)
    m2 <- Rfast::colmeans(x2, parallel = parallel)
    f1 <- Rfast::colVars(x1, parallel = parallel) / n1
    f2 <- Rfast::colVars(x2, parallel = parallel) / n2
    fac <- f1 + f2
    dof <- fac^2 / ( f1^2 / (n1 - 1) + f2^2 / (n2 - 1) )
    stat <- ( m1 - m2 ) / sqrt(fac)
	
	if ( alternative == "unequal" ) {
      if ( logged ) {
        pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE ) 
	  } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE) 
    } else if ( alternative == "greater" ) {
      pvalue <- pt( stat, lower.tail = FALSE, dof, log.p = logged )	
    } else if ( alternative == "less" ) {
      pvalue <- pt( stat, dof, log.p = logged )
    }
	
    #if ( logged ) {
    #  pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE )
    #} else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )  
    result <- cbind(stat, pvalue, dof)

  } else {
    n <- dim(x)[1]
    if ( is.null(y) ) {
      z <- x[ ina == 1, ] - x[ ina == 2, ]
    } else  z <- x - y    
    m <- Rfast::colmeans(z, parallel = parallel)
    s <- Rfast::colVars(z, std = TRUE, parallel = parallel)
    stat <- sqrt(n) * m / s
    if ( logged ) {
      pvalue <- log(2) + pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE )  
    } else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE )    	
    result <- cbind(stat, pvalue)
  }

  result
}


#[export]
ttests.pairs <- function(x, logged = FALSE) {
  
  n <- dim(x)[1]
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x) / n
  fac <- outer(s, s, "+")
  down <- outer( s^2/(n - 1), s^2/(n - 1), "+" )
  dof <- fac^2 / down 
  stat <- outer(m, m, "-") / sqrt(fac)
  if ( logged ) {
    pvalue <- log(2) + pt( abs(stat), dof, lower.tail = FALSE, log.p = TRUE )
  } else  pvalue <- 2 * pt( abs(stat), dof, lower.tail = FALSE )

  if ( is.null( colnames(x) ) ) {
    rownames(stat) <- rownames(pvalue) <- rownames(dof) <- 
    colnames(stat) <- colnames(pvalue) <- colnames(dof) <- paste("Var", 1:dim(x)[2])
  } else {
    rownames(stat) <- rownames(pvalue) <- rownames(dof) <- 
    colnames(stat) <- colnames(pvalue) <- colnames(dof) <- colnames(x)
  }

  list(stat = stat, pvalue = pvalue, dof = dof)  
}


#[export]
allttests <- function (x, y = NULL, ina, logged = FALSE) {
    if (is.null(y)) {
        x1 <- x[ina == 1, ]
        x2 <- x[ina == 2, ]
        n1 <- sum(ina == 1)
        n2 <- length(ina) - n1
    }
    else {
        x1 <- x
        n1 <- dim(x1)[1]
        x2 <- y
        n2 <- dim(x2)[1]
    }
    m1 <- Rfast::colmeans(x1)
    m2 <- Rfast::colmeans(x2)
    f1 <- Rfast::colVars(x1)/n1
    f2 <- Rfast::colVars(x2)/n2
    fac <- outer(f1, f2, "+")
    down <- outer(f1^2/(n1 - 1), f2^2/(n2 - 1), "+")
    dof <- fac^2/down
    difa <- outer(m1, m2, "-")
    stat <- difa/sqrt(fac)
    if (logged) {
        pvalue <- log(2) + pt(abs(stat), dof, lower.tail = FALSE, 
            log.p = TRUE)
    }
    else pvalue <- 2 * pt(abs(stat), dof, lower.tail = FALSE)
    if (is.null(colnames(x))) {
        rownames(stat) <- rownames(pvalue) <- rownames(dof) <- colnames(stat) <- colnames(pvalue) <- colnames(dof) <- paste("Var", 
            1:dim(x)[2])
    }
    else {
        rownames(stat) <- rownames(pvalue) <- rownames(dof) <- colnames(stat) <- colnames(pvalue) <- colnames(dof) <- colnames(x)
    }
    list(stat = stat, pvalue = pvalue, dof = dof)
}


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

   z1 <- sample(z1, R * n1, replace = TRUE)
   dim(z1) <- c(n1, R)
   z2 <- sample(z2, R * n2, replace = TRUE)
   dim(z2) <- c(n2, R)
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
   

#[export]
exact.ttest2 <- function(x, y) {
  n1 <- length(x)
  n2 <- length(y)
  n <- n1 + n2
  z <- c(x, y)
  m1 <- sum(x)/n1
  m2 <- sum(y)/n2
  fac <- Rfast::Var(x)/n1 + Rfast::Var(y)/n2
  tobs <- abs(m1 - m2)/sqrt(fac)
  st <- n1 * m1 + n2 * m2
  st2 <- sum(z^2)
  perms <- Rfast::comb_n(z, n1)
  sx <- Rfast::colsums(perms)
  x2 <- Rfast::colsums(perms^2)
  v1 <- (x2 - sx^2/n1)/(n1 - 1)
  pm1 <- sx / n1

  pm2 <- (st - sx) 
  y2 <- st2 - x2
  v2 <- (y2 - pm2^2/n2)/(n2 - 1) 
  pm2 <- pm2/n2
  tp <-  abs(pm1 - pm2)/sqrt(v1/n1 + v2/n2)
  res <- c( dim(perms)[2], tobs, sum(tp > tobs) / dim(perms)[2] )
  names(res) <- c("permutations", "stat", "p-value" )
  res
}



#[export]
percent.ttest <- function(x, y, logged = FALSE) {
  n1 <- length(x)
  n2 <- length(y)
  p1 <- sum(x)/n1
  p2 <- sum(y)/n2
  j22 <- p2 * (1 - p2) * n2
  j11 <- j22 + p1 * (1 - p1) * n1
  vb <- j11 / (j22 * j11 - j22^2)
  dof <- n1 + n2 - 2
  phi <- sum( (x - p1)^2/p1/(1 - p1) ) + sum( (y - p2)^2/p2/(1 - p2) )
  phi <- phi / dof
  b <- log( p2/(1 - p2) ) - log( p1/(1 - p1) )
  stat <- b^2 / (vb * phi)
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  res <- c(phi, stat, pval)
  names(res) <- c("phi", "stat", "p-value")
  res
}


#[export]
percent.ttests <- function(x, y, logged = FALSE) {
  n1 <- dim(x)[1]
  n2 <- dim(y)[1]
  p1 <- Rfast::colmeans(x)
  p2 <- Rfast::colmeans(y)
  j22 <- p2 * (1 - p2) * n2
  j11a <- p1 * (1 - p1) 
  j11 <- j22 + j11a * n1
  vb <- j11 / (j22 * j11 - j22^2)
  dof <- n1 + n2 - 2
  phi <- Rfast::rowsums( (t(x) - p1)^2 / j11a ) + Rfast::rowsums( (t(y) - p2)^2 /j22 ) * n2
  phi <- phi / dof
  b <- log( p2/(1 - p2) ) - log( p1/(1 - p1) )
  stat <- b^2 / (vb * phi)
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  res <- cbind(phi, stat, pval)
  colnames(res) <- c("phi", "stat", "p-value")
  res
}