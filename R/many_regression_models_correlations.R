#[export]
allbetas <- function(y, x, pvalue = FALSE, logged = FALSE) {

  n <- dim(x)[1]
  denom <- n - 1
  my <- sum(y) / n
  m <- Rfast::colmeans(x)
  r <- ( Rfast::eachcol.apply(x, y) - n * my * m )/denom

  sx <- Rfast::colVars(x, suma = n * m)
  be <- r / sx
  a <- my - be * m

  if ( !pvalue )  {
    result <- cbind(a, be)
    if ( is.null( colnames(x) ) ) {
      rownames(result) <- paste("X", 1:ncol(x), sep = "" )
    } else    rownames(result) <- colnames(x)

  } else {
    sy <- Rfast::Var(y)
    rho2 <- r^2 / (sx * sy)
    dof <- n - 2
    stat <- r * dof / (1 - rho2)
    pvalue <- pf( stat, 1, n - 2, lower.tail = FALSE, log.p = logged)
    result <- cbind(a, be, r / sqrt( sy * sx ), stat, pvalue)
    if ( is.null( colnames(x) ) ) {
      rownames(result) <- paste("X", 1:ncol(x), sep = "" )
    } else  rownames(result) <- colnames(x)
  }

  result
}


#[export]
correls <- function(y, x, type = "pearson", a = 0.05, rho = 0) {
  n <- length(y)

  if (type == "pearson") {
    r <- as.vector( cor(y, x) ) ## the correlation value between y and all the xs
    zh0 <- 0.5 * log( (1 + rho) / (1 - rho) )  ## Fisher's transformation for Ho
    zh1 <- 0.5 * log( (1 + r) / (1 - r) )  ## Fisher's transformation for H1
    se <- 1/sqrt(n - 3)  ## standard error for Fisher's transformation of Ho

  } else if (type == "spearman") {
    r <- as.vector( cor(Rfast::Rank(y), Rfast::colRanks(x) ) ) ## the correlation value between y and all the xs
    zh0 <- 0.5 * log( (1 + rho) / (1 - rho) )  ## Fisher's transformation for Ho
    zh1 <- 0.5 * log( (1 + r) / (1 - r) )  ## Fisher's transformation for H1
    se <- 1.029563 / sqrt(n - 3)  ## standard error for Fisher's transformation of Ho
  }

  test <- (zh1 - zh0) / se ## test statistic
  pvalue <-  2 * pt( abs(test), n - 3, lower.tail = FALSE )  ## p-value
  b1 <- zh1 - qt(1 - a/2, n - 3) * se
  b2 <- zh1 + qt(1 - a/2, n - 3) * se
  ca <- cbind(b1 ,b2)
  ela <- exp( 2 * ca )
  ci <- ( ela - 1 ) / ( ela + 1 )  ## confidence intervals
  res <- cbind(r, ci, test, pvalue)
  colnames(res) <- c( 'correlation', paste( c( a/2 * 100, (1 - a/2) * 100 ), "%", sep = ""), 'z-stat', 'p-value' )

  if ( is.null(colnames(x)) ) {
    rownames(res) <- paste("X", 1:ncol(x), sep = "")
  } else  rownames(res) <- colnames(x)

  res
}


#[export]
expregs <- function(y, x, di, tol = 1e-09, logged = FALSE) {
    
  dm <- dim(x)
  a0 <- sum(y) /dm[1] 
  di2 <- 1 - di
  sdi <-  - sum(di)
  lam <-  - sdi / sum(y)
  ini <-  - 2 * sdi * log(lam) + 2 * sdi
  con <- log( dm[1] )
  D <- dm[2]
  lik <- numeric(D)
  yyhat0 <- y * a0
  dera20 <-  - sum(yyhat0)  
  dera0 <- sdi - dera20
  
  for (j in 1:D) {
    X <- x[, j] 
    X2 <- X^2
    s2 <-  - sum(di * X)
    yyhat <- yyhat0
    dera2 <-  dera20  
    dera <- dera0
    derab <-  - sum(X * yyhat) 
    derb <-  s2 - derab
    derb2 <-  - sum(X2 * yyhat) 
    aold <- c(1/a0, 0)
    anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
  
    while ( sum( abs(anew - aold) ) > tol ) {
      aold <- anew
      a <- anew[1]     ;      b <- anew[2] 
      yyhat <- y * exp( -a - b * X)
      dera2 <-  - sum(yyhat)  
      dera <- sdi - dera2
      derab <-  - sum(X * yyhat) 
      derb <-  s2 - derab
      derb2 <-  - sum(X2 * yyhat) 
      anew <- aold - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 ) 
    }
    a <- anew[1]    ;     b <- anew[2]
    lik[j] <-  a * sdi + b * s2 + dera2
  }
    
  stat <-  2 * lik - ini
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}


#[export]
gammaregs <- function(y, x, tol = 1e-07, logged = FALSE, maxiters = 100) {
  
    dm <- dim(x)
    n <- dm[1]
    D <- dm[2]
    devi <- phi <- numeric(D)
    ind <- 1:D 

    sx <- n
    ly <- log(y)
    sy <- sum(y)
    be <- sum(ly)/n
    m <- exp(-be)
    der2 <- sy * m
    d1 <- sum(ly) + n * log(m) - der2
    der <-  - der2 + sx
    be <- be - der/der2
    m <- exp(-be)
    d2 <- sum(ly) + n * log(m) - der2
    while ( abs(d2 - d1) > tol) {
      d1 <- d2
      der2 <- sy * m
      der <-  - der2 + sx
      be <- be - der/der2
      m <- exp(-be)
      d2 <- sum(ly) + n * log(m) - der2
    }
    common <- y * m
    phi <- sum( (common - 1)^2 )/(n - 1)
    ini <-  - 2 * d2 - 2 * n
    b1 <- be 
   
    sx <- Rfast::colsums(x)
    x2 <- x^2 
    for (j in ind) {
      d1 <- ini
      X <- x[, j]
      sxj <- sx[j]
      x2j <- x2[, j]
      m <- exp(-b1)
      dera2 <- sy * m
      dera <-  - dera2 + n
      derab <- sum(common * X)
      derb <-  - derab + sxj
      derb2 <- sum(common * x2j)
      be <- c(b1, 0) - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
      a <- be[1]   ;    b <- be[2]
      m <- exp(- a - b * X)
      com <- y * m
      d2 <-  - 2 * sum( log(com) ) + 2 * sum(com) - 2 * n
      i <- 2
      while ( abs(d2 - d1) > tol  &  i < maxiters ) {
        i <- i + 1
	    d1 <- d2
        dera2 <- sum(com)
        dera <-  - dera2 + n
        derab <- sum(com * X)
        derb <-  - derab + sxj
        derb2 <- sum(y * x2j * m)
        be <- be - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
        a <- be[1]   ;    b <- be[2]
        m <- exp(- a - b * X)
        com <- y * m
        d2 <-  - 2 * sum( log(com) ) + 2 * sum(com) - 2 * n
      }
      com <- y * m
      phi[j] <- sum( (com - 1)^2 ) / (n - 2)
      devi[j] <-  d2
    }     
    stat <- (ini - devi)/phi
    pvalue <- pf(stat, 1, n - 2, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue)
}


#[export]
geom.regs <- function(y, x, tol = 1e-07, type = 1, logged = FALSE, parallel = FALSE, maxiters = 100) {
  mod <- .Call(Rfast_geom_regs,y, x, tol,logged, type,parallel,maxiters)
  colnames(mod) <- c("stat", "pval")
  mod
}


#[export]
groupcorrels <- function(y, x, type ="pearson", ina) {
  p <- dim(x)[2]
  a <- Rfast::sort_unique(ina)
  d <- length(a) 
  res <- matrix(nrow = d, ncol = p)
  for (i in 1:d) res[i, ] <- cor(y[ina == i], x[ina == i, ], method = type)  
  rownames(res) <- a
  res
}


#[export]
invgauss.regs <- function(y, x, tol = 1e-08, logged = FALSE, maxiters = 100) {
  
  dm <- dim(x)
  n <- dm[1]
  sy <- sum(1 / y)
  ly <- log(y)
  y2 <- 2 * y
  me <- sy/n
  D <- dm[2]
  phi <- dev <- numeric(D)
  con <- log(n)
  d0 <- n/Rfast::invgauss.mle(y)$param[2]
  be <- Rfast::allbetas( ly, x )

  for (j in 1:D) {
    X <- x[, j]
	X2 <- X^2 
    a <- be[j, 1]   ;   b <- be[j, 2]
    mi <- exp(a + b * X)
    mi2 <- mi^2
    dera <- sum( (mi - y) / mi2)
    dera2 <- sum( (y2 - mi ) / mi2 )
    derb <- sum(X * (mi - y)/mi2 )
    derb2 <- sum(X2 * (y2 - mi) / mi2 )
    derab <- sum(X * (y2 - mi) / mi2 )
    be1 <- c(a, b)
    be2 <- be1 - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
	i <- 2  
    while ( sum( abs(be1 - be2) ) > tol  &  i < maxiters ) {
	  i <- i + 1
      be1 <- be2
      a <- be1[1]   ;   b <- be1[2]
      mi <- exp(a + b * X)
	  mi2 <- mi^2
      dera <- sum( (mi - y) / mi2)
      dera2 <- sum( (y2 - mi ) / mi2 )
      derb <- sum(X * (mi - y)/mi2 )
      derb2 <- sum(X2 * (y2 - mi) / mi2 )
      derab <- sum(X * (y2 - mi) / mi2 )
      be2 <- be1 - c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
    }
    
    dev[j] <- n * ( me - sum( 1 / mi) / n )
    phi[j] <- sum( (y - mi)^2 / ( mi^3) )/(n - 2) 
  } 
  stat <- (d0 - dev)/phi
  pvalue <- pf(stat, 1, n -2, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}


#[export]
logistic_only <- function(x,y,tol = 1e-09,b_values = FALSE) {
	if(b_values){
		x<-.Call(Rfast_logistic_only_b,x,y,tol)
		rownames(x)<-c("deviance","constant","slope")
		return (x)
	}
	.Call(Rfast_logistic_only,x,y,tol)
}


#[export]
multinom.regs <- function(y, x, tol = 1e-08, logged = FALSE, parallel = FALSE, maxiters = 100) {
  if ( !is.numeric(y) )  y <- as.numeric(y)
  mod <- .Call("Rfast_multinom_regs",PACKAGE = "Rfast",y, x, tol,logged, parallel, maxiters )
  colnames(mod) <- c("stat", "pval")
  mod
}


#[export]
mvbetas <- function (y, x, pvalue = FALSE) {

    n <- dim(y)[1]
    my <- Rfast::colmeans(y)
    mx <- sum(x)/n
    sx <- ( sum(x^2) - sum(x)^2 / n ) / (n - 1)
    denom <- n - 1    
    r <- ( Rfast::eachcol.apply(y, x) - n * mx * my ) /denom
    be <- r/sx
    a <- my - be * mx
    
    if ( !pvalue ) {
      result <- cbind(a, be)
      if ( is.null( colnames(y) ) ) {
         rownames(result) <- paste("Y", 1:ncol(y), sep = "")
      }  else  rownames(result) <- colnames(x)

    } else {
      sy <- Rfast::colVars(y, suma = n * my, std = TRUE)
      rho <- r/(sqrt(sx) * sy)
      sqdof <- sqrt(n - 2)
      ta <- rho * sqdof/sqrt(1 - rho^2)
      pvalue <- 2 * pt(abs(ta), n - 2, lower.tail = FALSE)
      result <- cbind(a, be, rho, pvalue)
      if ( is.null( colnames(y) ) ) {
        rownames(result) <- paste("Y", 1:ncol(y), sep = "")
      } else  rownames(result) <- colnames(x)
    }
    result
}


#[export]
normlog.regs <- function (y, x, tol = 1e-08, logged = FALSE, parallel = FALSE, 
    maxiters = 100)  {
    
   dm <- dim(x)
   n <- dm[1]
   ly <- log(y + 0.1)
   con <- Rfast::Var(y) * (n - 1)
   be <- allbetas( ly, x )
   mod <- .Call(Rfast_normlog_regs, y, x, be, con, tol, logged, parallel, maxiters)
   colnames(mod) <- c("stat", "pvalue")
   mod
}


#[export]
poisson_only <- function(x,y,tol = 1e-09,b_values=FALSE) {
	if(b_values){
		x<-.Call(Rfast_poisson_only_b,x,y,sum(y*log(y),na.rm=TRUE),tol)
		rownames(x)<-c("log-lik","constant","slope")
		return (x)
	}
	.Call(Rfast_poisson_only,x,y,sum(y*log(y),na.rm=TRUE),tol)
}

#[export]
prop.regs <- function(y, x, varb = "quasi", tol = 1e-09, logged = FALSE, maxiters = 100) {
  stat <- .Call(Rfast_prop_regs,x,y, tol, varb,maxiters)
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pval)
}


#[export]
qpois.regs <- function (x, y, tol = 1e-09, logged = FALSE) {
  ylogy <- sum(y * log(y), na.rm = T)
  stat <- .Call(Rfast_qpois_regs,x, y, tol, ylogy, mean(y))
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  ret <- cbind(stat, pval)
  colnames(ret) <- c("stat", "pval")
  ret
}


#[export]
colanovas <- function(y, x, logged = FALSE) {
  if ( is.data.frame(x) ) {
	x <- Rfast::data.frame.to_matrix(x)
  }
  n <- dim(x)[1]
  b <- sum(y)^2/n
  sy2 <- sum(y^2)
  a <- .Call('Rfast_col_anovas', PACKAGE = 'Rfast',y,x)
  k <- Rfast::colrange(x,cont=FALSE)
  mst <- (a - b) / (k - 1)
  mse <- (sy2 - a) / (n - k)
  fa <- mst / mse
  pvalue <- pf(fa, k - 1, n - k, lower.tail = FALSE, log.p = logged)
  tab <- cbind(fa, pvalue)
  colnames(tab) <- c("F stat", "p-value")
  tab
}


#[export]
regression <- function (x, y, poia = NULL, logged = FALSE) {
    if (is.matrix(x)) {
        n <- length(y)
        rho <- as.vector(cor(y, x))
        sqdof <- sqrt(n - 2)
        stat <- rho * sqdof/sqrt(1 - rho^2)
        if (logged) {
            pvalue <- log(2) + pt(abs(stat), n - 2, lower.tail = FALSE, 
                log.p = TRUE)
        }
        else pvalue <- 2 * pt(abs(stat), n - 2, lower.tail = FALSE)
    }
    else {
        
        if ( is.null(poia) )  poia <- Rfast::which.is(x)
        if ( length(poia) == 0 ) {
            n <- length(y)
            rho <- as.vector(cor(y, x))
            sqdof <- sqrt(n - 2)
            stat <- rho * sqdof/sqrt(1 - rho^2)
            if (logged) {
                pvalue <- log(2) + pt(abs(stat), n - 2, lower.tail = FALSE, 
                  log.p = TRUE)
            }
            else pvalue <- 2 * pt(abs(stat), n - 2, lower.tail = FALSE)
        }
        else {
            dm <- dim(x)
            n <- dm[1]
            D <- dm[2]
            stat <- numeric(D)
            pvalue <- numeric(D)
            n <- length(y)
            rho <- as.vector(cor(y, x[, -poia]))
            sqdof <- sqrt(n - 2)
            stat[-poia] <- rho * sqdof/sqrt(1 - rho^2)
            if (logged) {
                pvalue[-poia] <- log(2) + pt(abs(stat[-poia]), 
                  n - 2, lower.tail = FALSE, log.p = TRUE)
            }
            else pvalue[-poia] <- 2 * pt(abs(stat[-poia]), n - 
                2, lower.tail = FALSE)
            mod <- Rfast::colanovas(y, x[, poia, drop = FALSE], 
                logged = logged)
            stat[poia] <- mod[, 1]
            pvalue[poia] <- mod[, 2]
        }
    }
    cbind(stat, pvalue)
}


#[export]
spml.regs <- function(y, x, tol = 1e-07, logged = FALSE, maxiters = 100, parallel = FALSE) {
  res <- .Call(Rfast_spml_regs, as.matrix(y),x,tol,logged,maxiters,parallel)
  colnames(res) <- c("statistic", "p-value")
  res
}


#[export]
univglms <- function (y, x, oiko = NULL, logged = FALSE) {
    dm <- dim(x)
    n <- dm[1]
    d <- dm[2]
    if (is.null(oiko)) {
	    y <- as.numeric(y)  
        if ( length( Rfast::sort_unique(y) ) == 2 ) {
            oiko = "binomial"
        }
        else if ( sum( Rfast::Round(y) - y) == 0 ) {
            oiko = "poisson"
        }
        else oiko = "normal"
    }
    if (oiko == "binomial") {
        p <- sum(y)/n
        ini <-  - 2 * ( n * p * log(p) + (n - n * p) * log(1 - p) )
        mod <- Rfast::logistic_only(x, y)
        stat <- ini - mod
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    }
    else if (oiko == "poisson") {
        m <- sum(y)/n
        ini <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
        mod <- Rfast::poisson_only(x, y)
        stat <- ini - mod
        pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
    }
    else if (oiko == "normal") {
        rho <- as.vector(cor(y, x))
        sqdof <- sqrt(n - 2)
        stat <- rho * sqdof/sqrt(1 - rho^2)
        if (logged) {
            pval <- log(2) + pt(abs(stat), n - 2, lower.tail = FALSE, log.p = TRUE)
        }
        else pval <- 2 * pt(abs(stat), n - 2, lower.tail = FALSE)
    }
    else if (oiko == "quasipoisson") {
        m <- sum(y)/n
        ini <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
        mod <- Rfast::quasi.poisson_only(x, y)
        stat <- (ini - mod[1, ])/mod[2, ]
        pval <- pf(stat, 1, n - 2, lower.tail = FALSE, log.p = logged)
    }
    result <- cbind(stat, pval)
    colnames(result) <- c("stat", "pvalue")
    if (is.null(colnames(x))) {
        rownames(result) <- paste("Var", 1:d, sep = "")
    }
    else row.names(result) <- colnames(x)
    result
}


#[export]
univglms2 <- function (y, x, oiko = NULL, logged = FALSE) {
    dm <- dim(x)
    n <- dm[1]
    d <- dm[2]
    if (is.null(oiko)) {
	    y <- as.numeric(y)  
        if (length(Rfast::sort_unique(y)) == 2) {
            oiko <- "binomial"
        }
        else if (sum(round(y) - y) == 0) {
            oiko <- "poisson"
        }
        else oiko <- "normal"
    }
    if (oiko == "binomial") {
        poia <- Rfast::which.is(x)
        x <- Rfast::data.frame.to_matrix(x)
        if (length(poia) == 0) {
            p <- sum(y)/n
            ini <-  - 2 * ( n * p * log(p) + (n - n * p) * log(1 - p) )
            mod <- Rfast::logistic_only(x, y)
            stat <- ini - mod
            pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
        }
        else {
            stat <- numeric(d)
            pval <- numeric(d)
            p <- sum(y)/n
            ini <-  - 2 * ( n * p * log(p) + (n - n * p) * log(1 - p) )
            mod <- Rfast::logistic_only(x[, -poia, drop = FALSE], y)
            stat[-poia] <- ini - mod
            pval[-poia] <- pchisq(stat[-poia], 1, lower.tail = FALSE, log.p = logged)
            x <- cbind(y, x[, poia] - 1)
            dc <- Rfast::colrange(x, cont = FALSE)
            mod <- Rfast::g2tests(x, 2:(length(poia) + 1), 1, dc)
            stat[poia] <- mod$statistic
            pval[poia] <- pchisq(mod$statistic, mod$df, lower.tail = FALSE, log.p = logged)
        }
        result <- cbind(stat, pval)
    }
    else if (oiko == "poisson") {
        poia <- Rfast::which.is(x)
        x <- Rfast::data.frame.to_matrix(x)
        if (length(poia) == 0) {
            m <- sum(y)/n
            ini <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
            mod <- Rfast::poisson_only(x, y)
            stat <- ini - mod
            pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
        }
        else {
            stat <- numeric(d)
            pval <- numeric(d)
            m <- sum(y)/n
            ini <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
            mod <- Rfast::poisson_only(x[, -poia, drop = FALSE], y)
            stat[-poia] <- ini - mod
            pval[-poia] <- pchisq(stat[-poia], 1, lower.tail = FALSE, log.p = logged)
            d1 <- numeric(d)
            k <- Rfast::colrange(x[, poia, drop = FALSE], cont = FALSE)
            for (i in poia) {
                ina <- x[, i]
                ni <- tabulate(ina)
                ni <- ni[ni > 0]
                si <- rowsum(y, ina)
                mi <- si/ni
                d1[i] <- sum(si * log(mi))
            }
            d0 <- n * m * log(m)
            stat[poia] <- 2 * d1[poia] - 2 * d0
            pval[poia] <- pchisq(stat[poia], k - 1, lower.tail = FALSE, log.p = logged)
        }
        result <- cbind(stat, pval)
    }
    else if (oiko == "quasipoisson") {
        poia <- Rfast::which.is(x)
        x <- Rfast::data.frame.to_matrix(x)
        if (length(poia) == 0) {
            m <- sum(y)/n
            ini <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
            mod <- Rfast::quasi.poisson_only(x, y)
            stat <- (ini - mod[1, ])/mod[2, ]
            pval <- pf(stat, 1, n - 2, lower.tail = FALSE, log.p = logged)
        }
        else {
            stat <- numeric(d)
            pval <- numeric(d)
            sy <- sum(y)
            m <- sy/n
            ini <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * log(m)
            if (length(poia) < d) {
                mod <- Rfast::quasi.poisson_only(x[, -poia, drop = FALSE], y)
                stat[-poia] <- (ini - mod[1, ])/mod[2, ]
                pval[-poia] <- pf(stat[-poia], 1, n - 2, lower.tail = FALSE, log.p = logged)
            }
            d0 <- 2 * sy * log(sy/n)
            for (i in poia) {
                ina <- x[, i]
                ni <- tabulate(ina)
                ni <- ni[ni > 0]
                k <- length(ni)
                si <- rowsum(y, ina)
                mi <- si/ni
                d1 <- Rfast::colsums(si * log(mi))
                up <- (2 * d1 - d0)/(k - 1)
                yi2 <- rowsum(y^2, ina)/mi
                phi <- ( Rfast::colsums(yi2) - sy ) / (n - k)
                stat[i] <- up/phi
                pval[i] <- pf(stat[i], k - 1, n - k, lower.tail = FALSE, log.p = logged)
            }
        }
        result <- cbind(stat, pval)
    }
    else if (oiko == "normal") {
        result <- Rfast::regression(x, y, logged = logged)
    }
    colnames(result) <- c("stat", "pvalue")
    result
}


#[export]
quasi.poisson_only <- function(x,y,tol = 1e-09, maxiters = 100) {
	.Call(Rfast_quasi_poisson_only,x,y,sum(y*log(y),na.rm=TRUE),tol,maxiters)
}

  

