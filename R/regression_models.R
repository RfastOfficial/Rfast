#[export]
ar1 <- function(y, method = "cmle") {
  N <- length(y)
  if ( method == "cmle" ) {
    dera2 <- N - 1
    derab <- sum( y[-N] )
    derb2 <- sum( y[-N]^2 )
    dera <- derab
    derb <- sum( y[-N] * y[-1])   
    cphi <- c( derb2 * dera - derab * derb, - derab * dera + dera2 * derb ) / ( dera2 * derb2 - derab^2 )
    s <- sum( (y[-1] - cphi[1] - cphi[2] * y[-N])^2 )/dera2
    param <- c(cphi, s)
    names(param) <- c("constant", "phi", "sigma")
  } else if ( method == "yw" ) {   
    m <- sum(y)/N
    z <- y - m
    phi <- sum( z[-N] * z[-1] ) / sum(z^2)
    sigma <- (1 - phi^2) * sum(z^2)/(N - 2)
    param <- c(m, phi, sigma)
    names(param) <- c("mean", "phi", "sigma")
  } else if ( method == "ols" ) {
    m <- sum(y)/N
    z <- y - m
    y1 <- z[-1] 
    x1 <- z[-N]
    phi <- as.vector( cov(y1, x1) )/ Var(x1)
    constant <-  mean(y1) - phi * mean(x1)
    sigma <- sum( (y1 - constant - phi * x1 )^2 ) / (N - 1)
    param <- c(constant, phi, sigma)
    names(param) <- c("constant", "phi", "sigma")
  }
  param 
}


#[export]
bc <- function(x, low = -1, up = 1) {
  lx <- log(x)
  slx <- sum( lx )
  n2 <-  - 0.5 * length(x)
  vlx <- Rfast::Var(lx)
  ell <- function(lambda, x, vlx, slx, n2) {
    if ( abs(lambda) < 1e-12 ) {
      s <- vlx
    } else {
      y <- x^lambda
      s <- Rfast::Var(y) / lambda^2
    }
    n2 * log(s) + lambda * slx
  }
  optimise(ell, c(low, up), x = x, vlx = vlx, slx = slx, n2 = n2, 
           maximum = TRUE, tol = 1e-06)$maximum
}  


#[export]
gammareg <- function(y, x, tol = 1e-07, maxiters = 100) {
  X <- model.matrix( y~., data.frame(x) )
  sx <- Rfast::colsums(X)
  dm <- dim(X)
  n <- dm[1]  ;  p <- dm[2]

  mod <- Rfast::gammacon(y, tol = tol) 
  m <- exp(mod$be)
  be <- c( mod$be, numeric(p - 1) )
  d1 <-  - sum(y / m) + n * mod$be
  con <- y / m
  der <- Rfast::eachcol.apply(X, con) - sx
  der2 <-  - crossprod( con * X, X) 
  be <- be - solve(der2, der)
  xbe <-  X %*% be
  m <- as.vector( exp( xbe ) )
  con <- y / m
  d2 <-  - sum(con) - sum(xbe)

  i <- 2
  while ( abs(d1 - d2) > tol  & i < maxiters) {
    i <- i + 1
    d1 <- d2
    der <- Rfast::eachcol.apply(X, con) - sx
    der2 <-  - crossprod( con * X, X) 
    be <- be - solve(der2, der)
    xbe <-  X %*% be
    m <- as.vector( exp( xbe ) )
    con <- y / m
    d2 <-  - sum(con) - sum(xbe)
  }
  phi <- sum( (con - 1)^2 ) / (n - p)
  devi <-  - 2 * sum( log(con) ) + sum(con) - n
  info <- c(i, devi, phi)
  names(info) <- c("iters", "deviance", "phi")
  names(be) <- colnames(X)
  list(info = info, be = be)
}



#[export]
glm_logistic <- function (x, y, full = FALSE, tol = 1e-09,maxiters = 100) {
    x <- model.matrix(y ~ ., data.frame(x))
    mod <- .Call(Rfast_glm_logistic, x, y,tol,maxiters)
	rownames(mod$be) <- colnames(x)
    res <- list(be = mod$be, devi = mod$deviance)
    if (full) {
        be <- mod$be
        se <- chol2inv(chol(mod$der2))
        se <- sqrt(diag(se))
        wald <- be/se
        pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
        info <- cbind(be, se, wald, pval)
        colnames(info) <- c("estimate", "std error", "Wald stat", 
            "p-value")
        rownames(info) <- colnames(x)
        res <- list(info = info, devi = mod$deviance)
    }
    res["iter"] <- mod$iter
    res
}


#[export]
glm_poisson <- function (x, y, full = FALSE,tol = 1e-09) {
    x <- model.matrix(y ~ ., data.frame(x))
    mod <- .Call(Rfast_glm_poisson, x, y, sum(y * log(y), na.rm = TRUE),tol)
	rownames(mod$be) <- colnames(x)
    res <- list(be = mod$be, devi = mod$deviance)
    if (full) {
        be <- mod$be
        se <- chol2inv(chol(mod$L2))
        se <- sqrt(diag(se))
        wald <- be/se
        pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
        info <- cbind(be, se, wald, pval)
        colnames(info) <- c("estimate", "std error", "Wald stat", 
            "p-value")
        rownames(info) <- colnames(x)
        res <- list(info = info, devi = mod$deviance)
    }
    res
}


#[export]
invgauss.reg <- function (y, x, tol = 1e-07, maxiters = 100) {
    X <- model.matrix(~., data.frame(x))
    dm <- dim(X)
    ly <- log(y)
    y2 <- 2 * y
    n <- dm[1]
    sy <- sum(1/y)
    be1 <- solve( crossprod(X), crossprod(X, ly) )
    mi <- exp(as.vector(X %*% be1))
    f <- Rfast::colsums((mi - y)/mi^2 * X)
    f2 <- crossprod((y2/mi^2 - 1/mi) * X, X)
    be2 <- be1 - solve(f2, f)
    i <- 2
    while ( sum(abs(be1 - be2)) > tol & i < maxiters ) {
        i <- i + 1
        be1 <- be2
        mi <- as.vector(exp(X %*% be1))
        f <- Rfast::colsums((mi - y)/mi^2 * X)
        f2 <- crossprod((y2/mi^2 - 1/mi) * X, X)
        be2 <- be1 - solve(f2, f)
    }
    lambda <- 1/(sy/n - sum(1/mi)/n)
    loglik <- n/2 * log(lambda/2/pi) - 1.5 * sum(ly) - lambda/2 * 
        sum( (-y/mi^2 + sy/n) )
    deviance <- n/lambda
    phi <- sum((y - mi)^2/(mi^3))/(n - 2)
	names(be2) <- colnames(X)
    list(iters = i, loglik = loglik, deviance = deviance, phi = phi, 
        be = be2)
}


#[export]
lmfit <- function(x, y, w = NULL) {
  if ( is.null(w) ) {
    be <- solve( crossprod(x), crossprod(x, y) )
  } else  be <- solve( crossprod(x, w * x), crossprod(x, w * y) )
  e <- as.vector( y - x %*% be )
  rownames(be) <- colnames(x)
  list(be = be, residuals = e)
}


#[export]
logistic.cat1 <- function(y, x, logged = FALSE) {
  N <- Rfast::Table(y, x) 
  cj <- log( N[2, ] / N[1, ] )
  be <- c(cj[1], cj[-1] - cj[1])
  n <- length(y)
  Nj <- Rfast::colsums( N )
  pj <- N[2, ] / Nj
  p <- sum(N[2, ]) / n
  D0 <-  - 2 * n * sum( p * log(p) + (1 - p) * log(1 - p) )
  D1 <-  - 2 * sum( Nj * pj * log(pj) + Nj * (1 - pj) * log(1 - pj) )
  se <- pj * (1 - pj) * Nj 
  se <- c( 1/se[1], 1 / se[-1] + 1/se[1] )
  se <- sqrt(se)
  stat <- be/se
  pval <- pchisq(stat^2, 1, lower.tail = FALSE)
  mat <- cbind(be, se, stat^2, pval)
  colnames(mat) <- c("estimate", "std. error", "Wald", "p-value")
  devs <- c(D0, D1)
  stat <- D0 - D1
  pvalue <- pchisq(stat, length(be) - 1, lower.tail = FALSE, log.p = logged)
  devs <- c(D0, D1, stat, pvalue)
  names(devs) <- c("null deviance", "residual deviance", "stat", "p-value")
  list(info = mat, devs = devs)
}
  

#[export]
multinom.reg <- function(y, x, tol = 1e-07, maxiters = 50) {
  Y <- Rfast::design_matrix(y)[, -1]
  if ( !is.matrix(Y) ) {
    mod <- Rfast::glm_logistic(x, y, tol = tol, maxiters = maxiters)
    res <- list(iters = mod$iter, loglik = 0.5 * mod$devi, be = mod$be)
  } else {
    X <- model.matrix(y~., data.frame(x))
    p <- dim(X)[2] 
    dm <- dim(Y)
    n <- dm[1]
    d <- dm[2] 
    m <- Rfast::colmeans(Y)
    b0 <- Rfast::Log(m / (1 - m) )  
    b1 <- matrix( c(b0, numeric(p * d - d) ), nrow = p, ncol = d, byrow = TRUE)
    e <- Y - rep(m, rep(n, d) )
    id <- matrix(1:c(p * d), ncol = d)
    der <- numeric(d * p)
    der2 <- matrix(0, p * d, p * d)
    for (i in 1:d) { 
      der[id[, i]] <- Rfast::colsums( e[, i] * X )
      for (j in i:d) {
        if (i != j) {
          der2[id[, i], id[, j]] <- der2[id[, j], id[, i]] <-  - crossprod(m[i] * m[j] * X, X) 
        } else  der2[id[, i], id[, i]] <- crossprod(m[i] * (1 - m[i]) * X, X)
      }
   }  
    b2 <- b1 + solve(der2, der)
    k <- 2
    res <- try(while (sum(abs(b2 - b1)) > tol  &  k < maxiters) {
      k <- k + 1
      b1 <- b2
      m1 <- exp(X %*% b1) 
      m <- m1 / (rowsums(m1) + 1)
      e <- Y - m
      for (i in 1:d) { 
      der[id[, i]] <- Rfast::colsums( e[, i] * X )
        for (j in i:d) {
          if (i != j) {
            der2[id[, i], id[, j]] <- der2[id[, j], id[, i]] <-   - crossprod(m[, i] * m[, j] * X, X) 
          } else  der2[id[, i], id[, i]] <- crossprod(m[, i] * (1 - m[, i]) * X, X)
        }
      }  
      b2 <- b1 + solve(der2, der)
    }, silent = TRUE)
    if (inherits(res,"try-error"))   b2 <- b1
    colnames(b2) <- paste("Y", 1:d, sep = "")
    rownames(b2) <- colnames(X)
    Y1 <- design_matrix(y, ones = FALSE)
    m <- cbind(1, m1) 
    m <- m / rowsums(m)
    loglik <-  - sum(Y1 * log(Y1/m), na.rm = TRUE )
    res <- list(iters = k, loglik = loglik, be = b2)
  }
  res
}


#[export]
normlog.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
  x <- model.matrix(y ~ ., data.frame(x))
  mod <- .Call(Rfast_normlog_reg,y, x, tol, maxiters)
  rownames(mod$be) <- colnames(x)
  mod
}


#[export]
poisson.cat1 <- function(y, x, logged = FALSE) {
  ni <- tabulate(x)
  si <- rowsum(y, x)
  mi <- si/ni
  n <- length(y)
  k <- length(ni)
  be <- log(mi)
  be <- c(be[1], be[-1] - be[1])
  s <- matrix(0, k, k)
  diag(s) <- si
  s[1, ] <- si
  s[, 1] <- si
  s[1, 1] <- sum(si)
  se <- sqrt( diag( chol2inv(chol(s)) ) )
  #########
  cl1 <- sum( si * log(mi) )
  cl0 <- s[1, 1] * log( s[1, 1]/n )  
  stat <- 2 * cl1 - 2 * cl0
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  ########
  stat <- be^2/se^2
  pval <- pchisq(stat, 1, lower.tail = FALSE)
  mat <- cbind(be, se, stat, pval)
  colnames(mat) <- c("estimate", "std. error", "Wald", "p-value")
  list(info = mat, res = res)
}
 

#[export]
prop.reg <- function (y, x, varb = "quasi", tol = 1e-09, maxiters = 100) {
    X <- model.matrix(~., data.frame(x))
    L <- .Call(Rfast_prop_reg, X, y, tol, maxiters)
    der2 <- L$der2
    u <- y - as.vector(L$p)
    be <- as.vector(L$be)
	names(be) <- colnames(X)
    Ainv <- spdinv(der2)
    phi <- NULL
    if (varb == "quasi") {
        B <- crossprod(u * X)
        vb <- Ainv %*% B %*% Ainv
    }
    else {
        dm <- dim(X)
        dof <- dm[1] - dm[2]
        phi <- sum((y - L$p)^2/L$p/(1 - L$p))/dof
        vb <- phi * Ainv
    }
    info <- cbind(be, sqrt(diag(vb)), be^2/diag(vb))
    info <- cbind(info, pchisq(info[, 3], 1, lower.tail = FALSE))
    rownames(info) <- colnames(X)
    colnames(info) <- c("Estimate", "Std. error", "Wald", "p-value")
    list(iters = L$i, varb = vb, phi = phi, info = info)
}


################################
#### Spatial median regression
#### Tsagris Michail 10/2014
#### Biman Chakraborty (2003) On multivariate quantile regression
#### Journal of Statistical Planning and Inference
#### http://www.stat.nus.edu.sg/export/sites/dsap/research/documents/tr01_2000.pdf
################################

#[export]
spatmed.reg <- function(y, x, tol = 1e-07) {
  x <- model.matrix(y ~ ., data.frame(x) )
  p <- dim(x)[2]
  d <- dim(y)[2]

  B1 <- solve(crossprod(x), crossprod(x, y) )
  est <- y - x %*% B1
  ww <- sqrt( Rfast::rowsums(est^2) )
  ## ww <- sqrt( Rfast::rowsums( est^2 ) )
  z <- x / ww
  B2 <- solve( crossprod(z, x), crossprod(z, y) )
  i <- 2
  while ( sum( abs(B2 - B1) ) > tol ) {
    i <- i + 1
    B1 <- B2
    est <- y - x %*% B1
    ww <- sqrt( Rfast::rowsums(est^2) )
	## ww <- sqrt( Rfast::rowsums( est^2 ) )
    ela <- which( ww == 0 )
    z <- x / ww
    if ( length(ela) > 0 )  z[ela, ] <- 0
    B2 <- solve( crossprod(z, x), crossprod(z, y) )
  }
  be <- B2
  
  if ( is.null(colnames(y)) ) {
    colnames(be) <- paste("Y", 1:d, sep = "")
  } else  colnames(be) <- colnames(y)
  rownames(be)  <- colnames(x)

  list(iters = i, be = be)
}


#[export]
spml.reg <- function(y, x, tol = 1e-07, seb = FALSE, maxiters = 100) {
  x <- model.matrix(~., data.frame(x) )
  y <- as.matrix(y)
  mod <- .Call(Rfast_spml_reg, y, x, tol, seb, maxiters)
  colnames(mod$be) <- c("Cosinus of y", "Sinus of y")
  rownames(mod$be) <- colnames(x)   
 
  if ( seb ) {
    colnames(mod$seb) <- c("Cosinus of y", "Sinus of y")
    rownames(mod$seb) <- colnames(x)    
  } else   mod$seb <- NULL
  
  mod
}


#[export]
weib.reg <- function (y, x, tol = 1e-07, maxiters = 100) {
    x <- model.matrix(y ~ ., data.frame(x))
    mod <- .Call(Rfast_weib_reg, y, x, tol, 
        maxiters)
    rownames(mod$be) <- colnames(x)
    list(iters = mod$iters, loglik = mod$loglik, shape = mod$shape, be = mod$be)
}


#[export]
qpois.reg <- function (x, y, full = FALSE, tol = 1e-09, maxiters = 100) {
    x <- model.matrix(y ~ ., data.frame(x))
    mod <- .Call(Rfast_qpois_reg, x, y, sum(y * log(y), na.rm = TRUE), 
        tol,maxiters)
	rownames(mod$be) <- colnames(x)	
    res <- list(be = mod$be, devi = mod$deviance, varb = mod$phi * 
        spdinv(mod$L2), phi = mod$phi)
    if (full) {
        be <- mod$be
        varb <- mod$phi * spdinv(mod$L2)
        info <- cbind(be, sqrt(diag(varb)), be^2/diag(varb))
        info <- cbind(info, pchisq(info[, 3], 1, lower.tail = FALSE))
        rownames(info) <- colnames(x)
        colnames(info) <- c("Estimate", "Std. error", "Wald", 
            "p-value")
        res <- list(info = info, devi = mod$deviance, varb = varb, 
            phi = mod$phi)
    }
    res
}

