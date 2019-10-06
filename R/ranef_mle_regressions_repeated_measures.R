
#[export]
colrint.regbx <- function(y, x, id) {
  mod <- Rfast::lmfit(cbind(1, x), y) 
  be <- mod$be
  e <- mod$residuals
  N <- dim(y)[1]
  n <- Rfast::sort_unique.length(id)
  d <- N / n
  f <- 1 - 1/d
  seid <- rowsum(e^2, id)
  myid <- rowsum(y, id)/ d
  my <- Rfast::colsums(myid) / n
  com <- ( t(myid) - my )^2
  sml <- Rfast::rowsums( t(seid) - d * com )/N/f
  dml <- Rfast::rowsums(com) / n /sml  - 1/ d
  tauml <- dml * sml
  neg <- which( tauml < 0)
  if ( length(neg) >0 ) {
    sml[neg] <- sml[neg] + tauml[neg]
    tauml[neg] <- 0
  } 
  com <- rowsum(e, id)
  ranef <- tauml /( tauml + sml/d ) * t(com) / d
  sz <- Rfast::colsums(seid)
  sz2 <- Rfast::colsums(com^2) 
  loglik <- n * d * log(sml) + n * log1p(d * tauml / sml) + sz/sml - tauml / (sml^2 + d * tauml * sml) * sz2
  loglik <-  -0.5 * loglik - n * d/2 * log(2 * pi) 
  dev <-  -2 * loglik
  info <- cbind(tauml, sml, loglik, dev, dev + 4 * log(N) )
  colnames(info) <- c("sigma_tau", "sigma_errors", "log-likelihood", "deviance", "BIC")
  list(info = info, be = be, ranef = ranef)
}


#[export]
colvarcomps.mle <- function (x, id, ranef = FALSE, tol = 1e-08, maxiters = 100, 
    parallel = FALSE) {
    if ( Rfast::Var( tabulate(id) ) != 0) {
        res <- .Call("Rfast_colrint_mle", PACKAGE = "Rfast", 
            x, id, ranef, tol, maxiters, parallel)
    }
    else {
        N <- length(id)
        n <- Rfast::sort_unique.length(id)
        d <- N/n
        f <- 1 - 1/d
        mxid <- rowsum(x, id)/d
        mx <- Rfast::colsums(mxid)/n
        ex <- t(t(x) - mx)
        sexid <- rowsum(ex^2, id)
        co <- t((t(mxid) - mx)^2)
        sx2 <- Rfast::colsums(co)
        sml <- Rfast::colsums(sexid - d * co)/N/f
        dml <- sx2/n/sml - 1/d
        tauml <- dml * sml
        poia <- which(dml < 0)
        sml[poia] <- sml[poia] + tauml[poia]
        tauml[poia] <- 0
        sx <- Rfast::colsums(sexid)
        sx2 <- sx2 * d^2
        loglik <- N * log(sml) + n * log1p(d * tauml/sml) + sx/sml - 
            tauml/(sml^2 + d * tauml * sml) * sx2
        loglik <- -0.5 * loglik - N/2 * log(2 * pi)
        info <- cbind(tauml, sml, loglik)
        colnames(info) <- c("sigma_tau", "sigma_errors", "log-likelihood")
        res <- list(info = info)
        if (ranef) {
            com <- t(rowsum(ex, id)/d)
            ranef <- tauml/(tauml + sml/d) * com
            res <- list(info = info, ranef = t(ranef))
        }
    }
    res
}


#[export]
colvarcomps.mom <- function(x, id, parallel = FALSE) {
  k <- Rfast::sort_unique.length(id) 
  ni <- tabulate(id)
  ni <- ni[ni > 0]
  sam <- length(ni)
  n <- dim(x)[1]
  sx2 <- Rfast::colsums(x^2, parallel = parallel)
  m <- rowsum(x, id)
  a <- Rfast::colsums(m^2/ni)
  b <- Rfast::colsums(m)^2/n
  df1 <- k - 1
  df2 <- n - k   
  mst <- (a - b)/df1       ### /(k - 1)
  mse <- (sx2 - a)/df2   ##  /(n - k)
  fa <- mst / mse
  ranvar <- (mst - mse)/sam  
  ranvar[ranvar <= 0 ] <- 0
  rat <- ranvar / (ranvar + mse)
  L <- (fa *  qf(0.025, df2, df1) - 1) / sam
  U <- (fa * qf(0.975, df2, df1) - 1 ) / sam
  l <- L / (1 + L)
  u <- U / (1 + U)
  res <- cbind(ranvar, mse, rat, l, u)
  colnames(res) <- c("ranvar", "MSE", "ratio", "2.5% lower", "97.5% upper")
  res
}


#[export]
rint.mle <- function(x, ina, ranef = FALSE, tol = 1e-09, maxiters = 100) {
  res <- .Call(Rfast_rint_mle,x,ina,ranef,tol,maxiters)
  names(res$info) <- c("sigma_tau", "sigma_errors", "log-lik")
  res      
}  


#[export]
rint.reg <- function (y, x, id, tol = 1e-08, ranef = FALSE, maxiters = 100) {
    x <- cbind(1, x)
    mod <- .Call(Rfast_rint_reg, x, y, id, tol, ranef, maxiters)
    names(mod$info) <- c("iters", "sigma_tau", "sigma_errors", "log-lik", 
        "deviance", "BIC")
    mod
}

# rint.reg <- function (y, x, id, tol = 1e-08) {
#     X <- cbind(1, x)
#     dm <- dim(X)
#     n <- dm[1]
#     p <- dm[2]
#     xx <- crossprod(X)
#     sx <- rowsum(X, id)
#     sxy <- crossprod(X, y)
#     sy <- as.vector( rowsum(y, id) )
#     ni <- tabulate(id)
#     mx <- sx/ni
#     my <- sy/ni
#     funa <- function(d, n, ni, S, hi2) sum(log1p(ni * d)) + n * 
#         log(S - d * sum(ni^2 * hi2/(1 + ni * d)))
#     mod <- .lm.fit(X, y)
#     b1 <- mod$coefficients
#     S <- sum(mod$residuals^2)
#     hi2 <- (my - mx %*% b1)^2
#     mod <- optimise(funa, c(0, 50), n = n, ni = ni, S = S, hi2 = hi2, 
#         tol = tol)
#     d <- mod$minimum
#     b2 <- solve(xx - d * crossprod(sx/(1 + ni * d), sx), sxy - 
#         d * crossprod(sx, sy/(1 + ni * d)))
#     i <- 2
#     while (sum(abs(b2 - b1)) > tol) {
#         i <- i + 1
#         b1 <- b2
#         S <- sum((y - X %*% b1)^2)
#         hi2 <- (my - mx %*% b1)^2
#         mod <- optimise(funa, c(0, 50), n = n, ni = ni, S = S, 
#             hi2 = hi2, tol = tol)
#         d <- mod$minimum
#         b2 <- solve(xx - d * crossprod(sx/(1 + ni * d), sx), 
#             sxy - d * crossprod(sx, sy/(1 + ni * d)))
#     }
#     se <- (S - d * sum(ni^2 * hi2/(1 + ni * d)))/n
#     tau <- d * se
#     loglik <- -0.5 * mod$objective - 0.5 * n * (log(2 * pi) - 
#         log(n) + 1)
#     dev <- -2 * loglik
#     seb <- sqrt(diag(solve(xx - d * crossprod(sx/(1 + ni * d), 
#         sx)) * se))
#     info <- c(i, tau, se, loglik, dev, dev + (p + 2) * log(n))
#     names(info) <- c("iters", "sigma_tau", "sigma_errors", "log-lik", 
#         "deviance", "BIC")
#     list(info = info, be = b2, seb = seb)
# }


#[export]
rint.regbx <- function(y, x, id) {
  mod <- Rfast::lmfit(cbind(1, x), y) 
  be <- mod$be
  e <- mod$residuals
  N <- length(y)
  n <- Rfast::sort_unique.length(id)
  d <- N / n
  f <- 1 - 1/d
  seid <- rowsum(e^2, id)
  myid <- rowsum(y, id)/ d
  my <- sum(myid) / n
  com <- (myid - my)^2
  sml <- sum(seid - d * com )/N/f
  dml <- sum(com) / n /sml  - 1/ d
  tauml <- dml * sml
  com <- rowsum(e, id)
  ranef <- tauml /( tauml + sml/d ) * com / d
  sz <- sum(seid)
  sz2 <- sum(com^2) 
  loglik <- n * d * log(sml) + n * log1p(d * tauml / sml) + sz/sml - tauml / (sml^2 + d * tauml * sml) * sz2
  loglik <-  -0.5 * loglik - n * d/2 * log(2 * pi) 
  dev <-  -2 * loglik
  info <- c(tauml, sml, loglik, dev, dev + 4 * log(N) )
  names(info) <- c("sigma_tau", "sigma_errors", "log-likelihood", "deviance", "BIC")
  list(info = info, be = be, ranef = ranef)
}


#[export]
rint.regs <- function(y, x, id, tol = 1e-08, logged = FALSE, parallel = FALSE, maxiters = 100) {
  mod <- .Call(Rfast_rint_regs, x, y, id, tol, logged, parallel, maxiters)
  colnames(mod) < c("stat", "pval")
  mod
}


#[export]
rm.anova <- function(y, logged = FALSE) {
  dm <- dim(y)  
  d <- dim(y)[2]
  n <- dim(y)[1]
  ina <- rep(1:n, each = d)
  xi <- rep(1:d, n)
  yi <- Rfast::rowmeans(y)
  yj <- Rfast::colmeans(y)
  yt <- mean(yi)
  sst <- n * sum( (yj - yt)^2 )
  yi <- rep(yi, each = d)
  yj <- rep(yj, n)
  ssr <- sum( (as.vector( t(y) ) - yi - yj + yt)^2)
  dft <- d - 1
  dfs <- n - 1
  dfr <- dft * dfs
  mst <- sst/dft
  msr <- ssr/dfr
  stat <- mst/msr
  pvalue <- pf(stat, dft, dfr, lower.tail = FALSE, log.p = logged)
  c(stat, pvalue)
}


#[export]
rm.anovas <- function(y, x, logged = FALSE) {
  d <- length(x)
  n <- dim(y)[1]/d
  ina <- rep(1:n, each = d)
  xi <- rep(x, n)
  yi <- rowsum(y, ina)/d
  yj <- rowsum(y, xi)/n
  yt <- Rfast::colmeans(y)
  sst <- n * Rfast::colsums( (yj - yt)^2 )
  #sss <- d * colsums( (yi - yt)^2 )
  yi <- rep(yi, each = d)
  yj <- rep(yj, n)
  ssr <- Rfast::colsums( (y - yi - yj + yt)^2 )
  dft <- d - 1
  dfs <- n - 1
  dfr <- dft * dfs
  mst <- sst / dft
  #mss <- ssr / dfs
  msr <- ssr / dfr
  stat <- mst/msr
  pvalue <- pf(stat, dft, dfr, lower.tail = FALSE, log.p = logged) 
  cbind(stat, pvalue)
}


#[export]
rm.lines <- function(y, x, logged = FALSE) {
  z <- x - mean( x )
  xi <- z / sum( z^2 )
  d <- length(x)
  n <- dim(y)[1]/d
  xi <- rep(xi, n)
  ina <- rep(1:n, each = d)
  be <- rowsum(xi * y, ina)
  s <- Rfast::colVars(be, std = TRUE)
  stat <- sqrt(n) * Rfast::colmeans(be) / s
  if ( logged ) {
    pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE, log.p = TRUE ) 
  } else  pvalue <- 2 * pt( abs(stat), n - 1, lower.tail = FALSE) 
  cbind(stat, pvalue)
}


#[export]
varcomps.mle <- function(x, ina, tol = 1e-09) {
	mat <- .Call(Rfast_varcomps_mle, x, ina, Rfast::sort_unique.length(ina), tol)
	syina <- mat$syina
	mat <- mat$mat
	d <- mat[4]
	mat <- mat[1:3]
	ranef <- mat[1]/(mat[1] + mat[2]/d)*syina/d  ## mat[1]*syina/(mat[1]*d+mat[2])
	names(mat) <- c("sigma_tau", "sigma_errors", "log-likelihood")
	list(info = mat, ranef = ranef)
}


#[export]
varcomps.mom <- function(x, ina) {
  ni <- tabulate(ina)
  k <- length(ni)
  n <- length(x)
  sx2 <- sum(x^2)
  m <- Rfast::group(x, ina)
  a <- sum(m^2/ni)
  b <- sum(m)^2/n
  df1 <- k - 1
  df2 <- n - k   
  mst <- (a - b)/df1       ### /(k - 1)
  mse <- (sx2 - a)/df2   ##  /(n - k)
  fa <- mst / mse
  ranvar <- (mst - mse)/k  
  ranvar[ranvar <= 0 ] <- 0
  rat <- ranvar / (ranvar + mse)
  L <- (fa *  qf(0.025, df2, df1) - 1) / k
  U <- (fa * qf(0.975, df2, df1) - 1 ) / k
  l <- L / (1 + L)
  u <- U / (1 + U)
  res <- c(ranvar, mse, rat, l, u)
  names(res) <- c("ranvar", "MSE", "ratio", "2.5% lower", "97.5% upper")
  res
}
