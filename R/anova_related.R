#[export]
anova1 <- function (x, ina, logged = FALSE) {
  ina <- as.numeric(ina)
  k <- max(ina)
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  n <- length(x)
  sx2 <- sum(x^2)
  m <- Rfast::group(x, ina)
  a <- sum(m^2/ni)
  b <- sum(m)^2/n
  mst <- (a - b) / (k - 1)
  mse <- (sx2 - a) / (n - k)
  fa <- mst / mse
  pvalue <- pf(fa, k - 1, n - k, lower.tail = FALSE, log.p = logged)
  tab <- c(fa, pvalue)
  names(tab) <- c("F stat", "p-value")
  tab
}


#[export]
anovas <- function (x, ina, logged = FALSE) {
    ina <- as.numeric(ina)
    k <- max(ina)
    ni <- tabulate(ina)
	ni <- ni[ni > 0]
    n <- dim(x)[1]
    sx2 <- Rfast::colsums(x^2)
    m <- rowsum(x, ina)
    a <- Rfast::colsums(m^2/ni)
    b <- Rfast::colsums(m)^2/n
    mst <- (a - b) / (k - 1)
    mse <- (sx2 - a) / (n - k)
    fa <- mst / mse
    pvalue <- pf(fa, k - 1, n - k, lower.tail = FALSE, log.p = logged)
    tab <- cbind(fa, pvalue)
    colnames(tab) <- c("F value", "p-value")
    if (!is.null(colnames(x))) 
        rownames(tab) <- colnames(x)
    tab
}



#[export]
ancova1 <- function(y, ina, x, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  a <- length(ni) 
  N <- length(ina)
  sy <- sum(y)
  sx <- sum(x)
  com <- sy^2/N
  com2 <- sx * sy / N
  syy <- sum(y^2) - com
  sxx <- sum(x^2) - sx^2/N
  sxy <- sum(x * y) - com2
  tyy <- sum( Rfast::group(y, ina)^2/ni ) - com
  txx <- sum( Rfast::group(x, ina)^2/ni ) - sx^2/N
  txy <- sum( Rfast::group(x, ina)/ni * Rfast::group(y, ina) ) - com2
  eyy <- syy - tyy
  exx <- sxx - txx
  exy <- sxy - txy
  b <- exy / exx
  sse <- eyy - exy^2/exx
  sse2 <- syy - sxy^2/sxx
  dof <- N - a - 1
  mse <- sse / dof
  ftreat <- (sse2 - sse)/(a - 1) / mse
  fb <- sxy^2/sxx/mse
  pvaltreat <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  pvalb <- pf(fb, 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- c(ftreat, fb, pvaltreat, pvalb)
  names(res) <- c("Ftreat", "Fbeta", "pvalue-treat", "pvalue-beta")
  res
} 


#[export]
ancovas <- function(y, ina, x, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  a <- length(ni) 
  N <- length(ina)
  sy <- Rfast::colsums(y)
  sx <- sum(x)
  com <- sy^2/N
  com2 <- sx * sy / N
  syy <- Rfast::colsums(y^2) - com
  sxx <- sum(x^2) - sx^2/N
  sxy <- Rfast::colsums(x * y) - com2
  tyy <- Rfast::colsums( rowsum(y, ina)^2/ni ) - com
  txx <- sum( Rfast::group(x, ina)^2/ni ) - sx^2/N
  txy <- Rfast::colsums( Rfast::group(x, ina)/ni * rowsum(y, ina) ) - com2
  eyy <- syy - tyy
  exx <- sxx - txx
  exy <- sxy - txy
  b <- exy / exx
  sse <- eyy - exy^2/exx
  sse2 <- syy - sxy^2/sxx
  dof <- N - a - 1
  mse <- sse / dof
  ftreat <- (sse2 - sse)/(a - 1) / mse
  fb <- sxy^2/sxx/mse
  pvaltreat <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  pvalb <- pf(fb, 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- cbind(ftreat, fb, pvaltreat, pvalb)
  colnames(res) <- c("Ftreat", "Fbeta", "pvalue-treat", "pvalue-beta")
  res
} 



#[export]
block.anova <- function(x, treat, block, logged = FALSE) {
  a <- Rfast::sort_unique.length(treat) 
  b <- Rfast::sort_unique.length(block)
  N <- length(x)
  com <- sum(x)^2/N
  sst <- sum(x^2) - com
  ssa <- sum( Rfast::group(x, treat)^2 ) / b - com
  ssb <- sum( Rfast::group(x, block)^2 ) / a - com
  dof <- (a - 1) * (b - 1) 
  mse <- (sst - ssa - ssb) / dof
  ftreat <- ssa / (a - 1)/mse
  pval <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- c(ftreat, pval)
  names(res) <- c("F-stat", "p-value")
  res
}  



#[export]
block.anovas <- function(x, treat, block, logged = FALSE) {
  a <- Rfast::sort_unique.length(treat)
  b <- Rfast::sort_unique.length(block)
  N <- dim(x)[1]
  com <- Rfast::colsums(x)^2/N
  sst <- Rfast::colsums(x^2) - com
  ssa <- Rfast::colsums( rowsum(x, treat)^2 ) / b - com
  ssb <- Rfast::colsums( rowsum(x, block)^2 ) / a - com
  dof <- (a - 1) * (b - 1) 
  mse <- (sst - ssa - ssb) / dof
  ftreat <- ssa / (a - 1)/mse
  pval <- pf(ftreat, a - 1, dof, lower.tail = FALSE, log.p = logged) 
  res <- cbind(ftreat, pval)
  colnames(res) <- c("F-stat", "p-value")
  res
}  



#[export]
colpoisson.anovas <- function (y, x, logged = FALSE) {
  dm <- dim(x)
  n <- dm[1]
  sy <- sum(y)
  d0 <- 2 * sy * log(sy/n)
  stat <- numeric(dm[2])
  pvalue <- numeric(dm[2])
  phi <- numeric(dm[2])
  ######  C++
  for ( i in 1:dm[2] ) {
    ina <- x[, i]
    ni <- tabulate(ina)
    k <- length(ni)
    si <- rowsum(y, ina)
    mi <- si/ni
    d1 <- sum( si * log(mi) )
    stat[i] <- 2 * d1 - d0
    pvalue[i] <- pchisq(stat[i], k - 1, lower.tail = FALSE, log.p = logged)
  }
  ######
  res <- cbind(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}



#[export]
colquasipoisson.anovas <- function (y, x, logged = FALSE) {
  dm <- dim(x)
  n <- dm[1]
  sy <- sum(y)
  d0 <- 2 * sy * log(sy/n)
  d1 <- numeric(dm[2])
  yi2 <- y^2
  stat <- numeric(dm[2])
  pvalue <- numeric(dm[2])
  phi <- numeric(dm[2])
  ######  C++
  for ( i in 1:dm[2] ) {
    ina <- x[, i]
    ni <- tabulate(ina)
    k <- length(ni)
    si <- Rfast::group(y, ina)
    mi <- si/ni
    d1 <- sum(si * log(mi))
    up <- 2 * d1 - d0
    yi2 <- Rfast::group(yi2, ina)/mi
    phi[i] <- (sum(yi2) - sy) / (n - k)    
    stat[i] <- up / phi[i]
    pvalue[i] <- pf(stat[i], k - 1, n - k, lower.tail = FALSE, log.p = logged)
  }
  #######
  res <- cbind(stat, pvalue, phi)
  names(res) <- c("stat", "p-value", "phi")
  res
}



#[export]
geom.anova <- function (y, ina, type = 1, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  if (type == 1) {   
     si <- Rfast::group(y, ina)
     pi <-  1/(1 + si/ni )
     ell1 <- sum( ni * log(pi) + si * log(1 - pi) )
     prob <- 1/(1 + sum(si)/sum(ni))
     ell0 <- 2 * sum( dgeom(y, prob, log = TRUE) )
     stat <- 2 * ell1 - ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }  else {
     si <- Rfast::group(y, ina)
     pi <- ni/si
     ell1 <- sum( ni * log(pi) + (ni/pi - ni) * log(1 - pi) )
     n <- sum(ni)
     prob <- n/sum(si)
     ell0 <- n * log(prob) + (n/prob - n) * log(1 - prob)
     stat <- 2 * ell1 - 2 * ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}



#[export]
geom.anovas <- function (y, ina, type = 1, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  n <- sum(ni)
  if (type == 1) {   
     ell1 <- 0
     si <- rowsum(y, ina)
     pi <-  1/(1 + si/ni )
     ell1 <- Rfast::colsums( ni * log(pi) + si * log(1 - pi) )
     prob <- 1/(1 + Rfast::colsums(si)/n )
     ell0 <- n * log(prob) + Rfast::colsums(si) * log(1 - prob)
     stat <- 2 * ell1 - 2 * ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }  else {
     si <- rowsum(y, ina)
     pi <- ni/si
     ell1 <- Rfast::colsums( ni * log(pi) + (ni/pi - ni) * log(1 - pi) )
     prob <- n/sum(si)
     ell0 <- n * log(prob) + (n/prob - n) * log(1 - prob)
     stat <- 2 * ell1 - 2 * ell0
     pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)     
  }
  res <- cbind(stat, pvalue)
  colnames(res) <- c("stat", "p-value")
  res
}


#[export]
poisson.anova <- function(y, ina, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  n <- sum(ni)
  si <- Rfast::group(y, ina)
  mi <- si/ni
  d1 <- sum( si * log(mi) )
  d0 <- sum(si) * log(sum( si)/n )  
  stat <- 2 * d1 - 2 * d0
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  res <- c(stat, pvalue)
  names(res) <- c("stat", "p-value")
  res
}


#[export]
poisson.anovas <- function(y, ina, logged = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  n <- sum(ni)
  si <- rowsum(y, ina)
  mi <- si/ni
  d1 <- Rfast::colsums( si * log(mi) )
  d0 <- Rfast::colsums(si) * log( Rfast::colsums(si)/n )  
  stat <- 2 * d1 - 2 * d0
  pvalue <- pchisq(stat, k - 1, lower.tail = FALSE, log.p = logged)
  res <- cbind(stat, pvalue)
  colnames(res) <- c("stat", "p-value")
  res
}


#[export]
quasipoisson.anova <- function (y, ina, logged = FALSE) {
    ni <- tabulate(ina)
	ni <- ni[ni > 0]
    k <- length(ni)
    n <- sum(ni)
    si <- Rfast::group(y, ina)
    sy <- sum(si)
    mi <- si/ni
    d1 <- sum(si * log(mi))
    d0 <- sy * log(sy/n)
    up <- ( 2 * d1 - 2 * d0 ) / (k - 1)
    yi2 <- Rfast::group(y^2, ina)/mi
    phi <- ( sum(yi2) - sy ) / (n - k)   
    stat <- up /  phi
    pvalue <- pf(stat, k - 1, n - k, lower.tail = FALSE, log.p = logged)
    res <- c(stat, pvalue, phi)
    names(res) <- c("stat", "p-value", "phi")
    res
}


#[export]
quasipoisson.anovas <- function (y, ina, logged = FALSE) {
    ni <- tabulate(ina)
	ni <- ni[ni > 0]
    k <- length(ni)
    n <- sum(ni)
    si <- rowsum(y, ina)
    sy <- Rfast::colsums(si)
    mi <- si/ni
    d1 <- Rfast::colsums(si * log(mi))
    d0 <- sy * log(sy/n)
    up <- ( 2 * d1 - 2 * d0 ) / (k - 1)
    yi2 <- rowsum(y^2, ina)/mi
    phi <- ( Rfast::colsums(yi2) - sy ) / (n - k)
    stat <- up / phi
    pvalue <- pf(stat, k - 1, n - k, lower.tail = FALSE, log.p = logged)
    cbind(stat, pvalue, phi)
}


#[export]
anova_propreg <- function(mod, poia = NULL) {
  b <- mod$info[, 1]
  vb <- mod$varb
  if ( is.null(poia) ) {
    poia <- 2:length(b)
    stat <- b[poia] %*% solve(vb[poia, poia], b[poia] )
  } else if ( length(poia) == 1 ){ 
    stat <- mod$info[poia, 3]
  } else  stat <- b[poia] %*% solve(vb[poia, poia], b[poia] )
  dof <- length(poia)
  pvalue <- pchisq(stat, dof, lower.tail = FALSE)
  res <- c(stat, pvalue, dof) 
  names(res) <- c("statistic", "p-value", "degrees of freedom")
  res
}


#[export]
anova_qpois.reg <- function(mod, poia = NULL) {
   b <- mod$info[, 1]
   vb <- mod$varb
   if (is.null(poia)) {
     poia <- 2:length(b)
     stat <- b[poia] %*% solve(vb[poia, poia], b[poia])
   } else if (length(poia) == 1) {
     stat <- mod$info[poia, 3]
   } else stat <- b[poia] %*% solve(vb[poia, poia], b[poia])
   dof <- length(poia)
   pvalue <- pchisq(stat, dof, lower.tail = FALSE)
   res <- c(stat, pvalue, dof)
   names(res) <- c("statistic", "p-value", "degrees of freedom")
   res
}


#[export]
anova_quasipois.reg <- function(mod0, mod1, n) {
   phi <- mod1$phi   
   df1 <- length(mod1$be) - length(mod0$be) 
   stat <- (mod0$devi - mod1$devi) / df1 / phi
   df2 <- n - length(mod1$be)
   pvalue <- pf(stat, df1, df2, lower.tail = FALSE)
   res <- c(stat, pvalue, df1, df2)
   names(res) <- c("statistic", "p-value", "Numerator dof", "Denominator dof")
   res
}


