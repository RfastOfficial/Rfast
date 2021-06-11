#[export]
bic.corfsreg <- function (y, x, tol = 2) {
    dm <- dim(x)
    n <- dm[1]
    p <- dm[2]
    con <- n * log(2 * pi) + n
    logn <- log(n)
    x <- Rfast::standardise(x, center = TRUE, scale = TRUE)
    y <- ( y - mean(y) ) / Rfast::Var(y)
    model <- NULL
    tool <- sela <- numeric( min(n, p) )
    tool[1] <- n * log( Var(y) * (n - 1) /n) + 2 * logn 
    sela[1] <- 0
	oop <- options(warn = -1)
    on.exit( options(oop) )
    yx <- Rfast::eachcol.apply(x, y)
    sel <- which.max(abs(yx))
    r <- yx[sel] / (n - 1)
    z <- x[, sel, drop = FALSE]
    x[, sel] <- 0
    model <- .lm.fit(as.matrix(z), y)
    tool[2] <- n * log(sum(model$residuals^2)/n) + 3 * logn 
    sela[2] <- sel
    k <- 2         
    while ( k < n - 19 & tool[k - 1] - tool[k] > tol ) {
      m <- n - 3 - k
      k <- k + 1
      #e1 <- .lm.fit(z, y)$residuals
      #e2 <- .lm.fit(z, x)$residuals
      res <- .lm.fit(z, cbind(y, x))$residuals
      e1 <- res[, 1]
      e2 <- res[, -1]
      yx.z <- cor(e1, e2)
      sel <- which.max( abs(yx.z) )
      z <- cbind(z, x[, sel])
      x[, sel] <- 0
      model <- .lm.fit(z, y)
      tool[k] <- n * log( sum( model$residuals^2 )/n ) + (k + 1) * logn
      sela[k] <- sel
    }  
    res <- cbind(sela[1:c(k-1)], tool[1:c(k-1)] + con)
    colnames(res) <- c("sel", "bic")
    res
}


#[export]
bic.fs.reg <- function(y, x, tol = 2, type = "logistic") {
  ret <- .Call(Rfast_bic_fs_reg,y, x, tol, type)
  colnames(ret) <- c("vars", "bic")
  ret
}


#[export]
bs.reg <- function(y, x, alpha = 0.05, type = "logistic") {
  .Call(Rfast_bs_reg,y, x, alpha, type)
}


#[export]
cor.fbed <- function(y, x, ystand = TRUE, xstand = TRUE, alpha = 0.05, K = 0) {
  quan <- 1 - alpha/2
  sig <- log(alpha)
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  ind <- 1:p
  sela <- NULL
  card <- 0
  sa <- NULL
  pva <- NULL

  runtime <- proc.time()
  if (xstand)   x <- Rfast::standardise(x)
  if (ystand)   y <- ( y - mean(y) ) / Rfast::Var(y, std = TRUE)
  oop <- options(warn = -1)
  on.exit( options(oop) )
  yx <- Rfast::eachcol.apply(x, y) / (n - 1)
  n.tests <- p
  stat <- abs( 0.5 * log( (1 + yx) / (1 - yx) ) * sqrt(n - 3) ) 
  critvalue <- qt(quan, n - 3)
  s <- which(stat > critvalue)
  
  if ( length(s) > 0 ) {
    sel <- which.max(stat) 
    sela <- sel
    s <- s[ - which(s == sel) ]
    r <- yx[sel]
    pv <- log(2) + pt(stat[sel], n - 2, lower.tail = FALSE, log.p = TRUE)   
    pva <- pv 
    sa <- stat[sel] 
    stat <- numeric(p)
    z <- x[, sel]
    
    if ( length(s) > 0 ) {
      xz <- as.vector( cor(z, x[, ind[s] ]) )
      n.tests <- n.tests + length( ind[s] )
      yx.z <- abs( ( yx[ ind[s] ] - xz * r ) / sqrt(1 - xz^2) / sqrt(1 - r^2) ) 
      stat[ ind[s] ] <- abs( 0.5 * log( (1 + yx.z) / (1 - yx.z) ) * sqrt(n - 4) )
      critvalue <- qt(quan, n - 4)
      s <- which( stat > critvalue )
      if ( length(s) > 0 ) {
        sel <- which.max(stat)
        pv <- log(2) + pt(stat[sel], n - 4, lower.tail = FALSE, log.p = TRUE) 
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pv)
        sela <- c(sela, sel[sel > 0] )
        s <- s[ - which(s == sel) ]
        z <- cbind(z, x[, sel])
      }  ## end if ( length(s) > 0 )
      ######################
      while ( sum(s > 0) > 0 )  {
        stat <- numeric(p)
        m <- n - 3 - length(sela)
        er <- .lm.fit( z, cbind(y, x[, ind[s]]) )$residuals
        r <- cor(er[, 1], er[, -1])
        n.tests <- n.tests + length( ind[s] )
        stat[ ind[s] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
        critvalue <- qt(quan, m)
        s <- which( stat > critvalue )
        sel <- which.max(stat) * ( length(s) > 0 )
        pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pv)
        sela <- c(sela, sel[sel > 0] )
        s <- s[ - which(s == sel) ]
        z <- cbind(z, x[, sel])
      }  ## end while( sum(s > 0) > 0 )
    }  ##  end if ( length(s) > 0 )
  card <- sum(sela > 0)

  if ( K == 1) {
    stat <- numeric(p)
    m <- n - 3 - length(sela)
    #e1 <- .lm.fit(z, y)$residuals
    #e2 <- .lm.fit(z, x[, ind[-sela] ])$residuals
    er <- .lm.fit( z, cbind(y, x[, ind[-sela]]) )$residuals
    #r <- cor(e2, e1) 
    r <- cor(er[, 1], er[, -1])
    n.tests[2] <- length( ind[-sela] )
    stat[ ind[-sela] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
    critvalue <- qt(quan, m)
    s <- which( stat > critvalue )
    sel <- which.max(stat) * ( length(s) > 0 )
    pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
    sa <- c(sa, stat[sel]) 
    pva <- c(pva, pv)
    sela <- c(sela, sel[sel > 0] )
    s <- s[ - which(s == sel) ]
    z <- cbind(z, x[, sel])     
    while ( sum(s > 0) > 0 ) {
      stat <- numeric(p)
      m <- n - 3 - length(sela)
      #e1 <- .lm.fit(z, y)$residuals
      #e2 <- .lm.fit(z, x[, ind[-sela] ])$residuals
      er <- .lm.fit( z, cbind(y, x[, ind[s]]) )$residuals
      #r <- cor(e2, e1) 
      r <- cor(er[, 1], er[, -1])
      n.tests[2] <- n.tests[2] + length( ind[s] )
      stat[ ind[ s ] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
      critvalue <- qt(quan, m)
      s <- which( stat > critvalue )
      sel <- which.max(stat) * ( length(s) > 0 )
      pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pv)
      sela <- c(sela, sel[sel > 0] )
      s <- s[ - which(s == sel) ]
      z <- cbind(z, x[, sel])     
    } ## end while ( sum(s > 0) > 0 ) 
    card <- c(card, sum(sela > 0) )  
  } ## end  if (K == 1)
 
  if ( K > 1  )  {
    stat <- numeric(p)
    m <- n - 3 - length(sela)
    #e1 <- .lm.fit(z, y)$residuals
    #e2 <- .lm.fit(z, x[, ind[-sela] ])$residuals
    er <- .lm.fit( z, cbind(y, x[, ind[-sela]]) )$residuals
    #r <- cor(e2, e1) 
    r <- cor(er[, 1], er[, -1])
    n.tests[2] <- length( ind[-sela] )
    stat[ ind[-sela] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
    critvalue <- qt(quan, m)
    s <- which( stat > critvalue )
    sel <- which.max(stat) * ( length(s) > 0 )
    pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
    sa <- c(sa, stat[sel]) 
    pva <- c(pva, pv)
    sela <- c(sela, sel[sel > 0] )
    s <- s[ - which(s == sel) ]
    z <- cbind(z, x[, sel])     
    while ( sum(s > 0) > 0 ) {
      stat <- numeric(p)
      m <- n - 3 - length(sela)
      #e1 <- .lm.fit(z, y)$residuals
      #e2 <- .lm.fit(z, x[, ind[-sela] ])$residuals
      er <- .lm.fit( z, cbind(y, x[, ind[s]]) )$residuals
      #r <- cor(e2, e1) 
      r <- cor(er[, 1], er[, -1])
      n.tests[2] <- n.tests[2] + length( ind[s] )
      stat[ ind[s] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
      critvalue <- qt(quan, m)
      s <- which( stat > critvalue )
      sel <- which.max(stat) * ( length(s) > 0 )
      pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
      sa <- c(sa, stat[sel]) 
      pva <- c(pva, pv)
      sela <- c(sela, sel[sel > 0] )
      s <- s[ - which(s == sel) ]
      z <- cbind(z, x[, sel])     
    } ## end while ( sum(s > 0) > 0 ) 
    card <- c(card, sum(sela > 0) )  

    vim <- 1
    while ( vim < K  & card[vim + 1] - card[vim] > 0 ) {
      stat <- numeric(p)
      vim <- vim + 1
      m <- n - 3 - length(sela)
      #e1 <- .lm.fit(z, y)$residuals
      #e2 <- .lm.fit(z, x[, ind[-sela] ])$residuals
      er <- .lm.fit( z, cbind(y, x[, ind[-sela]]) )$residuals
      #r <- cor(e2, e1) 
      r <- cor(er[, 1], er[, -1])
      n.tests[vim + 1] <- length( ind[-sela] )
      stat[ ind[-sela] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
      critvalue <- qt(quan, m)
      s <- which( stat > critvalue )
      sel <- which.max(stat) * ( length(s) > 0 )
      pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
      sa <- c(sa, stat[sel]) 
      s <- s[ - which(s == sel) ]
      pva <- c(pva, pv)
      sela <- c(sela, sel[sel > 0] )
      z <- cbind(z, x[, sel])     
      while ( sum(s > 0) > 0 ) {
        stat <- numeric(p)
        m <- n - 3 - length(sela)
        #e1 <- .lm.fit(z, y)$residuals
        #e2 <- .lm.fit(z, x[, ind[-sela] ])$residuals
        er <- .lm.fit( z, cbind(y, x[, ind[s]]) )$residuals
        #r <- cor(e2, e1) 
        r <- cor(er[, 1], er[, -1])
        n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
        stat[ ind[s] ] <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )    
        critvalue <- qt(quan, m)
        s <- which( stat > critvalue )
        sel <- which.max(stat) * ( length(s) > 0 )
        pv <- log(2) + pt(stat[sel], m, lower.tail = FALSE, log.p = TRUE) 
        sa <- c(sa, stat[sel]) 
        pva <- c(pva, pv)
        sela <- c(sela, sel[sel > 0] )
        s <- s[ - which(s == sel) ]
        z <- cbind(z, x[, sel])     
      } ## end while ( length(s) > 0 )
      card <- c(card, sum(sela > 0) )  
    }  ## end while ( vim < K )
  }  ## end if ( K > 1)
  } ## end if ( length(s) > 0 )
     
  runtime <- proc.time() - runtime
  len <- sum( sela > 0 )
  if (len > 0) {
    res <- cbind(sela[1:len], sa[1:len], pva[1:len] )
    info <- matrix(nrow = length(card), ncol = 2)
    info[, 1] <- card
    info[, 2] <- n.tests
  } else {
    res <- matrix(c(0, 0, 0), ncol = 3)
    info <- matrix(c(0, p), ncol = 2)
  }  
  colnames(res) <- c("Vars", "stat", "log p-value")
  rownames(info) <- paste("K=", 1:length(card)- 1, sep = "")
  colnames(info) <- c("Number of vars", "Number of tests")
  list(runtime = runtime, res = res, info = info)
}


#[export]
cor.fsreg <- function(y, x, ystand = TRUE, xstand = TRUE, threshold = 0.05, tolb = 2, tolr = 0.02, stopping = "BIC") {
  threshold <- log(threshold)
  dm <- dim(x)
  n <- dm[1]
  p <- dm[2]
  con <- n * log(2 * pi) + n  
  logn <- log(n)
  if (xstand)   x <- Rfast::standardise(x)
  if (ystand)   y <- ( y - mean(y) ) / Rfast::Var(y, std = TRUE)
  oop <- options(warn = -1)
  on.exit( options(oop) )
  yx <- Rfast::eachcol.apply(x, y) / (n - 1)
  sel <- which.max( abs(yx) )
  r <- yx[sel] 
  stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(n - 3) )  ## 
  pv <- log(2) + pt(stat, n - 3, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
  model <- NULL
  #############################
  ###### BIC stopping criterion
  #############################
  if ( stopping == "BIC" ) {
    info <- cbind(0, 0, 0)
    tool <- numeric( min(n, p) ) 
    if ( pv < threshold ) {
      info <- cbind(sel, pv, stat)
      z <- x[, sel, drop = FALSE]
      model <- .lm.fit(z, y)
      tool[1] <- n * log( sum(model$residuals^2)/n ) + 3 * logn 
    } else  info <- rbind(info, c(0, 0, 0))  
    if ( !is.null(model) ) {
      xz <- as.vector( cor(z, x) )
	  xz[sel] <- 0
      yx.z <- ( yx - xz * r ) / sqrt(1 - xz^2) / sqrt(1 - r^2)
      sel <- which.max( abs(yx.z) )
      r <- yx.z[sel]
      stat <- 0.5 * abs( log( (1 + r) / (1 - r) ) ) * sqrt(n - 4)    
      pv <- log(2) + pt(stat, n - 4, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
      if ( pv < threshold )  {
        model <- .lm.fit( cbind(z, x[, sel]), y )
        tool[2] <-  n * log( sum(model$residuals^2)/n ) + 4 * logn
        if ( tool[1] - tool[2] > tolb ) {
          info <- rbind(info, c(sel, pv, stat) )
	    z <- cbind(z, x[, info[1:2, 1] ])
          x[, sel] <- 0  
        } else  info <- rbind(info, c(0, 0, 0)) 
      } else  info <- rbind(info, c(0, 0, 0)) 

    }
    k <- 2
  
    if ( info[2, 1] > 0 ) { 
      while ( info[k, 2] < threshold  &  k < n - 20  &  tool[ k - 1 ] - tool[ k ] > tolb  &  k < p )  {
        sela <- info[, 1]
        m <- n - 3 - k 
        k <- k + 1	
        e1 <- model$residuals
        e2 <- .lm.fit(z, x)$residuals
        ## yx.z <- colsums(e1 * e2) / sqrt( Rfast::colsums(e2^2) * sum(e1^2) ) 
        yx.z <- cor(e2, e1) 
        sel <- which.max( abs(yx.z) )  
        r <- yx.z[sel]
        stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )  ## 
        pv <- log(2) + pt(stat, m, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
        if ( pv < threshold )  {
          z <- cbind(z, x[, sel])
          model <- .lm.fit(z, y )
          tool[k] <- n * log( sum(model$residuals^2)/n ) + (k + 2) * logn
          if ( tool[k - 1] - tool[k] > tolb ) {
            info <- rbind(info, c(sel, pv, stat) )
            x[, sel] <- 0  
          } else  info <- rbind(info, c(0, 0, 0)) 
        } else  info <- rbind(info, c(0, 0, 0)) 

      } 
    }
    
    info <- cbind(info, tool[1:k] + con)
    colnames(info)[4] <- "bic"
    info <- info[1:c(k-1), , drop = FALSE]   

  #############################
  ###### Adjusted R-squared stopping criterion
  #############################
  } else if ( stopping == "ar2" ) {
    tool <- numeric( min(n, p) ) 
    info <- cbind(0, 0, 0)
    if ( pv < threshold ) {
      info <- cbind(sel, pv, stat)
      z <- x[, sel, drop = FALSE]
      model <- .lm.fit(z, y)
      r2 <- 1 - cor(y, model$residuals )^2
      tool[1] <- 1 - (1 - r2) * (n - 1 ) / ( n - 2)
    } else  info <- rbind(info, c(0, 0, 0))
    if ( !is.null(model) ) {
	  options(warn = -1)
      xz <- as.vector( cor(z, x) )
      yx.z <- ( yx - xz * r ) / sqrt(1 - xz^2) / sqrt(1 - r2)
      sel <- which.max( abs(yx.z) )
      r <- yx.z[sel]
      stat <- 0.5 * abs( log( (1 + r) / (1 - r) ) ) * sqrt(n - 4)    
      pv <- log(2) + pt(stat, n - 4, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
      if ( pv < threshold )  {
        z <- cbind(z, x[, sel ])
        model <- .lm.fit(z, y ) 
        r2 <- 1 - cor(y, model$residuals )^2
        tool[2] <- 1 - (1 - r2) * (n - 1) / ( n - 3)
          if ( tool[2] - tool[1] > tolr ) {
          info <- rbind(info, c(sel, pv, stat) )
	    z <- cbind(z, x[, info[1:2, 1] ])
          x[, info[1:2, 1] ] <- 0  
        } else  info <- rbind(info, c(0, 0, 0))  
       } else  info <- rbind(info, c(0, 0, 0)) 

    }
    k <- 2
  
    if ( info[2, 1] > 0 ) { 
      while ( info[k, 2] < threshold  &  k < n - 20  &  tool[ k ] - tool[ k - 1 ] > tolr  &  k < p )  {
        sela <- info[, 1]
        m <- n - 3 - k
        k <- k + 1
        e1 <- model$residuals
        e2 <- .lm.fit(z, x)$residuals
        ## yx.z <- colsums(e1 * e2) / sqrt( colsums(e2^2) * sum(e1^2) ) 
        yx.z <- cor(e1, e2) 
        sel <- which.max( abs(yx.z) )  
        r <- yx.z[sel]
        stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )  ## 
        pv <- log(2) + pt(stat, m, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
        if ( pv < threshold )  {
          z <- cbind(z, x[, sel])
          model <- .lm.fit(z, y ) 
          r2 <- 1 - cor(y, model$residuals )^2
          tool[k] <- 1 - (1 - r2) * (n - 1) / ( n - k - 1)
          if ( tool[k] - tool[k - 1] > tolr ) {
            info <- rbind(info, c(sel, pv, stat) )
            x[, sel] <- 0 
          } else  info <- rbind(info, c(0, 0, 0)) 
        } else  info <- rbind(info, c(0, 0, 0)) 

      } 
    }
    info <- cbind(info, tool[1:k])
    colnames(info)[4] <- "adjusted R2"
    info <- info[1:c(k-1), , drop = FALSE]   

  #############################
  ###### BIC and adjusted R-squared stopping criterion
  #############################
  } else if ( stopping == "BICR2" ) {
    info <- cbind(0, 0, 0, 0)
    toolb <- toolr <- numeric( min(n, p) )
    if ( pv < threshold ) {
      info <- cbind(sel, pv, stat)
      z <- x[, sel, drop = FALSE]
      model <- .lm.fit(z, y)   
      r2 <- 1 - cor(y, model$residuals )^2
      toolr[1] <- 1 - (1 - r2) * (n - 1 ) / ( n - 2)
      toolb[1] <- n * log( sum(model$residuals^2)/n ) + con + 3 * logn 
    } else  info <- rbind(info, c(0, 0, 0))
      if ( !is.null(model) ) {
        xz <- as.vector( cor(z, x) )
        yx.z <- ( yx - xz * r ) / sqrt(1 - xz^2) / sqrt(1 - r2)
        sel <- which.max( abs(yx.z) )
        r <- yx.z[sel]
        stat <- 0.5 * abs( log( (1 + r) / (1 - r) ) ) * sqrt(n - 4)    
        pv <- log(2) + pt(stat, n - 4, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
        if ( pv < threshold )  {
          z <- cbind(z, x[, sel])
          model <- .lm.fit(z, y ) 
          r2 <- 1 - cor(y, model$residuals )^2
          toolr[2] <- 1 - (1 - r2) * (n - 1) / ( n - 3)
          toolb[2] <-  n * log( sum( model$residuals^2)/n ) + con + 4 * logn
          if ( toolb[1] - toolb[2] > tolb  &  toolr[2] - toolr[1] > tolr ) {
            info <- rbind(info, c(sel, pv, stat) )
	      z <- cbind(z, x[, info[1:2, 1] ])
            x[, info[1:2, 1] ] <- 0 
          } else  info <- rbind(info, c(0, 0, 0))  
        } else  info <- rbind(info, c(0, 0, 0)) 
    }
    k <- 2
  
    if ( info[2, 1] > 0 ) { 
      while ( info[k, 2] < threshold  &  k < n - 20 &  toolb[ k - 1 ] - toolb[ k ] > tolb  & toolr[ k ] - toolr[ k - 1 ] > tolr  &  k < p )  {
        sela <- info[, 1]
        m <- n - 3 - k
        k <- k + 1
        e1 <- model$residuals
        e2 <- .lm.fit(z, x)$residuals
        ## yx.z <- Rfast::colsums(e1 * e2) / sqrt( Rfast::colsums(e2^2) * sum(e1^2) ) 
        yx.z <- cor(e1, e2) 
        sel <- which.max( abs(yx.z) )  
        r <- yx.z[sel]
        stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) * sqrt(m) )  ## 
        pv <- log(2) + pt(stat, m, lower.tail = FALSE, log.p = TRUE)  ## logged p-values  
        if ( pv < threshold )  {
          z <- cbind(z, x[, sel])
          model <- .lm.fit(z, y ) 
          r2 <- 1 - cor(y, model$residuals )^2
          toolr[k] <- 1 - (1 - r2) * (n - 1) / ( n - k - 1)
          toolb[k] <-  n * log( sum(model$residuals^2)/n ) + con + (k + 2) * logn
          if ( toolb[k - 1] - toolb[k] > tolb  &  toolr[k] - toolr[k - 1] > tolr ) {
            info <- rbind(info, c(sel, pv, stat) )
            x[, sel] <- 0 
          } else  info <- rbind(info, c(0, 0, 0)) 
        } else  info <- rbind(info, c(0, 0, 0)) 

      } 
    }

    info <- cbind(info, toolb[1:k], toolr[1:k])
    colnames(info)[4:5] <- c("bic", "adjusted R2")
    info <- info[1:c(k-1), , drop = FALSE]   

  }

  info
}


#[export]
fs.reg <- function(y,ds,sig = 0.05,tol = 2,type = "logistic") {
	x <- .Call(Rfast_fs_reg,y,ds,sig,tol,type)
	colnames(x) <- c("vars","log.pval","stat","bic")
	x
}


#[export]
omp <- function (y, x, xstand = TRUE, tol = qchisq(0.95, 1) + log(length(y)), type = "logistic") {
    tic <- proc.time()
    dm <- dim(x)
    d <- dm[2]
    n <- dm[1]
    ind <- 1:d
    if (xstand)   x <- Rfast::standardise(x)
    phi <- NULL
    oop <- options(warn = -1)
    on.exit( options(oop) )

  a <- try( 
    if ( type == "logistic" ) {
        p <- sum(y)/n
        rho <- -2 * (n * p * log(p) + (n - n * p) * log(1 - p))
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- Rfast::glm_logistic(x[, sel], y)
        est <- exp(-mod$be[1] - x[, sel] * mod$be[2])
        res <- y - 1/(1 + est)
        rho[2] <- mod$devi
        ind[sel] <- 0
        i <- 2
        while ( (rho[i - 1] - rho[i]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- Rfast::glm_logistic(x[, sela], y)
            est <- as.vector(exp(-mod$be[1] - x[, sela] %*% mod$be[-1]))
            res <- y - 1/(1 + est)
            rho[i] <- mod$devi
            ind[sela] <- 0
        }
    }
    else if ( type == "poisson" ) {
        m <- sum(y)/n
        rho <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * 
            log(m)
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- Rfast::glm_poisson(x[, sel], y)
        res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
        rho[2] <- mod$devi
        ind[sel] <- 0
        i <- 2
        while ( (rho[i - 1] - rho[i]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- Rfast::glm_poisson(x[, sela], y)
            res <- y - as.vector(exp(mod$be[1] + x[, sela] %*% 
                mod$be[-1]))
            rho[i] <- mod$devi
            ind[sela] <- 0
        }
    }
    else if ( type == "quasipoisson" ) {
        m <- sum(y)/n
        rho <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * n * m * 
            log(m)
        phi <- 1
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- try( Rfast::qpois.reg(x[, sel], y), silent = TRUE )
        if ( identical(class(mod), "try-error") ) {
          rho[2] <- rho[1]
          phi[2] <- phi[1]
        } else {
          phi[2] <- mod$phi
          res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
          rho[2] <- mod$devi
          ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i - 1] - rho[i])/phi[i] > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- try( Rfast::qpois.reg(x[, sela], y), silent = TRUE )
            if ( identical(class(mod), "try-error") ) {
              rho[i] <- rho[i - 1]
              phi[i] <- phi[i - 1]
            } else {
              res <- y - as.vector( exp(mod$be[1] + x[, sela] %*% mod$be[-1]) )
              rho[i] <- mod$devi
              phi[i] <- mod$phi
              ind[sela] <- 0
            }
        }
    }
    else if ( type == "quasibinomial" ) {
        p <- sum(y)/n
        y0 <- 1 - y
        rho <- 2 * sum(y * log(y/p), na.rm = TRUE) + 2 * sum(y0 * 
            log(y0/(1 - p)), na.rm = TRUE)
        phi[1] <- 1
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- try( Rfast::prop.reg(y, x[, sel], varb = "glm"), silent = TRUE )
        if ( identical(class(mod), "try-error") ) {
          rho[2] <- rho[1]
          phi[2] <- phi[1]
        } else {
          phi[2] <- mod$phi
          est <- exp(-mod$info[1, 1] - x[, sel] * mod$info[2, 1])
          p <- 1/(1 + est)
          res <- y - p
          rho[2] <- 2 * sum(y * log(y/p), na.rm = TRUE) + 2 * sum(y0 * 
            log(y0/(1 - p)), na.rm = TRUE)
          ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i - 1] - rho[i])/phi[i] > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- try( Rfast::prop.reg(y, x[, sela], varb = "glm"), silent = TRUE )
            if ( identical(class(mod), "try-error") ) {
               rho[i] <- rho[i - 1]
               phi[i] <- phi[i - 1]
            } else {
              est <- as.vector( exp(-mod$info[1, 1] - x[, sela] %*% mod$info[-1, 1]) )
              p <- 1/(1 + est)
              res <- y - p
              rho[i] <- 2 * sum(y * log(y/p), na.rm = TRUE) + 2 * sum(y0 * log(y0/(1 - p)), na.rm = TRUE)
              phi[i] <- mod$phi
              ind[sela] <- 0
            }
        }
    }
    else if ( type == "normlog" ) {
        ini <- Rfast::normlog.mle(y)
        rho <- sum( (y - ini$param[1])^2 )
        phi[1] <- 1
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- try( Rfast::normlog.reg(y, x[, sel]), silent = TRUE )
        if ( identical(class(mod), "try-error") ) {
          rho[2] <- rho[1]
          phi[2] <- phi[1]
        } else {   
          res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
          rho[2] <- mod$deviance
          phi[2] <- mod$devi/(n - 2)
          ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i - 1] - rho[i])/phi[i] > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- try( Rfast::normlog.reg(y, x[, sela]), silent = TRUE )
            if ( identical(class(mod), "try-error") ) {
              rho[i] <- rho[i - 1]
              phi[i] <- phi[i-1]
            } else {
              res <- y - as.vector( exp(mod$be[1] + x[, sela] %*% mod$be[-1]) )
              rho[i] <- mod$deviance
              phi[i] <- mod$deviance/(n - length(mod$be))
              ind[sela] <- 0
            }
        }
    }
    else if ( type == "gamma" ) {
        ini <- Rfast::gammacon(y)
        rho <- ini$deviance
        phi[1] <- 1
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- try( Rfast::gammareg(y, x[, sel]), silent = TRUE )
        if ( identical(class(mod), "try-error") ) {
          rho[2] <- rho[1]
          phi[2] <- phi[1]
        } else {
          res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
          rho[2] <- mod$info[2]
          phi[2] <- mod$info[3]
          ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i - 1] - rho[i])/phi[i] > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- try( Rfast::gammareg(y, x[, sela]), silent = TRUE )
            if ( identical(class(mod), "try-error") ) {
              rho[i] <- rho[i - 1]
              phi[i] <- phi[i - 1]
            } else {
              res <- y - as.vector(exp(mod$be[1] + x[, sela] %*% mod$be[-1]))
              rho[i] <- mod$info[2]
              phi[i] <- mod$info[3]
              ind[sela] <- 0
            }
        }
    }
    else if ( type == "weibull" ) {
        ini <- Rfast::weibull.mle(y)
        rho <- 2 * ini$loglik
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max(abs(ela))
        sela <- sel
        names(sela) <- NULL
        mod <- try( Rfast::weib.reg(y, x[, sel]), silent = TRUE )
        if ( identical(class(mod), "try-error") ) {
          rho[2] <- rho[1]
        } else {
          res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
          rho[2] <- 2 * mod$loglik
          ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i] - rho[i - 1]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- try( Rfast::weib.reg(y, x[, sela]), silent = TRUE )
            if ( identical(class(mod), "try-error") ) {
              rho[i] <- rho[i - 1]
            } else {
              res <- y - as.vector( exp(mod$be[1] + x[, sela] %*% mod$be[-1]) )
              rho[i] <- 2 * mod$loglik
              ind[sela] <- 0
            }
        }
    }
    else if ( type == "mv" ) {
        p <- dim(y)[2]
        con <- n * p * log(2 * pi) + n * p
        mod <- Rfast::mvnorm.mle(y)
        rho <- -2 * mod$loglik
        res <- Rfast::eachrow(y, mod$mu, oper = "-")
        ela <- numeric(d)
        ela <- Rfast::eachcol.apply(x, res[, 1], indices = ind, 
            oper = "*", apply = "sum")^2
        for (i in 2:p) ela <- ela + Rfast::eachcol.apply(x, res[, 
            i], indices = ind, oper = "*", apply = "sum")^2
        sel <- which.max(ela)
        sela <- sel
        names(sela) <- NULL
        res <- .lm.fit(cbind(1, x[, sela]), y)$residuals
        rho[2] <- con + n * log(det(crossprod(res)/(n - 1)))
        ind[sel] <- 0
        i <- 2
        while ( (rho[i] - rho[i - 1]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res[, 1], indices = ind[ind > 
                0], oper = "*", apply = "sum")^2
            for (j in 2:p) r[ind] <- r[ind] + Rfast::eachcol.apply(x, 
                res[, j], indices = ind[ind > 0], oper = "*", 
                apply = "sum")^2
            sel <- which.max(r)
            sela <- c(sela, sel)
            res <- .lm.fit(cbind(1, x[, sela]), y)$residuals
            rho[i] <- con + n * log(det(crossprod(res)/(n - i - 
                1)))
            ind[sela] <- 0
        }
    }
    else if ( type == "multinomial" ) {
	  y1 <- Rfast::design_matrix(y)
	  p <- dim(y1)[2] - 1
	  mod <- Rfast::multinom.mle(y1)
        rho <- mod$loglik
	  y1 <- y1[, -1, drop = FALSE] 
        for (j in 1:p)  y1[, j] <- as.numeric(y1[, j])
        res <- Rfast::eachrow(y1, mod$prob[-1], oper = "-")
        ela <- numeric(d)
        ela <- Rfast::eachcol.apply(x, res[, 1], indices = ind, 
            oper = "*", apply = "sum")^2
        for (i in 2:p) ela <- ela + Rfast::eachcol.apply(x, res[, 
            i], indices = ind, oper = "*", apply = "sum")^2
        sel <- which.max(ela)
        sela <- sel
        names(sela) <- NULL
        mod <- try( Rfast::multinom.reg(y, x[, sela]), silent = TRUE )
        if (identical(class(mod), "try-error")) {
            rho[2] <- rho[1]
        }
        else {
            est <- exp(cbind(1, x[, sel]) %*% mod$be)
            est <- est/(Rfast::rowsums(est) + 1)
            res <- y1 - est
            rho[2] <-  -2 * mod$loglik
            ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i] - rho[i - 1]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res[, 1], indices = ind[ind > 
                0], oper = "*", apply = "sum")^2
            for (j in 2:p) r[ind] <- r[ind] + Rfast::eachcol.apply(x, 
                res[, j], indices = ind[ind > 0], oper = "*", 
                apply = "sum")^2
            sel <- which.max(r)
            sela <- c(sela, sel)
            mod <- try( Rfast::multinom.reg(y, x[, sela]), silent = TRUE )
            if (identical(class(mod), "try-error")) {
                rho[i] <- rho[i - 1]
            }   else {
                rho[i] <-  -2 * mod$loglik
                ind[sela] <- 0
                est <- exp(cbind(1, x[, sela]) %*% mod$be)
                est <- est/(Rfast::rowsums(est) + 1)
                res <- y1 - est
            }
        }
    }
    runtime <- proc.time() - tic
    len <- length(sela)
    info <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(info) <- c("Selected Vars", "Deviance")
    if (!is.null(phi)) 
        phi <- phi[1:len]
    list(runtime = runtime, phi = phi, info = info)
}


#[export]
ompr <- function (y, x, ystand = TRUE, xstand = TRUE, method = "BIC", tol = 2) {

    dm <- dim(x)
    d <- dm[2]
    n <- dm[1]
    ind <- 1:d
    runtime <- proc.time()
    if (ystand)  {
      m <- sum(y)/n
      y <- (y - m)/Rfast::Var(y, std = TRUE)
    }
    if (xstand)   x <- Rfast::standardise(x)
	
    if (method == "sse") {
        rho <- Rfast::Var(y) * (n - 1)
        r <- Rfast::eachcol.apply(x, y)
        epe <- which( is.na(r) )
        ind[epe] <- 0
        sel <- which.max(abs(r))
        sela <- sel
        res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals
        rho[2] <- sum(res^2)
        ind[sel] <- 0
        r[sel] <- 0
        i <- 2
        while ( (rho[i - 1] - rho[i])/(rho[i - 1]) > tol & i < n ) {
            i <- i + 1
            r[sela] <- 0
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            res <- .lm.fit(x[, sela], y)$residuals
            rho[i] <- sum(res^2)
            ind[sela] <- 0
        }
        len <- length(sela)
        info <- cbind(c(0, sela[-len]), rho[1:len])
        colnames(info) <- c("Vars", "|sse|")

    } else if (method == "BIC") {
        con <- n * log(2 * pi) + n
        rho <- n * log(Var(y) * (n - 1)/n) + 2 * log(n)
        r <- Rfast::eachcol.apply(x, y)
        epe <- which( is.na(r) )
        ind[epe] <- 0
        sel <- which.max(abs(r))
        sela <- sel
        res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals
        rho[2] <- n * log(sum(res^2)/n) + 3 * log(n)
        ind[sel] <- 0
        r[sel] <- 0
        i <- 2
        while (rho[i - 1] - rho[i] > tol & i < n) {
            i <- i + 1
            r[sela] <- 0
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            res <- .lm.fit(x[, sela], y)$residuals
            rho[i] <- n * log(sum(res^2)/n) + (i + 1) * log(n)
            ind[sela] <- 0
        }
        len <- length(sela)
        info <- cbind(c(0, sela[-len]), rho[1:len] + con)
        colnames(info) <- c("Vars", "BIC")

    } else if (method == "ar2") {
        down <- Var(y) * (n - 1)
        rho <- 0
        r <- Rfast::eachcol.apply(x, y)
        epe <- which( is.na(r) )
        ind[epe] <- 0
        sel <- which.max(abs(r))
        sela <- sel
        res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals
        r2 <- 1 - sum(res^2)/down
        rho[2] <- 1 - (1 - r2) * (n - 1)/(n - 2)
        ind[sel] <- 0
        r[sel] <- 0
        i <- 2
        while (rho[i] - rho[i - 1] > tol & i < n) {
            i <- i + 1
            r[sela] <- 0
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            res <- .lm.fit(x[, sela], y)$residuals
            r2 <- 1 - sum(res^2)/down
            rho[i] <- 1 - (1 - r2) * (n - 1)/(n - i - 1)
            ind[sela] <- 0
        }
        len <- length(sela)
        info <- cbind(c(0, sela[-len]), rho[1:len])
        colnames(info) <- c("Vars", "adjusted R2") 

    } else if ( method == "pvalue" ) {
        r <- Rfast::eachcol.apply(x, y) / (n - 1)
        epe <- which( is.na(r) )
        ind[epe] <- 0
        sel <- which.max( abs(r) )
        sela <- sel
        stat <- abs( 0.5 * log( (1 + r[sel])/(1 - r[sel]) ) * sqrt(n - 3) )
        pv <- 2 * pt(stat, n - 3, lower.tail = FALSE)
        if ( pv < tol ) {
          ind[sel] <- 0
          r[sel] <- 0
          res <- .lm.fit(x[, sel, drop = FALSE], y)$residuals 
          sse2 <- sum(res^2)
          i <- 1
          while ( pv[i] < tol  &  i < n - 10 ) {
            i <- i + 1
            sse1 <- sse2
            r[sela] <- 0
            r[ind] <- Rfast::eachcol.apply( x, res, indices = ind[ind > 0] )/Rfast::Var(res, std = TRUE)/(n-1)
            sel <- which.max( abs(r) )
            sela <- c(sela, sel)
            res <- .lm.fit(x[, sela], y)$residuals 
            sse2 <- sum(res^2)
            dof <- n - length(sela)
            stat <- ( sse1 - sse2 ) / ( sse2 / dof  )
            pv[i] <- pf(stat, 1, dof, lower.tail = FALSE)
            ind[sela] <- 0
          }
        }
        len <- length(sela)
        info <- cbind(sela[-len], pv[-len])
        colnames(info) <- c("Vars", "p-value")
    }
    runtime <- proc.time() - runtime
    list(runtime = runtime, info = info)
}

