#[export]
rmdp <- function(y, alpha = 0.05, itertime = 100, parallel = FALSE) {
  ## y is the data
  ## alpha is the significance level
  ## itertime is the number of iterations for the first step of the algorithm
  dm <- dim(y)
  n <- dm[1]  ## sample size
  p <- dm[2]  ## dimensionality
  h <- round(n/2) + 1  ## subset of data to be used for the location 
  ## and scatter estimators
  init_h <- 2  ## initial sample size for the first step of the algorithm
  delta <- alpha/2
  runtime <- proc.time()
  #######################
  ####################### 
  id <- replicate( itertime, sample.int(n, 2) ) - 1
  final_vec <- as.vector(.Call(Rfast_rmdp,y,h,id,itertime, parallel))
  #######################
  ####################### 
  submcd <- seq(1, n)[final_vec != 0]
  mu_t <- Rfast::colmeans( y[submcd, ] ) 
  var_t <- Rfast::colVars( y[submcd, ] )
  sama <- ( t(y) - mu_t )^2 / var_t
  disa <- Rfast::colsums(sama) 
  disa <- disa * p / Rfast::Median(disa) 
  b <- Rfast::hd.eigen(y[submcd, ], center = TRUE, scale = TRUE)$values
  tr2_h <- sum(b^2)
  tr2 <- tr2_h - p^2 / h
  cpn_0 <- 1 + (tr2_h) / p^1.5
  w0 <- (disa - p) / sqrt( 2 * tr2 * cpn_0 ) < qnorm(1 - delta)
  nw <- sum(w0)
  sub <- seq(1, n)[w0]
  mu_t <- Rfast::colmeans( y[sub, ] ) 
  var_t <- Rfast::colVars( y[sub, ])
  sama <- ( t(y) - mu_t )^2 / var_t
  disa <- Rfast::colsums(sama)
  b <- Rfast::hd.eigen(y[sub, ], center = TRUE, scale = TRUE)$values
  tr2_h <- sum(b^2) 
  tr2 <- tr2_h - p^2 / nw
  scal <- 1 + exp( - qnorm(1 - delta)^2 / 2 ) / (1 - delta) * sqrt( tr2) / p / sqrt(pi)
  disa <- disa / scal
  cpn_1 <- 1 + (tr2_h) / p^1.5
  dis <- (disa - p) / sqrt(2 * tr2 * cpn_1 ) 
  wei <- dis < qnorm(1 - alpha)
  ####
  runtime <- proc.time() - runtime
  list(runtime = runtime, dis = dis, wei = wei)
}