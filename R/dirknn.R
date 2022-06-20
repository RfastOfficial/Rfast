#[export]
dirknn <- function(xnew, x, y, k, type = "C", parallel = FALSE) {
  .Call(Rfast_dir_knn,t(xnew),t(x),y,k,type,parallel)
}


#[export]
dirknn.cv <- function(y, x, k = 5:10, type = "C", folds = NULL, nfolds = 10, 
stratified = TRUE, seed = NULL, parallel = FALSE, pred.ret = FALSE) {
  
  crit <- matrix(nrow = nfolds, ncol = length(k))
  if ( is.null( folds ) )   folds <- makefolds(y, nfolds = nfolds, stratified = stratified, seed = NULL)
  preds <- list() 
  y <- as.numeric(y)
  for (i in 1:nfolds) {  
    ytest <- y[ folds[[i] ] ]  ##
    xtest <- x[folds[[ i ]], ]  ## xnew
    xtrain <- x[ - folds[[ i ]], ]  ## x
    ytrain <- y[ - folds[[ i ]] ]  ## y
	if (pred.ret) { 
      preds[[ i ]] <- Rfast::dirknn(xtest, xtrain, ytrain, k = k, type = type, parallel = parallel)  
	  be <- preds[[ i ]]- ytest  ## afaireis apo kath sthlh to ytes
	} else  be <- Rfast::dirknn(xtest, xtrain, ytrain, k = k, type = type, parallel = parallel) - ytest 
    if ( type == "C" ) {
	  crit[i, ] <- Rfast::colmeans(be == 0)  ## pososto swstn se tkathe sthlh
	} else  crit[i, ] <- Rfast::colmeans(be^2)  
  }  
  if ( !pred.ret )  preds <- NULL
  list(preds = preds, crit = Rfast::colmeans(crit) )
}



makefolds <- function(ina, nfolds = 10, stratified = TRUE, seed = NULL) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if ( is.null(seed) )  set.seed(seed)

  if ( !stratified ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    ep <- sample( length(ina) )
    nr <- round( length(ina)/nfolds )
    mat <- matrix( ep[1:(nr * nfolds) ], ncol = nfolds )
    mat[ -c( 1:length(ina) ) ] <- NA
    for ( i in 1:nfolds ) runs[[ i ]] <- mat[, i]
    rem <- ep[ - c(1:(nr * nfolds)) ]
    ela <- sample(nfolds, length(rem))
    if ( length(ela) > 0 )  for ( i in 1:length(ela) )  runs[[ ela[i] ]] <- c( runs[[ i ]], rem[ i ] )
  } else {
    labs <- unique(ina)
    run <- list()
    for (i in 1:length(labs)) {
      names <- which( ina == labs[i] )
      run[[i]] <- sample(names)
    }
    run <- unlist(run)
    for ( i in 1:length(ina) ) {
      k <- i %% nfolds
      if ( k == 0 )  k <- nfolds
      runs[[k]] <- c( runs[[ k ]], run[i] )
    }
  }
  for (i in 1:nfolds)  {
    if ( any( is.na(runs[[ i ]]) ) )  runs[[ i ]] <- runs[[ i ]][ !is.na(runs[[ i ]]) ]
  }
  if ( length(runs[[ nfolds ]]) == 0 ) runs[[ nfolds ]] <- NULL
  runs
}