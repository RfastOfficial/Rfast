#[export]
dirknn.cv <- function(y, x, k = 5:10, type = "C", folds = NULL, nfolds = 10, 
stratified = TRUE, seed = FALSE, parallel = FALSE, pred.ret = FALSE) {
  
  crit <- matrix(nrow = nfolds, ncol = length(k))
  if ( is.null( folds ) )   folds <- makefolds(y, nfolds = nfolds, stratified = stratified, seed = FALSE)
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


makefolds <- function(target, nfolds = 10, stratified = TRUE, seed = FALSE) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if (seed)  set.seed(1234)

  if ( !stratified ) {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    mat <- matrix(sample(length(target)), ncol = nfolds)
    for (i in 1:c(nfolds - 1)) runs[[i]] <- mat[, i]
    names <- prod(dim(mat)) - length(target)
    runs[[nfolds]] <- mat[1:c(nrow(mat) - names), nfolds]
  } else {
    labs <- unique(target)
    run <- list()
    for (i in 1:length(labs)) {
      names <- which(target == labs[i])
      run[[i]] <- sample(names)
    }
    run <- unlist(run)
    for (i in 1:length(target)) {
      k <- i %% nfolds
      if (k == 0)  k <- nfolds
      runs[[k]] <- c(runs[[k]], run[i])
    }
  }
  runs
}
