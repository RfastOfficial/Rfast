#[export]
dirknn <- function(xnew, x, y, k, type = "C", parallel = FALSE) {
  .Call(Rfast_dir_knn,t(xnew),t(x),y,k,type,parallel)
}


#[export]
dirknn.cv <- function(y, x, k = 5:10, type = "C", folds = NULL, nfolds = 10, 
stratified = TRUE, seed = NULL, parallel = FALSE, pred.ret = FALSE) {
  
  crit <- matrix(nrow = nfolds, ncol = length(k))
  if ( is.null( folds ) )   folds <- .makefolds(y, nfolds = nfolds, stratified = stratified, seed = NULL)
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


