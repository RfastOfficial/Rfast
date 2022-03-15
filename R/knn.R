#[export]
knn <- function(xnew, y, x, k, dist.type = "euclidean", type = "C", method = "average", freq.option = 0, mem.eff = FALSE) {
  .Call(Rfast_k_nn,xnew, y, x, k, dist.type, type, method, freq.option, mem.eff)
}



#[export]
knn.cv <- function(folds = NULL, nfolds = 10, stratified = FALSE, seed = FALSE, y, x, k, dist.type = "euclidean", type = "C", 
				   method = "average", freq.option = 0, pred.ret = FALSE, mem.eff = FALSE) {
  if (is.null(folds)) folds <- generateFolds(y, nfolds = nfolds, stratified = stratified, seed = seed) 
  .Call(Rfast_k_nn_cv,folds, y, x, k, dist.type, type, method, freq.option, pred.ret, mem.eff)
}

generateFolds <- function(target, nfolds = 10, stratified = TRUE, seed = NULL) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if ( is.null(seed) )  set.seed(seed)

  if (!stratified) {
    oop <- options(warn = -1)
    on.exit(options(oop))
    ep <- sample(length(target))
    nr <- round(length(target) / nfolds)
    mat <- matrix(ep[1:(nr * nfolds)], ncol = nfolds)
    for (i in 1:nfolds) runs[[i]] <- mat[, i]
    rem <- ep[-c(1:(nr * nfolds))]
    ela <- sample(nfolds, length(rem))
    if (length(ela) > 0) {
        for (i in 1:length(ela)) {
            runs[[ela[i]]] <- c(runs[[i]], rem[i] )
        }
    }
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
  for (i in 1:nfolds)  {
    if (any(is.na(runs[[i]]))) {
        runs[[i]] <- runs[[i]][!is.na(runs[[i]])]
    }
  }
  runs
}

