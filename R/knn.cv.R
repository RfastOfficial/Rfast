
#[export]
knn.cv <- function(folds = NULL, nfolds = 10, stratified = FALSE, seed = FALSE, y, x, k, dist.type = "euclidean", type = "C", 
				   method = "average", freq.option = 0, pred.ret = FALSE, mem.eff = FALSE) {
  if (is.null(folds)) folds <- generateFolds(y, nfolds = nfolds, stratified = stratified, seed = seed) 
  .Call(Rfast_k_nn_cv,folds, y, x, k, dist.type, type, method, freq.option, pred.ret, mem.eff)
}

generateFolds <- function(target, nfolds = 10, stratified = T, seed = F) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if (seed) set.seed(1234)
  
  if (!stratified) {
    options(warn = -1)
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