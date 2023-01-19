#[export]
knn <- function(xnew, y, x, k, dist.type = "euclidean", type = "C", method = "average", freq.option = 0, mem.eff = FALSE) {
  .Call(Rfast_k_nn,xnew, y, x, k, dist.type, type, method, freq.option, mem.eff)
}



#[export]
knn.cv <- function(folds = NULL, nfolds = 10, stratified = FALSE, seed = NULL, y, x, k, dist.type = "euclidean", type = "C", 
				   method = "average", freq.option = 0, pred.ret = FALSE, mem.eff = FALSE) {
  if (is.null(folds)) folds <- generateFolds(y, nfolds = nfolds, stratified = stratified, seed = seed) 
  .Call(Rfast_k_nn_cv,folds, y, x, k, dist.type, type, method, freq.option, pred.ret, mem.eff)
}

generateFolds <- function (ina, nfolds = 10, stratified = TRUE, seed = NULL) {
    names <- paste("Fold", 1:nfolds)
    runs <- sapply(names, function(x) NULL)
    if (!is.null(seed)) 
        set.seed(seed)
    oop <- options(warn = -1)
    on.exit(options(oop))
    if (!stratified) {
        rat <- length(ina)%%nfolds
        mat <- matrix(Rfast2::Sample.int(length(ina)), ncol = nfolds)
        mat[-c(1:length(ina))] <- NA
        for (i in 1:c(nfolds - 1)) runs[[i]] <- mat[, i]
        a <- prod(dim(mat)) - length(ina)
        runs[[nfolds]] <- mat[1:c(nrow(mat) - a), nfolds]
    }
    else {
        labs <- unique(ina)
        run <- list()
        for (i in 1:length(labs)) {
            names <- which(ina == labs[i])
            run[[i]] <- sample(names)
        }
        run <- unlist(run)
        for (i in 1:length(ina)) {
            k <- i%%nfolds
            if (k == 0) 
                k <- nfolds
            runs[[k]] <- c(runs[[k]], run[i])
        }
    }
    for (i in 1:nfolds) {
        if (any(is.na(runs[[i]]))) 
            runs[[i]] <- runs[[i]][!is.na(runs[[i]])]
    }
    if (length(runs[[nfolds]]) == 0) 
        runs[[nfolds]] <- NULL
    runs
}


