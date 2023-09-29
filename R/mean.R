#[export]
colmeans <- function(x,parallel = FALSE) {
	UseMethod("colmeans")
}

#[export s3]
colmeans.matrix <- function(x,parallel = FALSE,cores = 0) {
	.Call(Rfast_col_means,x,parallel,cores)
}

#[export s3]
colmeans.data.frame <- function(x,parallel = FALSE,cores = 0) {
	.Call(Rfast_col_means,x,parallel,cores)
}

#[export]
rowmeans <- function(x,cores = 0) {
  	as.vector(.Call(Rfast_row_means,x,cores))
}