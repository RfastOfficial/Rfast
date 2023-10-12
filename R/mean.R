#[export]
colmeans <- function(x,parallel = FALSE,cores = 0) {
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
rowmeans <- function(x) {
  	as.vector(.Call(Rfast_row_means,x))
}