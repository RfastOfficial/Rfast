#[export]
colmeans <- function(x,parallel = FALSE) {
	UseMethod("colmeans")
}

#[export s3]
colmeans.matrix <- function(x,parallel = FALSE) {
	.Call(Rfast_col_means,x,parallel)
}

#[export s3]
colmeans.data.frame <- function(x,parallel = FALSE) {
	.Call(Rfast_col_means,x,parallel)
}

#[export]
rowmeans <- function(x) {
  	as.vector(.Call(Rfast_row_means,x))
}