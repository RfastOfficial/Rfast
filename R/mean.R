#[export]
colmeans <- function(x,parallel = FALSE) {
	.Call(Rfast_col_means,x,parallel)
}

#[export]
rowmeans <- function(x) {
  	as.vector(.Call(Rfast_row_means,x))
}