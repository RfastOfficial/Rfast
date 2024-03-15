#[export]
rowMads <- function(x,method = "median", na.rm = FALSE,parallel = FALSE, cores = 0) {
	.Call(Rfast_row_mads,x,method,na.rm,parallel,cores)
}

#[export]
colMads <- function(x,method = "median", na.rm = FALSE,parallel = FALSE, cores = 0) {
	.Call(Rfast_col_mads,x,method,na.rm,parallel,cores)
}

#[export]
Mad <- function(x,method = "median",na.rm = FALSE) {
	.Call(Rfast_mad2,x,method,na.rm)
}